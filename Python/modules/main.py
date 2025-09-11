'''
Code for finding AR-quiver of extended module categories of linear Nakayama algebras

 0 -> 1 -> -> 2 -> ... -> n-1

Prints Latex code for the AR-quiver and compiles it into the folder outputLatexFiles. 
'''
from classes import *
from functions import *
from drawGraph import *
from tauFunction import tauInv
from datetime import datetime


import numpy as np

'''
Initial variables
'''

n=8 #Number of vertices
m=4 #how extended the module category is
rel = 3 #[(0,9),(3,12)] # l-integer if homogeneous relations Rad^l, or a list of minimal zero-relations given through the vertices they start and end in 
cutOffIterations = 200 #How many times do the while loop run before we give up?

#Note that quivers are zero-indexed, i.e. Q_n: 0 -> 1 -> ... -> (n-1)

yLevels = [-3,-2,-1,4,3,2,1,0] #Set the y-level which the tauOrbits of each projective is drawn
tikzScale = (4,2) #The x- and y-scale of the tikz diagram, Recommended: (4,2) works well for most instances when printing homology
nodeScale = 1 #The scale of each node in the tikz diagram
setOutputName = None #String with your prefered name for Latex-file
generateLatex = True #Set to False if you do not want to generate Latex-file
compileToPDF = True #Set to False if you do not want to automatically compile pdf
highlightConcentradedHomology = False #Set to False if you do not want to highlight nodes with concentrated homology
drawOnlyCircles = False #Set to False if you want the homologies printed as node labels
printDiagramTitle = False

'''
The main loop
'''

def mainLoop(n,rel,m,cutOff=100,yLevels=None,tikzScale=(1,1),nodeScale=1,outputName=None,outputLatex=True,compileLatex=True,hightlightConcHom=False,drawCircles=False,printLabel=True):
    SetOutputName = outputName
    if isinstance(rel, int):
        relationsInput = [(i,i+rel) for i in range(n-rel)]
        if outputName == None:
            SetOutputName = "HomogeneousNakayama_n="+str(n)+"_m="+str(m)+"_l="+str(rel)
    else:
        relationsInput = rel
        
    relations = checkRelations(n,relationsInput)
    if relations == None:
        print(f"\033[1;31m The relations are not in accepted form \033[0m")
        return None
    if yLevels==None or len(yLevels)<n:
        print("Using standard ylevels")
        yLevelsTauOrbits=range(n)[::-1]
    else:
        yLevelsTauOrbits=yLevels
    
    #Initial calculations
    

    dimProjectives,radicalOfProjectives = findProjectivesDim(n,relations) 
    dimInjectives = findInjectiveDim(n,relations)[::-1]

    homologiesOfProjectives = findHomologyOfProjectives(dimProjectives,n,m)
    projectiveModules = [extendedModule(x,n,m,dimProjectives,dimInjectives,radicalOfProjectives) for x in homologiesOfProjectives]
    modules = projectiveModules[:]

    
    #The loop constructing the AR-Quiver
    
    Next = []
    tempVect = np.array([1 for i in range(n)])
    for i in range(len(dimProjectives)):
        if np.dot(dimProjectives[i],tempVect)==1:
            Next.append(projectiveModules[i]) #Setting our beginning to the simple projectives
            projectiveModules[i].xcoord = 0 #Setting the simple projectives x-coordinate to 0 
    
    counter=1
    notFinishedARquiver = False
    continueLoop = True
    while len(Next)>0 and counter<cutOff and continueLoop:
        tempNext = []
        for M in Next:
            if (M.homologies<0).any():
                continueLoop=False
            if M.tauInv == None:
                tauInv(M,projectiveModules,radicalOfProjectives,dimProjectives,dimInjectives,n,m,modules)
            tempNext+=(M.irrTo)
        Next=list( dict.fromkeys(tempNext) )
        counter+=1
        if counter>=cutOff and len(Next)>0:
            print(f"\033[1;31m Reached cutoff-value before finishing.\033[0m")
            notFinishedARquiver = True
        if continueLoop == False and len(Next)>0:
            print("Error: Homologies are non-positive")
    calculatedModules=len(modules)

    
    #Generate tikz-code
    
    tauOrbits = []

    for i in range(len(projectiveModules)):
        temp=projectiveModules[i]
        orbit = [temp]
        counter=0
        temp.id=f"t-{counter}P{i}"
        while temp.injective == None and counter < cutOff and temp.tauInv !=None:
            counter +=1
            temp=temp.tauInv
            temp.id=f"t-{counter}P{i}"
            orbit.append(temp)
        tauOrbits.append(orbit)
    TikzNodes = ""
    TikzIrrArrows = ""
    TikzTauArrows = ""

    xMax=0
    for i in range(n):
        y=yLevelsTauOrbits[i]
        for j in range(len(tauOrbits[i])):
            M=tauOrbits[i][j]
            TikzNodes += drawNodes(M,M.xcoord,y,nodeScale,hightlightConcHom,drawCircles)
            TikzIrrArrows += drawIrrArrows(M)
            if M.xcoord!=None and M.xcoord>xMax:
                xMax=M.xcoord
            if M.tauInv != None:
                TikzTauArrows += drawTauArrow(M)

    xMidway=xMax/2
    yMax = max(yLevelsTauOrbits)+1
    TikzLabel = r'\node at (' + str(xMidway)+ ','+ str(yMax)+') [] '+r'{$'+str(m)+"$-mod of linear Nakayama with "+str(n)+" vertices and relations "+str(relations) +r'};'+'\n'
    if printLabel:
        stringToSave = preLatex + preTikz(tikzScale)+TikzLabel +TikzNodes+TikzIrrArrows+TikzTauArrows+postTikz+postLatex
    else:
        stringToSave = preLatex + preTikz(tikzScale)+TikzNodes+TikzIrrArrows+TikzTauArrows+postTikz+postLatex
    currentTime = datetime.today().strftime('%Y-%m-%d')
    if SetOutputName == None:
        texFileName = "_Nakayama_n="+str(n)+"_m="+str(m)
    else:
        texFileName = SetOutputName
    outputDirectory = "./Python/modules/outputLatexFiles/"
    fileName =outputDirectory+currentTime+texFileName+".tex"
    if outputLatex:
        with open(fileName,"w") as file:
            file.write(stringToSave)
        if compileLatex:
            compile_latex(fileName,outputDirectory)
    if notFinishedARquiver:
        print(f"\033[1;31m ERROR: The auslander reiten sequence was not fully calculated.\033[0m")
    print(f"\033[1;32m We found "+str(calculatedModules)+" "+str(m)+"-modules of the algebra.\033[m")
    return tauOrbits,modules



tauOrbits,modules=mainLoop(n,rel,m,cutOffIterations,yLevels,tikzScale,nodeScale,setOutputName,generateLatex,compileToPDF,highlightConcentradedHomology,drawOnlyCircles,printDiagramTitle)
