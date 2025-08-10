'''
Code for finding AR-quiver of extended module categories of linear Nakayama algebras

 0 -> 1 -> -> 2 -> ... -> n-1

Prints Latex code for the AR-quiver and compiles it into the folder outputLatexFiles. The nodes are the homologies of the modules, represented as matrices with the ith homology in the -ith row
 
Limitations:
- Radical check only works if the radical is indecomposable
- Finding projectives and injectives only works for linear nakayama
'''
from classes import *
from functions import *
from drawGraph import *
from datetime import datetime


import numpy as np

'''
Initial variables
'''

n=5 #Number of vertices
m=3 #how extended the module category is

rel = [(0,2),(2,4)] #list of minimal zero-relations, given through the vertices they start and end in
cutOffIterations = 100 #How many times do the while loop run before we give up?

#Note that quivers are zero-indexed, i.e. Q_n: 0 -> 1 -> ... -> (n-1)

yLevels = [1,2,3,4,5] #Set the y-level which the tauOrbits of each projective is drawn
tikzScale = (2,2) #The x- and y-scale of the tikz diagram
nodeScale = 0.5 #The scale of each node in the tikz diagram
setOutputName = None #String with your prefered name for Latex-file
generateLatex = True #Set to False if you do not want to generate Latex-file
compileToPDF = True #Set to False if you do not want to automatically compile pdf

'''
Defining the tau-function
'''

def tauInv(X,projectives,radicalOfProjectives,dimProjectives,dimInjectives,n,m,modules):
    if X.radicalOfProjective != None:
        #Note the following only works when the radicals are indecomposable
        X.irrTo.append(projectives[X.radicalOfProjective])
        projectives[X.radicalOfProjective].irrFrom.append(X)
        projectives[X.radicalOfProjective].xcoord = X.xcoord+1
    for N in [Y for Y in X.irrFrom if Y.injective == None]:
        X.irrTo.append(N.tauInv)
        N.tauInv.irrFrom.append(X)
    if X.injectiveShift[0] == None:
        tempHomologies = -X.homologies
        for N in X.irrTo:
            tempHomologies = tempHomologies+N.homologies
        tauInverse = extendedModule(tempHomologies,n,m,dimProjectives,dimInjectives,radicalOfProjectives)
        tauInverse.xcoord = X.xcoord+2
        X.tauInv = tauInverse
        modules.append(tauInverse)
    elif X.injective == None:
        tempHomologies = np.append(np.zeros((X.injectiveShift[0]+1,n)),[dimProjectives[X.injectiveShift[1]]],axis=0)
        tempHomologies = np.append(tempHomologies,np.zeros((m-X.injectiveShift[0]-2,n)),axis=0)
        tauInverse = extendedModule(tempHomologies,n,m,dimProjectives,dimInjectives,radicalOfProjectives)
        tauInverse.xcoord = X.xcoord+2
        X.tauInv = tauInverse
        modules.append(tauInverse)
    return 

'''
The main loop
'''

def mainLoop(n,rel,m,cutOff=100,yLevels=None,tikzScale=(1,1),nodeScale=1,outputName=None,outputLatex=True,compileLatex=True):
    relations = checkRelations(n,rel)
    if relations == None:
        print("The relations are not in accepted form")
        return None
    if yLevels==None or len(yLevels)<n:
        print("Using standard ylevels")
        yLevelsTauOrbits=yLevels
    else:
        yLevelsTauOrbits=yLevels
    
    #Initial calculations
    

    dimProjectives,radicalOfProjectives = findProjectivesDim(n,rel) 
    dimInjectives = findInjectiveDim(n,rel)[::-1]

    homologiesOfProjectives = findHomologyOfProjectives(dimProjectives,n,m)
    projectiveModules = [extendedModule(x,n,m,dimProjectives,dimInjectives,radicalOfProjectives) for x in homologiesOfProjectives]
    modules = projectiveModules[:]

    
    #The loop constructing the AR-Quiver
    
    Next = []
    tempVect = np.array([1 for i in range(n)])
    for i in range(len(dimProjectives)):
        if np.dot(dimProjectives[i],tempVect)==1:
            Next.append(projectiveModules[i]) #Setting our beginning to the simple projectives
            projectiveModules[i].xcoord = 0 #Setting the simple projectives x-coordinate to 0 TODO: Find how to manage for algebras with more than one simple projective
    
    counter=1
    notFinishedARquiver = False
    while len(Next)>0 and counter<cutOff:
        tempNext = []
        for M in Next:
            if M.tauInv == None:
                tauInv(M,projectiveModules,radicalOfProjectives,dimProjectives,dimInjectives,n,m,modules)
            tempNext+=(M.irrTo)
        Next=list( dict.fromkeys(tempNext) )
        counter+=1
        if counter>=cutOff and len(Next)>0:
            print("Reached cutoff-value before finishing.")
            notFinishedARquiver = True
    print(str(len(modules))+" calculated modules")
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

    for i in range(n):
        y=yLevelsTauOrbits[i]
        for j in range(len(tauOrbits[i])):
            M=tauOrbits[i][j]
            TikzNodes += drawNodes(M,M.xcoord,y,nodeScale)
            TikzIrrArrows += drawIrrArrows(M)
            if M.tauInv != None:
                TikzTauArrows += drawTauArrow(M)

    stringToSave = preLatex + preTikz(tikzScale) +TikzNodes+TikzIrrArrows+TikzTauArrows+postTikz+postLatex
    currentTime = datetime.today().strftime('%Y-%m-%d')
    if outputName == None:
        texFileName = "_Nakayama_n="+str(n)+"_m="+str(m)
    else:
        texFileName = outputName
    outputDirectory = "./Python/modules/outputLatexFiles/"
    fileName =outputDirectory+currentTime+texFileName+".tex"
    if outputLatex:
        with open(fileName,"w") as file:
            file.write(stringToSave)
        if compileLatex:
            compile_latex(fileName,outputDirectory)
    if notFinishedARquiver:
        print("ERROR: The auslander reiten sequence was not fully calculated.")
    print("We found "+str(calculatedModules)+" "+str(m)+"-modules of the algebra.")
    return tauOrbits

mainLoop(n,rel,m,cutOffIterations,yLevels,tikzScale,nodeScale,setOutputName,generateLatex,compileToPDF)
