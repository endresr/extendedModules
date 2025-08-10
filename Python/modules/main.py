'''
Code for finding AR-quiver of extended module categories of linear Nakayama algebras

 0 -> 1 -> -> 2 -> ... -> n-1

Prints Latex code for the AR-quiver. The nodes are the homologies of the modules, represented as matrices with the ith homology in the -ith row
 
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

n=9 #Number of vertices

m=3 #how extended the module category is

rel = [(0,4),(1,5),(2,6),(3,7),(4,8)] #list of minimal zero-relations, given through the vertices they start and end in
cutOffIterations = 10000 #How many times do the while loop run before we give up?

#Note that quivers are zero-indexed, i.e. Q_n: 0 -> 1 -> ... -> (n-1)

yLevelsTauOrbtis = [-1,6,-1,5,3,4,2,1,0] #Set the y-level which the tauOrbits of each projective is drawn
tikzScale = (1,2) #The x- and y-scale of the tikz diagram
nodeScale = 0.5 #The scale of each node in the tikz diagram


'''
Initial calcualations
'''

dimProjectives,radicalOfProjectives = findProjectivesDim(n,rel) 
dimInjectives = findInjectiveDim(n,rel)[::-1]

homologiesOfProjectives = findHomologyOfProjectives(dimProjectives,n,m)
projectiveModules = [extendedModule(x,n,m,dimProjectives,dimInjectives,radicalOfProjectives) for x in homologiesOfProjectives]
modules = projectiveModules[:]

'''
Defining the tau-function
'''

def tauInv(X,projectives,radicalOfProjectives,dimProjectives,dimInjectives,n,m):
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
The loop constructing the AR-quiver
'''

Next = []
tempVect = np.array([1 for i in range(n)])
for i in range(len(dimProjectives)):
    if np.dot(dimProjectives[i],tempVect)==1:
        Next.append(projectiveModules[i])
        projectiveModules[i].xcoord = 0

counter=1
while len(Next)>0 and counter<cutOffIterations:
    tempNext = []
    for M in Next:
        if M.tauInv == None:
            tauInv(M,projectiveModules,radicalOfProjectives,dimProjectives,dimInjectives,n,m)
        tempNext+=(M.irrTo)
    Next=list( dict.fromkeys(tempNext) )
    counter+=1
    if counter>=cutOffIterations and len(Next)>0:
        print("Reached cutoff-value before finishing.")

print(str(len(modules))+" calculated modules")
calculatedModules=len(modules)

#Generate tikz-code
tauOrbits = []

for i in range(len(projectiveModules)):
    temp=projectiveModules[i]
    orbit = [temp]
    counter=0
    temp.id=f"t-{counter}P{i}"
    while temp.injective == None and counter < calculatedModules:
        counter +=1
        temp=temp.tauInv
        temp.id=f"t-{counter}P{i}"
        orbit.append(temp)
    tauOrbits.append(orbit)




TikzNodes = ""
TikzIrrArrows = ""
TikzTauArrows = ""

for i in range(n):
    y=yLevelsTauOrbtis[i]
    for j in range(len(tauOrbits[i])):
        M=tauOrbits[i][j]
        TikzNodes += drawNodes(M,M.xcoord,y,nodeScale)
        TikzIrrArrows += drawIrrArrows(M)
        if M.tauInv != None:
            TikzTauArrows += drawTauArrow(M)

stringToSave = preLatex + preTikz(tikzScale) +TikzNodes+TikzIrrArrows+TikzTauArrows+postTikz+postLatex

currentTime = datetime.today().strftime('%Y-%m-%d')
#currentTime = datetime.today().strftime('%Y-%m-%d_%H%M') #If you want to have unique save each minute
fileName ="./Python/modules/outputLatexFiles/"+currentTime+"_Nakayama_n="+str(n)+"_m="+str(m)+".tex"
print(fileName)
with open(fileName,"w") as file:
    file.write(stringToSave)

outputDirectory = "./Python/modules/outputLatexFiles/"
compile_latex(fileName,outputDirectory)
