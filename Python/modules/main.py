'''
Code for finding AR-quiver of extended module categories of linear Nakayama algebras

Parts known to depend on the algebra being linear Nakayama:
- Radical check only works if the radical is indecomposable
'''
from classes import *
from functions import *
from drawGraph import *
from datetime import datetime

import numpy as np


n=5
m=10

rel = [(0,2),(2,4)]
cutOffIterations = 100 #How many times do the while loop run before we give up?

#Note that quivers are zero-indexed, i.e. Q_n: 0 -> 1 -> ... -> (n-1)




dimProjectives,radicalOfProjectives = findProjectivesDim(n,rel)
dimInjectives = findInjectiveDim(n,rel)[::-1]

homologiesOfProjectives = findHomologyOfProjectives(dimProjectives,n,m)
projectiveModules = [extendedModule(x,n,m,dimProjectives,dimInjectives,radicalOfProjectives) for x in homologiesOfProjectives]
modules = projectiveModules[:]

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


#Generate tikz-code
tauOrbits = []

for i in range(len(projectiveModules)):
    temp=projectiveModules[i]
    orbit = [temp]
    counter=0
    temp.id=f"t-{counter}P{i}"
    while temp.injective == None and counter < cutOffIterations:
        counter +=1
        temp=temp.tauInv
        temp.id=f"t-{counter}P{i}"
        orbit.append(temp)
    tauOrbits.append(orbit)

print(drawNodes(tauOrbits[4][5],1,2)+drawNodes(tauOrbits[4][6],1,2))
print(drawIrrArrows(tauOrbits[4][5]))
print(drawTauArrow(tauOrbits[4][6]))


TikzNodes = ""
TikzIrrArrows = ""
TikzTauArrows = ""

'''
for i in range(len(tauOrbits[5])):
    M=tauOrbits[5][i]
    TikzNodes += drawNodes(M,M.xcoord,3)
    TikzIrrArrows += drawIrrArrows(M)
    if M.tauInv != None:
        TikzTauArrows += drawTauArrow(M)
'''

for i in range(len(tauOrbits[4])):
    M=tauOrbits[4][i]
    TikzNodes += drawNodes(M,M.xcoord,0)
    TikzIrrArrows += drawIrrArrows(M)
    if M.tauInv != None:
        TikzTauArrows += drawTauArrow(M)

for i in range(len(tauOrbits[3])):
    M=tauOrbits[3][i]
    TikzNodes += drawNodes(M,M.xcoord,1)
    TikzIrrArrows += drawIrrArrows(M)
    if M.tauInv != None:
        TikzTauArrows += drawTauArrow(M)

for i in range(len(tauOrbits[2])):
    M=tauOrbits[2][i]
    TikzNodes += drawNodes(M,M.xcoord,-1)
    TikzIrrArrows += drawIrrArrows(M)
    if M.tauInv != None:
        TikzTauArrows += drawTauArrow(M)

for i in range(len(tauOrbits[1])):
    M=tauOrbits[1][i]
    TikzNodes += drawNodes(M,M.xcoord,-2)
    TikzIrrArrows += drawIrrArrows(M)
    if M.tauInv != None:
        TikzTauArrows += drawTauArrow(M)

for i in range(len(tauOrbits[0])):
    M=tauOrbits[0][i]
    TikzNodes += drawNodes(M,M.xcoord,2)
    TikzIrrArrows += drawIrrArrows(M)
    if M.tauInv != None:
        TikzTauArrows += drawTauArrow(M)

stringToSave = preLatex + preTikz +TikzNodes+TikzIrrArrows+TikzTauArrows+postTikz+postLatex

currentTime = datetime.today().strftime('%Y-%m-%d')
#currentTime = datetime.today().strftime('%Y-%m-%d_%H%M') #If you want to have unique save each minute
fileName ="./Python/modules/outputLatexFiles/"+currentTime+"_Nakayama_n="+str(n)+"_m="+str(m)+".tex"
print(fileName)
with open(fileName,"w") as file:
    file.write(stringToSave)