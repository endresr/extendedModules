'''
Code for finding AR-quiver of extended module categories of linear Nakayama algebras

Parts known to depend on the algebra being linear Nakayama:
- Radical check only works if the radical is indecomposable
'''
from classes import *
from functions import *


import numpy as np


n=5
m=3

relations = [(0,2),(2,4)]#list of minimal relations 

cutOffIterations = 40 #How many times do the while loop run before we give up?

#Note that quivers are zero-indexed, i.e. Q_n: 0 -> 1 -> ... -> (n-1)


n=5
m=3
rel = [(0,2),(2,4)]

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
    for N in [Y for Y in X.irrFrom if Y.injective == None]:
        X.irrTo.append(N.tauInv)
        N.tauInv.irrFrom.append(X)
    if X.injectiveShift[0] == None:
        tempHomologies = -X.homologies
        for N in X.irrTo:
            tempHomologies = tempHomologies+N.homologies
        tauInverse = extendedModule(tempHomologies,n,m,dimProjectives,dimInjectives,radicalOfProjectives)
        X.tauInv = tauInverse
        modules.append(tauInverse)
    elif X.injective == None:
        tempHomologies = np.append(np.zeros((X.injectiveShift[0]+1,n)),[dimProjectives[X.injectiveShift[1]]],axis=0)
        tempHomologies = np.append(tempHomologies,np.zeros((m-X.injectiveShift[0]-2,n)),axis=0)
        tauInverse = extendedModule(tempHomologies,n,m,dimProjectives,dimInjectives,radicalOfProjectives)
        X.tauInv = tauInverse
        modules.append(tauInverse)
    return 



Next = []
tempVect = np.array([1 for i in range(n)])
for i in range(len(dimProjectives)):
    if np.dot(dimProjectives[i],tempVect)==1:
        Next.append(projectiveModules[i])

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

for X in projectiveModules:
    orbit = [X]
    temp=X
    counter=1
    while temp.injective == None and counter < cutOffIterations:
        print(temp)
        temp=temp.tauInv
        print(temp)
        orbit.append(temp)
        counter +=1
    tauOrbits.append(orbit)

