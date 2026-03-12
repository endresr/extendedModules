from functions import *
from classes import *

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