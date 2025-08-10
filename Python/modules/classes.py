'''
Construct an extended module
'''

from functions import *
import numpy as np

class extendedModule:
    '''
    Constructing the modules

    Input:
        - homologies : an m x n numpy array where the ith row is the dimension vector of the ith homology
        - n : number of vertices in the underlying quiver
        - m : how extended are the modules
        - projectivesDimVect : the dimension vectors of the projectives in mod(Lambda)
        - injectivesDimVect : the dimension vectors of the injectives in mod(Lambda)
        - projectiveRadicals : an n x n numpy array where the ith row is the dimension vector of Rad(P_i) in mod(Lambda)

    Output: An object representing the extended module
    '''
    def __init__(self,homologies,n,m,projectivesDimVect,injectivesDimVect,projectiveRadicals):
        self.homologies = homologies
        self.id = None #Is added when we construct the Latex-code
        self.n = n
        self.xcoord = None #Is added dependent on the x-coordinate of the modules in self.irrTo
        self.concentratedHomology = isConcentratedHomology(self.homologies,self.n) #None if homology is not concentrated, i - int if concentrated in ith homology
        self.injectiveShift = isInjectiveShift(self.homologies,injectivesDimVect,n,self.concentratedHomology) #None, None if not a shift of injective; a,b if equal to I_b[a]
        if self.injectiveShift[0]==m-1:
            self.injective = self.injectiveShift[1]
        else:
            self.injective = None 
        self.projectiveShift = isProjectiveShift(self.homologies,projectivesDimVect,self.n,self.concentratedHomology) #None, None if not a shift of projective; a,b if equal to P_b[a]
        if self.projectiveShift[0]==0:
            self.projective = self.projectiveShift[1]
        else:
            self.projective = None
        self.irrFrom=[] #Is added during the construction loop of the AR-quiver
        self.radicalOfProjective = isRadicalOfProj(self.homologies,projectiveRadicals,self.n,self.concentratedHomology) #None if not the radical of projective, else it is the projective
        self.irrTo=[] #Is added during the construction loop of the AR-quiver
        self.tauInv=None #Is added during the construction loop of the AR-quiver

    def __repr__(self):
        return(str(self.homologies))
    def __str__(self):
        return(str(self.homologies))




    

