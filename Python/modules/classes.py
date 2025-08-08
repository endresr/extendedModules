'''
Construct an extended module
'''

from functions import *
import numpy as np

class extendedModule:
    def __init__(self,homologies,n,m,projectivesDimVect,injectivesDimVect,projectiveRadicals):
        self.homologies = homologies
        self.id = None
        self.n = n
        self.xcoord = None
        self.concentratedHomology = isConcentratedHomology(self.homologies,self.n)
        self.injectiveShift = isInjectiveShift(self.homologies,injectivesDimVect,n,self.concentratedHomology)
        if self.injectiveShift[0]==m-1:
            self.injective = self.injectiveShift[1]
        else:
            self.injective = None 
        self.projectiveShift = isProjectiveShift(self.homologies,projectivesDimVect,self.n,self.concentratedHomology)
        if self.projectiveShift[0]==0:
            self.projective = self.projectiveShift[1]
        else:
            self.projective = None
        self.irrFrom=[]
        self.radicalOfProjective = isRadicalOfProj(self.homologies,projectiveRadicals,self.n,self.concentratedHomology)
        self.irrTo=[]
        self.tauInv=None

    def __repr__(self):
        return(str(self.homologies))
    def __str__(self):
        return(str(self.homologies))




    

