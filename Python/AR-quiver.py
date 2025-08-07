import numpy as np

n=5
m=3

relations = [(1,3),(3,5)]#list of minimal relations
#print(relations[::-1])

def findProjectives(n,relations):
    projectives = []
    for a in range(1,n+1):
        upper = next((r[1] for r in relations if a<=r[0]),n+1)-1
        
        projectives.append(np.array([1 if a<=i and i<=upper else 0 for i in range(1,n+1)]))
    return projectives
#print(findProjectives(n,relations) )

def findSimples(n):
    simples = []
    for i in range(1,n+1):
        simples.append(np.array([1 if a==i else 0 for a in range(1,n+1)]))
    return simples
#print(findSimples(n))

print([i for i in range(n,0,-1)])

def findInjectives(n):
    injectives = []
    temp_rel = relations[::-1]
    for i in range(n,0,-1):
        lower = next((r[0] for r in temp_rel if i>=r[1]),0)+1
        injectives.append(np.array([1 if lower<=j and j<=i else 0 for j in range(1,n+1)]))
    return(injectives)
#print(findInjectives(n))

#if (np.zeros(2)==np.array([0,0])).all():
#    print("hei")

def isConcentratedHomology(M):
    '''
    This function checks if the extended module is concentrated in one homology, i.

    Input: M - Extended module

    Output: i, if true, otherwise None
    '''

def isInjectiveShift(M):
    '''
    This functions first check if the extended module is concentrated in a single homology i, 
    if so it checks if this homology is that of a shifted injective I_j. 

    Input: M - Extended module

    Output: (a,b) tuple
        a = i and b = j if homology is of shifted injective, otherwise a= None and b=None
    '''
    
def isRadicalOfProj(M):
    '''
    If module is concentrated in homology 0, then it checks if it is the radical of any projectives.

    Input: M - Extended Module
    Output: L
        L is list of projectives having M as radical or None
    '''

def almostSplitFrom(M):
    '''
    Finds the almost split sequence starting in M
        M -(f)-> E_1 + E_2 + ... + E_n -(g)-> tau- M=L

    - If M is shifted Injective I_i[j], then L=P_i[j+1] and E_k is found through the almost split sequences starting in tau E_k.
    Else
    - L is found through H^i(E_1)+H^i(E_2)+...+H^i(E_n)-H^i(M)
    where E_k is given by
        * If H^i(M)=0 for i!=0 and it is a radical of projectives P_l_1, P_l_2, ... , P_l_t, then E_k=P_l_k for 1<=k<=t, 
    and E_k for t+1<=k<=n is found through the almost split sequences starting in tau E_k
        * Otherwise E_k is given as above for 1<=k<=n
    '''
    irredMorIn = M.irrTo()
    

class extendedModule:
    def __init__(self,homologies,n,proj,inj,irrMorIn,irrMorOut):
        self.homologies = homologies
        self.n = n
        self.projective = proj
        self.injective = inj
        self.irrTo = irrMorIn #irreducible morphisms ending in module
        self.irrFrom = irrMorOut #irreducible morphisms starting in module
        if len(filter(lambda x:(x!=np.zeros(n)).all(),self.homologies)) !=1: