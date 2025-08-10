import numpy as np

def findProjectivesDim(n,relations):
    projectives = []
    projRadical = []
    for a in range(0,n):
        upper = next((r[1] for r in relations if a<=r[0]),n)-1
        proj=[1 if a<=i and i<=upper else 0 for i in range(0,n)]
        projRad=[1 if a<i and i<=upper else 0 for i in range(0,n)]
        projectives.append(np.array(proj))
        projRadical.append(np.array(projRad))
    return projectives,projRadical

def findHomologyOfProjectives(dimVecs,n,m):
    proj = []
    for v in dimVecs:
        tempProj = np.append([v],np.zeros((m-1,n)),axis=0)
        proj.append(tempProj)
    return proj


def findInjectiveDim(n,relations):
    injectives = []
    temp_rel = relations[::-1]
    for i in range(n)[::-1]:
        lower = next((r[0] for r in temp_rel if i>=r[1]),-1)+1
        tempInj = [1 if lower<=j and j<=i else 0 for j in range(0,n)]
        injectives.append(np.array(tempInj))
    return(injectives)

def findHomologyOfInjectives(dimVecs,n,m):
    inj = []
    for v in dimVecs:
        tempInj = np.append(np.zeros((m-1,n)),[v],axis=0)
        inj.append(tempInj)
    return inj

def isConcentratedHomology(homologies,n):
    '''
    This function checks if the extended module is concentrated in one homology, i.

    Input: M - Extended module

    Output: i, if true, otherwise None
    '''
    v=np.array([1 for i in range(n)])
    test =  np.flatnonzero(np.dot(homologies,v)>0)
    if len(test)==1:
        return test[0]
    else:
        return None



def isInjectiveShift(homologies,I,n,homDegree):
    '''
    This functions first check if the extended module is concentrated in a single homology i, 
    if so it checks if this homology is that of a shifted indecomposable injective I_j. 

    Input: M - Extended module, I - array of injectives in mod(Lambda) with rows equal to the dimensionvector of each injective
        underlying assumption is that the index of the rows correspond to the index of the injective

    Output: (a,b) tuple
        a = i and b = j if homology is of shifted injective, otherwise a= None and b=None
    '''
    if homDegree==None:
        return None,None
    else:
        tempOut = np.flatnonzero((homologies[homDegree]==I).all(1))
        if len(tempOut)>0:
            return homDegree,tempOut[0]
        else:
            return None,None
        
def isProjectiveShift(homologies,P,n,homDegree):
    '''
    This functions checks if the extende module is a shifted indecomposable projective 
    and in that case returns the homology, i, and the index of the projective, j.

    Input: M - Extended module, I - array of indec projectives in mod(Lambda) with rows equal to the dimensionvector of each projective
        underlying assumption is that the index of the rows correspond to the index of the projective

    Output: (a,b) tuple
        a = i and b = j if it is shifted projective, else a = None and b = None
    '''
    if homDegree == None:
        return None, None
    else: 
        tempOut = np.flatnonzero((homologies[homDegree]==P).all(1))
        if len(tempOut)>0:
            return homDegree,tempOut[0]
        else:
            return None,None
    
def isRadicalOfProj(homologies,projRadical,n,homDegree):
    '''
    If module is concentrated in homology 0, then it checks if it is the radical of any projectives.

    Input: M - Extended Module, projRadical - array of the dim vector of radicals
    Output: L
        L is list of projectives having M as radical or None
    '''
    if homDegree == 0:
        tempOut = np.flatnonzero((homologies[0]==projRadical).all(1))
        if len(tempOut)>0:
            return tempOut[0]
        else:
            return None
    else:
        return None
    
def checkRelations(n,rel):
    '''
    Sort the relations and check if they are minimal
    '''
    sortedRel= sorted(rel,key=lambda x:x[0])
    print(sortedRel)
    tempMax = [b[1] for b in sortedRel]
    if tempMax != sorted(tempMax) or tempMax[-1]>=n or sortedRel[0][0]<0:
        return None
    else:
        return sortedRel
