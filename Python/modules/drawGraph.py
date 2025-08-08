import numpy as np

def bmatrix(a):
    """
    Returns a LaTeX bmatrix string from a NumPy array.
        Function given by Google AI
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')

    # Convert array elements to strings and join with '&' for columns
    lines = np.array2string(a, max_line_width=np.inf, precision=4, suppress_small=True).replace('[', '').replace(']', '').splitlines()
    
    # Format each line as a LaTeX row
    rv = [r'\begin{bsmallmatrix}']
    rv += [' ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv += [r'\end{bsmallmatrix}']
    
    return '\n'.join(rv)

def drawNodes(M,x,y):
    node =r'\node (' + M.id +') at (' + str(x)+ ','+ str(y)+') [] '+r'{$'+bmatrix(M.homologies[::-1].T.astype(np.int16))+r'$};'+'\n'
    return node

def drawIrrArrows(M):
    arrows =''
    for N in M.irrTo:
        arrows += r'\draw[-latex] ('+ M.id+') -- ('+ N.id+');\n'
    return arrows

def drawTauArrow(M):
    return r'\draw[dashed] ('+ M.id+')--('+ M.tauInv.id+');\n'

preLatex = r"\documentclass[margin=2mm]{standalone}"+"\n"+r"\usepackage{tikz}"+"\n"+r"\usepackage{amsmath,amssymb,mathtools}"+"\n"+r"\begin{document}"+"\n"
postLatex = r"\end{document}"

preTikz = r"\begin{tikzpicture}"+"\n"
postTikz = r"\end{tikzpicture}"