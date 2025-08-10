import numpy as np

'''
Helper functions for printing the Latex-code
'''

import subprocess # To compile Latex code from Python


def compile_latex(tex_file_path,outputDirectory):
    """
    Compiles a .tex file to PDF using pdflatex.
        Given by Google's AI
    """
    temp = '-output-directory '+outputDirectory+' '+tex_file_path
    try:
        # Run pdflatex command
        subprocess.run(['pdflatex', '-output-directory', outputDirectory,"-interaction=batchmode",tex_file_path], check=True)
        print(f"Successfully compiled {tex_file_path} to PDF.")
    except subprocess.CalledProcessError as e:
        print(f"Error compiling {tex_file_path}: {e}")
    except FileNotFoundError:
        print("Error: pdflatex command not found. Ensure LaTeX is installed and in your PATH.")

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

def drawNodes(M,x,y,scale=1):
    node =r'\node (' + M.id +') at (' + str(x)+ ','+ str(y)+') [scale='+str(scale)+'] '+r'{$'+bmatrix(M.homologies[::].astype(np.int16))+r'$};'+'\n'
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

def preTikz(scale=(1,1)):
    return r"\begin{tikzpicture}[xscale="+str(scale[0])+",yscale="+str(scale[1])+"]\n"
postTikz = r"\end{tikzpicture}"