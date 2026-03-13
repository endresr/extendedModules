'''
Helper functions for printing the Latex-code
'''
from datetime import datetime
import subprocess
import shutil
import os
import numpy as np
from pathlib import Path

def compile_latex_with_latexmk(filename_without_extension, output_dir=None,input_dir=None):
    """
    Compiles a .tex file using latexmk and optionally moves the output to a specified directory.
    
    Args:
        filename_without_extension (str): The name of the LaTeX file without the .tex extension.
        output_dir (str, optional): The destination directory for the final PDF.
    """
    if input_dir is None:
        source_file = f"{filename_without_extension}.tex"
    else:
        source_file = f"{input_dir}{filename_without_extension}.tex"
    
    
    # Check if the source file exists
    if not os.path.exists(source_file):
        print(f"Error: {source_file} not found.")
        return

    # Use a temporary directory for compilation to keep the source directory clean
    # The 'latexmk' command uses '-output-directory' to place all auxiliary files there.
    if output_dir is None:
        temp_dir = 'temp_latex_build'
    else:
        temp_dir = output_dir+'temp_latex_build'    
    os.makedirs(temp_dir, exist_ok=True)
    
    command = ["latexmk", "-pdf", "-interaction=batchmode","-silent", f"-output-directory={temp_dir}", source_file]
    
    try:
        # Run latexmk
        subprocess.check_call(command)
        print(f"\033[1:32mSuccessfully compiled {source_file}.\033[0m")

        # If an output directory is specified, move the final PDF
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            pdf_filename = f"{filename_without_extension}.pdf"
            source_path = os.path.join(temp_dir, pdf_filename)
            destination_path = os.path.join(output_dir, pdf_filename)
            
            shutil.move(source_path, destination_path)
            print(f"Moved PDF to {destination_path}.")
        # Optional: Clean up the temporary directory after the process
        # Use latexmk's cleanup functionality for a complete job
        cleanup_command = ["latexmk", "-c", f"-output-directory={temp_dir}", source_file]
        subprocess.call(cleanup_command)
        if os.path.exists(temp_dir) and not os.listdir(temp_dir): # Only remove if empty
             shutil.rmtree(temp_dir)
    except subprocess.CalledProcessError as e:
        print(f"\033[1;91mLaTeX compilation failed for {source_file}.\033[0m")
        print(f"\033[1;91mError details:\033[0m {e}")
        print("\n The file might be too large. Try changing from printing cohomological dimension vectors, change to smaller scale factors or a smaller cutoff.")
    except Exception as e:
        print(f"\033[1;91mAn unexpected error occurred:\033[0m {e}")
        print("\n The file might be too large. Try changing from printing cohomological dimension vectors, change to smaller scale factors or a smaller cutoff.")

def bmatrix(a):
    """
    Returns a LaTeX bmatrix string from a NumPy array.
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

def moduloMatrix(M,k):
    m=M.shape
    tempMatrix = np.zeros((m[0],k))
    for i in range(m[0]):
        for j in range(m[1]):
            tempMatrix[i][j%k]+=M[i][j]
    return tempMatrix

def drawNodes(M,x,y,scale=1,highlightConcHom=False,drawOnlyCircles=False):
    node = r'\node (' + M.id +') at (' + str(x)+ ','+ str(y)+') [scale='+str(scale)
    circleNode = r',draw,fill,circle,inner sep=0pt,minimum size=4pt] {};' +'\n'
    matrixNode =  r'] {$'+bmatrix(M.homologies[::].astype(np.int64))+r'$};'+'\n'
    
    if M.concentratedHomology != None and highlightConcHom:
        if drawOnlyCircles:
            node = node + ",red" + circleNode
        else:
            node = node + ',rectangle,draw=red,line width=2pt'+matrixNode
    else:
        if drawOnlyCircles:
            node = node + circleNode
        else:
            node = node + matrixNode
    return node

def drawIrrArrows(M):
    arrows =''
    for N in M.irrTo:
        arrows += r'\draw[-latex] ('+ M.id+') -- ('+ N.id+');\n'
    return arrows

def drawTauArrow(M):
    return r'\draw[dashed] ('+ M.id+')--('+ M.tauInv.id+');\n'

preLatex = r"%&pdflatex"+"\n"+r"\documentclass[margin=2mm]{standalone}"+"\n"+r"\usepackage{tikz}"+"\n"+r"\usepackage{amsmath,amssymb,mathtools}"+"\n"+r"\begin{document}"+"\n"
postLatex = r"\end{document}"

def preTikz(scale):
    return r"\begin{tikzpicture}[xscale="+str(scale[0])+",yscale="+str(scale[1])+"]\n"
postTikz = r"\end{tikzpicture}"

def build_Quiver(n,m,rel,projectiveModules,cutoff,tikzScale,nodeScale,hightlightConcHom,drawCircles,yLevels):
        tauOrbits = []
        for i in range(len(projectiveModules)):
            temp=projectiveModules[i]
            orbit = [temp]
            counter=0
            temp.id=f"t-{counter}P{i}"
            while temp.injective == None and counter < cutoff and temp.tauInv !=None:
                counter +=1
                temp=temp.tauInv
                temp.id=f"t-{counter}P{i}"
                orbit.append(temp)
            tauOrbits.append(orbit)
        TikzNodes = ""
        TikzIrrArrows = ""
        TikzTauArrows = ""

        xMax=0
        for i in range(n):
            y=yLevels[i]
            for j in range(len(tauOrbits[i])):
                M=tauOrbits[i][j]
                TikzNodes += drawNodes(M,M.xcoord,y,nodeScale,hightlightConcHom,drawCircles)
                TikzIrrArrows += drawIrrArrows(M)
                if M.xcoord!=None and M.xcoord>xMax:
                    xMax=M.xcoord
                if M.tauInv != None:
                    TikzTauArrows += drawTauArrow(M)
        
        stringToSave = preLatex + preTikz(tikzScale)+TikzNodes+TikzIrrArrows+TikzTauArrows+postTikz+postLatex
        outputDirectory = './'#"./AR-quivers/"
        currentTime = datetime.today().strftime('%Y-%m-%d')
        #if rel is None:
        #    fileName = currentTime+"_"+str(m)+"-mod_Lambda-"+str(n)+"-custom_relations"
        #else:
        #    fileName = currentTime+"_"+str(m)+"-mod_Lambda-"+str(n)+"-"+str(rel)
        fileName ='preview'
        with open(outputDirectory+fileName+".tex","w") as file:
            file.write(stringToSave)
        compile_latex_with_latexmk(fileName,outputDirectory,outputDirectory)

def save_graph(oldName,fileName,outputDirectory):
    shutil.move('./'+oldName+'.tex',outputDirectory+fileName+".tex")
    shutil.move('./'+oldName+'.pdf',outputDirectory+fileName+".pdf")

def check_if_name_taken(filename,dir):
    filepath = Path(dir + filename)
    return filepath.exists()