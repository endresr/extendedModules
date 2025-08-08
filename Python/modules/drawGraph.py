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
    rv = [r'\begin{bmatrix}']
    rv += [' ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv += [r'\end{bmatrix}']
    
    return '\n'.join(rv)

print(bmatrix(np.array([[[1],[2],[3]],[[4],[5],[6]],[[7],[8],[9]]])))