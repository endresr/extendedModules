# Extended module categories of linear Nakayama algebras

This project provides python code to generate the postprojective component of the extended module categories of linear Nakayama algebras. The project was developed to assist in the exploration and illustration of the paper *A note on the \(m\)-extended module categories of Nakayama algebras* (https://arxiv.org/abs/2601.14843).

## Features:
- Construct a finite part of the postprojective component of the extended module category
- Print the component in a PDF through Latex
  - Representing the complexes either through their cohomological dimension vectors or as points

## Prerequisites
1. Python 3.8 or higher with Numpy installed
2. A tex installation

## Getting Started
1.	Clone the Repository:

```bash
git clone https://github.com/endresr/extendedModules.git
```
## Usage
Run the main program

```bash
python main.py
```

The program first asks you to input the extended module category you want to calculate.
1. First the amount of simples \(n\) in your Nakayama algebra
2. The extension \(m\) you want
3. The cutoff-value for the calculation. NB! If the resulting category is too large or infinite, a high cutoff-value might lead to the PDF not being generated
  - If the PDF is not generated, you can try the following fixes:
    - Change the cutoff-value
    - Change the scaling-factors
    - Print only nodes instead of the cohomological dimension vectors.
4. The relations of the algebra.
  - If you want homogeneous relations, you give a number less than the amount of simples
  - If you want other relations, input a list of tuples associated to the minimal relations
You are then asked if you want the PDF to be generated with standard configurations. If not, you may choose
5. The y-levels of the \(\tau\)-orbits as a list of numbers
6. The x- and y-scale as a tuple of numbers
7. The scale for the nodes in the tikz-diagram
8. If you want to print the nodes as cohomological dimension vectors or as filled circles.
The PDF is then generated and you get the following choices
1. Redraw using custom settings (You can then change the settings in steps 5.-8.)
2. Start from the beginning
3. Change the cutoff for the initial calculations
4. Quit

## Project Structure

```
extendedModules
├── modules/
│   ├── classes.py      # Extendedmodule class definition and the Tau-function
│   ├── functions.py    # Basic functions for computations
│   ├── drawGraph.py    # Code to construct the Latex file and compiling it
├── Archive/            # Old files 
│   └── 
├── main.py            # Main program
└── README.md
```

## Examples

## License

This project is licensed under the MIT License. See the LICENSE file for details.

