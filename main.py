import numpy as np
from modules.classes import *
from modules.functions import *
from modules.drawGraph import *

class ExtendedModuleCategoryOfNakayama:
    def __init__(self):
        self.n = None #Number of simples in the algebra
        self.m = None #How extended the extended module category is
        self.rel = None #The relations of the algebra
        self.cutoffIterations = None #How many iterations do we do until we either conclude it will continue forever or we are otherwise satisfied
        self.relLength = None #The length of the relations, if they are homogeneous
        self.ylevels = None #What y-level is the tau-orbits printed at
        self.tikzScale = None #Scaling of the resulting tikz-diagram
        self.nodeScale = None #The scale of each node in the tikz diagram
        self.outputName = None #String with your prefered name for Latex-file, if None a predefined name will be given
        self.highlightConcentradedHomology = False #Set to False if you do not want to highlight nodes with concentrated homology
        self.drawOnlyCircles = None #Set to False if you want the homologies printed as node labels

        self.projectiveModules = None
        self.tauOrbits = None
        self.modules = None
        self.customFileName = None
    
    def _get_user_answer_yes_no(self,prompt):
        while True:
            print(prompt)
            answer = input("\nEnter:")
            if answer.lower() =="y":
                return True
                
            elif answer.lower() !="n":
                print(f"\033[1;31m Error: Answer either 'y' or 'n'\033[0m ")
            else:
                return False

    
    def _get_validated_input(self, prompt, validator, error_msg):
        """Generic input validator"""
        while True:
            try:
                print(prompt)
                value = int(input("\nEnter: "))
                if validator(value):
                    return value
                print(f"\033[1;31m Error: {error_msg}\033[0m")    
            except ValueError:
                print(f"\033[1;31m Please enter a valid number\033[0m")
    
    def _get_user_input_start(self):
        """Retrieve the needed info from user"""
        while self.n is None:
            try:
                self.n =self._get_valid_n()
            except ValueError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
            except TypeError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
        while self.m is None:
            try:
                self.m =self._get_valid_m()
            except ValueError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
            except TypeError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
        while self.cutoffIterations is None:
            try:
                self.cutoffIterations =self._get_valid_cutoff()
            except ValueError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
            except TypeError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
        while self.rel is None:
            try:
                self.rel=self._get_valid_relations()
            except ValueError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
            except TypeError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
    
    def _get_valid_n(self):
        return self._get_validated_input(
            "\n\033[1;94mEnter n (Amount of simples of the algebra): \033[0m",
            lambda x: x>0,
            "n must be a positive integer"

        )
    def _get_valid_m(self):
        return self._get_validated_input(
            "\n\033[1;94mEnter m (How long can the complexes be): \033[0m",
            lambda x: x>0,
            "m must be a positive integer"

        )
    def _get_valid_cutoff(self):
        return self._get_validated_input(
            "\n\033[1;94mEnter cutoff for iterations\033[0m \n(If the component is infinite, the iterations will go on forever without cutoff): ",
            lambda x: x>0,
            "The cutoff must be a positive integer"

        )
    def _get_valid_relations(self):
        while True:
            answer = self._get_user_answer_yes_no("\n\033[1;94mDo you want homogeneous relations? (y/n)\033[0m")
            if answer:
                while True:
                    try:
                        print("\n\033[1;94mEnter relation length less than "+str(self.n)+"\033[0m")
                        value =int(input("\nEnter: "))
                        relationsInput = [(i,i+value) for i in range(self.n-value)]
                        self.relLength = value
                        break
                    except:
                        print(f"\033[1;31m Please enter a valid input\033[0m")
            else:
                while True:
                    try:
                        print("\n\033[1;94mEnter relations as a minimal set of generators given by start and end point, as comma separated values e.g. 1,3 2,4\033[0m")
                        value =input("\nEnter: ").split()
                        relationsInput = [tuple(map(int, elem.split(','))) for elem in value]
                        break
                    except:
                        print(f"\033[1;31m Please enter a valid input\033[0m")
            sortedRel = sorted(relationsInput,key=lambda x:x[0])
            tempMax = [b[1] for b in sortedRel]
            if tempMax == sorted(tempMax) and tempMax[-1]<self.n and sortedRel[0][0]>=0:
                return sortedRel
            print(f"\033[1;31m The relations are not in accepted form \033[0m")
            

    def _calculate_PostProjective(self):
        dimProjectives,radicalOfProjectives = findProjectivesDim(self.n,self.rel)
        dimInjectives = findInjectiveDim(self.n,self.rel)[::-1]

        homologiesOfProjectives = findHomologyOfProjectives(dimProjectives,self.n,self.m)
        self.projectiveModules = [extendedModule(x,self.n,self.m,dimProjectives,dimInjectives,radicalOfProjectives) for x in homologiesOfProjectives]
        modules = self.projectiveModules[:]

        #The loop constructing the AR-Quiver
    
        Next = []
        tempVect = np.array([1 for i in range(self.n)])
        for i in range(len(dimProjectives)):
            if np.dot(dimProjectives[i],tempVect)==1:
                Next.append(self.projectiveModules[i]) #Setting our beginning to the simple projectives
                self.projectiveModules[i].xcoord = 0 #Setting the simple projectives x-coordinate to 0 
    
        counter=1
        notFinishedARquiver = False
        continueLoop = True
        while len(Next)>0 and counter<self.cutoffIterations and continueLoop:
            tempNext = []
            for M in Next:
                if (M.homologies<0).any():
                    continueLoop=False
                if M.tauInv == None:
                    tauInv(M,self.projectiveModules,radicalOfProjectives,dimProjectives,dimInjectives,self.n,self.m,modules)
                tempNext+=(M.irrTo)
            Next=list( dict.fromkeys(tempNext) )
            counter+=1
            if counter>=self.cutoffIterations and len(Next)>0:
                print(f"\033[1;31m Reached cutoff-value before finishing.\033[0m")
                notFinishedARquiver = True
            if continueLoop == False and len(Next)>0:
                print("Error: Homologies are non-positive")
        print(f"\033[1;92mFinalized calculations. Found "+str(len(modules))+" modules.\033[0m ")
        self.modules=modules
    
    def _reset_startvalues(self):
        self.n = None 
        self.m = None 
        self.rel = None 
        self.cutoffIterations = None 
        self.relLength = None 
        self.modules = None
        self.projectiveModules = None
        self.tauOrbits = None
        
        
    def _reset_drawingvalues(self):
        self.ylevels = None 
        self.tikzScale = None 
        self.nodeScale = None 
        self.drawOnlyCircles = None 
        self.customFileName = None
        
    def _standard_drawingvalues(self):
        self.ylevels = range(self.n)[::-1]
        self.tikzScale = (2,4)
        self.nodeScale = 1
        self.drawOnlyCircles = False

    def _get_valid_ylevels(self):
        self.ylevels = None
        while self.ylevels is None:
            try:
                print("\n\033[1;94mSet the y-levels for which the tau-orbits are drawn in the AR-quiver.\033[0m \n Enter 'preset' if you want standard levels, or enter a comma-separated list of "+str(self.n)+" values.")
                value = input("\nEnter: ")
                if value=='preset':
                    self.ylevels=range(self.n)[::-1]
                else:
                    temp = [float(x) for x in value.split(",")]
                    if len(temp)!= self.n:
                        raise TypeError(f" \033[1;31m There are "+str(self.n)+"orbits which need a y-level.\033[0m")
                    else: 
                        self.ylevels=temp
            except ValueError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
            except TypeError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
    
    def _get_valid_tikz_scale(self):
        self.tikzScale = None
        while self.tikzScale is None:
            try:
                print("\n\033[1;94mSet the x- and y-scale of the tikz-diagram as a comma-separated list of values.\033[0m\n Recommended values are 2,4 if you print cohomological dimension vectors, .2,.4 if you only print dots.")
                value = input("\nEnter: ")
                temp = tuple(map(float, value.split(',')))
                if len(temp)!=2:
                    raise TypeError("Only two values are needed.")
                self.tikzScale =temp
            except ValueError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
            except TypeError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
    
    def _get_valid_node_scale(self):
        self.nodeScale = None
        while self.nodeScale is None:
            try:
                print("\n\033[1;94mEnter a scaling factor for the nodes in the quiver:\033[0m")
                self.nodeScale = float(input("\nEnter:"))
            except ValueError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
            except TypeError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
    
    def _get_valid_circles(self):
        self.drawOnlyCircles = None
        dimVec = self._get_user_answer_yes_no("\n\033[1;94mDo you want to print the cohomological dimension vectors? (y/n)\033[0m")
        if dimVec:
            self.drawOnlyCircles = False
        else: 
            self.drawOnlyCircles = True
    
    def _get_user_input_draw(self):
        standard = self._get_user_answer_yes_no("\n\033[1;94mDo you want to draw the AR-quiver with standard values? (y/n)\033[0m")
        if standard:
            self._standard_drawingvalues()
        else:
            self._get_valid_ylevels()
            self._get_valid_tikz_scale()
            self._get_valid_node_scale()
            self._get_valid_circles()
        
        
            
    def _draw_quiver(self):
        print("\n\033[1;93mRunning compiler\033[0m\n")
        build_Quiver(self.n,self.m,self.relLength,self.projectiveModules,self.cutoffIterations,self.tikzScale,self.nodeScale,self.highlightConcentradedHomology,self.drawOnlyCircles,self.ylevels)

    def intialize(self):
        print("\n Welcome! \n This script calculates the postprojective component up to a cutoff-point, of linear Nakayama algebras.")
        print("-----------------------")
        self._reset_startvalues()
        self._reset_drawingvalues()
        self._get_user_input_start()
        self._calculate_PostProjective()
        self._get_user_input_draw()
        self._draw_quiver()
    
    def display_menu(self):
        print("\n\033[1;96mMenu options:\033[0m")
        print("1. Redraw using custom settings")
        print("2. Start from beginning")
        print("3. Change cutoff for iterations")
        print("4. Save output files")
        print("5. Quit")
        return input("\n Enter your choice (1-5): ")

    def _redraw_menu(self):
        print("\n\033[1;96mWhat do you want to change?\033[0m")
        print("1. Change y-levels")
        print("2. Set x- and y-scales")
        print("3. Scale nodes")
        print("4. Print cohomological dimension vectors")
        print("5. Go through all settings")
        choice = None
        while choice is None:
            try:
                value = int(input("\n Enter your choice (1-5): "))
                if value>0 and value<6:
                    choice = value
            except:
                print("\n\033[1;91mInvalid choice. Please select 1-4.\033[0m")

        return choice

    def _redraw(self):
        self.customFileName = None
        menu_actions = {
            '1' : self._get_valid_ylevels,
            '2' : self._get_valid_tikz_scale,
            '3' : self._get_valid_node_scale,
            '4' : self._get_valid_circles,
            '5' : self._get_user_input_draw
        }
        choice = self._redraw_menu()
        menu_actions[str(choice)]()
        self._draw_quiver()

    def _new_cutoff(self):
        self.customFileName = None
        self.cutoffIterations = None
        while self.cutoffIterations is None:
            try:
                self.cutoffIterations =self._get_valid_cutoff()
            except ValueError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
            except TypeError as e:
                print(f"\033[1;31m Error: {str(e)}\033[0m")
        self._calculate_PostProjective()
        self._get_user_input_draw()
        self._draw_quiver()
        
    
    def _save_Pdf(self):
        if self.relLength is None:
            stdFileName = str(self.m)+"-mod_Lambda-"+str(self.n)+"-custom_relations"
        else:
            stdFileName = str(self.m)+"-mod_Lambda-"+str(self.n)+"-"+str(self.relLength)
        while self.customFileName is None:
            tempFilename = stdFileName
            standard = self._get_user_answer_yes_no("\n\033[1;94mDo you want to save with the following standard name? (y/n)\033[0m"+"\n  "+stdFileName)
            
            if not standard:
                while True:
                    print("\n\033[1;94mEnter a name for the saved files (without .pdf or .tex)\033[0m")
                    tempInput = input("\nEnter: ")
                    if ".tex" in tempInput or ".pdf" in tempInput:
                        print(f"\033[1;31m Error: Do not enter file extension (.tex or .pdf)\033[0m ")
                    else:
                        tempFilename = tempInput
                        break
            if not check_if_name_taken(tempFilename+".pdf","./AR-quivers/"):
                try:
                    save_graph("preview",tempFilename,"./AR-quivers/")
                    self.customFileName =tempFilename
                    print(f"\033[1;33m The files have been saved as: \033[0m" + "./AR-quivers/"+tempFilename)
                except Exception as e:
                    print(f"\033[1;91mAn unexpected error occurred:\033[0m {e}")
            else:
                overwrite = self._get_user_answer_yes_no(f"\033[1;31m The filename already exists. Overwrite? (y/n)\033[0m")
                if overwrite:
                    try:
                        save_graph("preview",tempFilename,"./AR-quivers/")
                        self.customFileName =tempFilename
                        print(f"\033[1;33m The files have been saved as: \033[0m" + "./AR-quivers/"+tempFilename)
                    except Exception as e:
                        print(f"\033[1;91mAn unexpected error occurred:\033[0m {e}")
    
    def _quit(self):
        if self.customFileName is None:
            quitWithoutSaving = self._get_user_answer_yes_no("\n\033[1;91mFinal version not saved. Continue without saving? (y/n)\033[0m")
            if not quitWithoutSaving:
                self._save_Pdf()
        print("\nThank you! Bye!")
    
    def run(self):
        self.intialize()

        menu_actions = {
            '1' : self._redraw,
            '2' : self.intialize,
            '3' : self._new_cutoff,
            '4' : self._save_Pdf,
            '5': self._quit
        }
        
        while True:
            choice = self.display_menu()
            if choice == '5':
                menu_actions[choice]()
                break
            action = menu_actions.get(choice)
            if action:
                action()
            else:
                print("\n\033[1;91mInvalid choice. Please select 1-5.\033[0m")

if __name__ == "__main__":
    calculator = ExtendedModuleCategoryOfNakayama()
    calculator.run()



            

                



    
    