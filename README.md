

———————— GENETIC ALGORITHM for BLENDSHAPE REFINEMENT ————————————
								Kyle Reed
							    University of Bath


The purpose of this code is to present a user interface in MAYA acting as a Genetic Algorithm to evolve
the appearance of both identity and smile models. 

To initiate the system please follow the following tasks:

1.) Unzip all folders 
2.) Move ‘Dataset’, ‘Blendshapes’, ‘Resources’ and ‘TedTargets’ folders into the main 
     ‘GeneticAlgorithmBlendshapeCode’ folder.
3.) Run ‘init.m’
	— This is the main MATLAB function to create datasets to be used in the MAYA interface
	— If a file is missing, uncomment sections of this script relating to that file. Some functions
	might take a VERY long time to finish.

—————————————

To run one of the EA scripts

1.) Move ‘GA_scripts’ into your MAYA ‘working root’
2.) Ensure that the Python Path for MAYA includes packages ‘numpy’ and ‘scipy’ 
	— If not, this can be resolved by adding the path to Maya.env file and ‘userSetup.py’
3.) Open ‘GA_scripts/mainMayaInterface.py’ in on of the Maya scripting windows.
4.) To run any 1-7 systems, first highlight the first paragraph which initialises global variables 
      and press the execute/play button
5.) Next highlight the system you are interested in and press the execute/play button


————————————— 

To use the interface:

1.) Check the checkboxes with greatest similarity
2.) Write in the number of THE greatest similarity
3.) Press ‘Next Generation’ button.
4.) Continue until convergence
5.) Save the resulting face by clicking ‘Save Face’ button.


—————————————

To run statistics after a solve:

1.) Data is saved in ‘resultingLogs’ folder in MATLAB workspace
2.) To visualise evolution, change the file name in ‘visualiseChosenEvolution.m’ function and run. 
