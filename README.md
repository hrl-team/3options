# 3 options normalization

This is the repository storing the data and code used to generate behavioural analyses for the paper:   
>__The functional form of value normalization in human reinforcement learning__   
Sophie Bavard, Stefano Palminteri      
[https://elifesciences.org/articles/83891](https://elifesciences.org/articles/83891)

## Data availability
All raw data matrices are provided in .mat and .csv formats. 

>[!NOTE]
>Variable and column names are described below, please contact me if any question arises.        
>sophie[dot]bavard[at]gmail[dot]com

## Behavioral analyses   
Run the *behavioral_analyses.m* script to produce the *data_fig_X.mat* files for each experiment.   
The script loads the raw matrix *data_expeX.mat* and extracts the variables of interest for the figures.

## Model fitting & model simulations
Participants' choices were modeled using likelihood maximization.
### Optimization
Run the *Model_fitting.m* script to produce the *Optimization_X.mat* files. The script loads the raw matrix *data_expeX.mat*, runs model fitting using the functions *function_model_simulations_X.m*, and saves fitted variables in the dedicated files.
### Simulations
Run the *Model_simulations.m* script to produce the *simulation_X.mat* files. The script loads the file *Optimization_X.mat*, runs model simulation using the functions *function_model_simulations_X.m*, and saves fitted variables in the dedicated files.
### Ex-ante simulations
Run the *Model_simulations_exante.m* script to produce the *exante_X.mat* files. The script does not use actual data, it simulates the behavior of agents using the functions *function_model_simulations_X.m*, and saves fitted variables in the dedicated files.

## Generate the figures   
Run the "*Figure...*" scripts to generate the figures.   
The scripts load the files *data_fig_X.mat*, *simulation_X.mat*, *exante_X.mat*.

## Functions   
Files *function_model_simulations_X.m* contain the algorithms used to fit and simulate the data for all models presented in the main text for experiment 1 and 2 or experiment 3 seperately.     
Files *structure_matrix_to_plotmatrix.m* and *vector_to_structure_matrix.m* were created and used to reorganize data structures. Created by Stefano Palminteri, 2018.   
Files *bars_datamodel.m* and *violinplotSB.m*  were created and used for visual purposes. Created by Sophie Bavard, 2021.

## Data   
The columns are ordered as follows:    
* COLUMN 1: participant number
* COLUMN 2: phase number
  * 0: training phase
  * 1: learning phase
  * 2: transfer phase
* COLUMN 3: trial number
* COLUMN 4: context number
  * expe 1 
    * 1 = WT (86/50/14)
    * 2 = WB (86/-/14)
    * 3 = NT (50/32/14)
    * 4 = NB (50/-/14)
  * expe 2
    * 1 = WT (86/50/14)
    * 2 = WB (86/-/14)
    * 3 = NT (86/68/50)
    * 4 = NB (86/-/50)
  * expe 3
    * 1 = WT (86/50/14)
    * 2 = WT (86/50/14)
    * 3 = NT (50/32/14)
    * 4 = NT (50/32/14)
* COLUMN 5: left option
* COLUMN 6: middle option
* COLUMN 7: right option
* COLUMN 8: choice left/middle/right
  * -1: left
  * 0: middle
  * 1: right
* COLUMN 9: choice accuracy
  * 0: incorrect
  * 1: correct
* COLUMN 10: outcome chosen option
* COLUMN 11: outcome unchosen option 1
* COLUMN 12: outcome unchosen option 2
* COLUMN 13: trial reaction time (ms)
* COLUMN 14: participant ID
