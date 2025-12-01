"""
Main program for running the POTPRIM model in different modes.

This script allows running the POTPRIM (Potential Priming) model in various modes:
- Normal: Standard model run with default parameters
- Sensitivity: Sensitivity analysis of model parameters 
- Bayesian: Bayesian calibration of model parameters
- Validation: Model validation against experimental data

The program can work with different datasets:
- Jilkova2022: Dataset from Jilkova et al. 2022 publication
- Jilkova2024: Dataset from Jilkova et al. 2024 publication

The script handles parameter initialization, data loading, model execution,
and result processing depending on the selected mode. For Bayesian mode,
it implements parameter optimization using Bayesian methods.

Usage:
    Run the script and specify mode and dataset variables at the top.
    Additional command line arguments are available for Bayesian mode:
    -p/--parallel: Run in parallel mode
    -f/--fields: Maximum number of fields to run
    -t/--tries: Number of optimization tries
    -d/--debug: Print debug information
"""

# stuff needed for Bayesian mode
import argparse
# import concurrent.futures
import json

import warnings

# import math
import os
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import random as ra
from scipy import stats
from scipy.stats import qmc
import BayesianFunctionsPotprim #needed in Bayesian mode
import MainFunctionsPotprim #needed in other modes, plotting, validation etc.
import sys
import csv
from datetime import datetime
import glob
import shutil
from matplotlib.lines import Line2D
from collections import defaultdict

# needed by all modes / normal, sensitivity and bayesian mode
from potPrimingMAOMfunction import *
from typing import Literal, get_args

# needed for sensitivity and bayesian
import copy

#this will stop execution right away when this error arises
warnings.simplefilter("error", RuntimeWarning)
warnings.simplefilter("error", DeprecationWarning)
warnings.simplefilter("error", FutureWarning)
warnings.simplefilter("error", ResourceWarning)
warnings.simplefilter("error", SyntaxWarning)

############## Modes #################################
# set allowed values for mode
modes = Literal["Normal", "Sensitivity", "Bayesian", "Hypercube"]
options = get_args(modes)

#%% set the mode to Normal, Sensitivity, Bayesian or Hypercube (Hypercube can be used for crossvalidation or scenarios simulation)
mode_ = "Bayesian"

# check if mode was set correctly, if not stop the run
assert mode_ in options, f'"{mode_}" is not in "{options}"'

############## Datasets #################################
# set allowed values for dataset
datasets = Literal["Jilkova2022", "Jilkova2024"]
options = get_args(datasets)

# set the dataset to Jilkova 2022 or Jilkova 2024
dataset_ = "Jilkova2022"

# check if dataset was set correctly, if not stop the run
assert dataset_ in options, f'"{mode_}" is not in "{options}"'

# if the Normal or Jilkova2022 mode was chosen, you can decide to turn on the Plotting
Plotting = True

# safety condition / do not allow plotting with Sensitivity or Bayesian mode
if mode_ == "Sensitivity" or mode_ == "Bayesian":
    Plotting = False

########### Folders setup ############################
# the_current_path = os.path.abspath(os.getcwd())
results_path = "./output"
sharable_path = "./output_Bayesian"

try:
    os.makedirs("./output_Bayesian")
except FileExistsError:
    # directory already exists
    pass

############## Read data #############################
# read the fixed parameter list

with open("datalistInput.json", "r", encoding="utf-8") as inputfileParam:
    AllParam = json.load(inputfileParam)
    
# inputfileParam = open(
#     "datalistInput.json"
# )  # input parameters (all) is always in same filenam
# AllParam = json.load(inputfileParam)


# Create variables
df_list = []
results_df = []  # temporary df to store returned dataframe

#%% define normal run which includes simple validation
#first define normal run as function so that it can be also used in Normal mode as well as at the end of Bayesian for immediate validation of the new calibration    
def normal_run(path_normal, 
               duration, 
               cols_measured_respSoil, 
               cols_measured_respSubstrate,
               cols_data_measured,
               cols_measured_PE,
               timepoints,
               AllParam, 
               Plotting,
               results_path
               ):
    
        # load input data - Treatments
        inputRun = pd.read_csv(path_normal, header=0, skiprows=0)    
        numTreatments = len(inputRun)
        #how many rows contain control treatments (these should come first in the input file)
        numControls = inputRun['treatment'].str.contains('control', case=False, na=False).sum()
        
        # create lists for respiration plot
        respSoil_mean_measure = []
        respSubstrate_mean_measure = []
        labels = []
        
        #%%--- run the model
        #run the model for all treatments/rows in treatment input file

        for treatment in range(numTreatments):
            treatmentVar = inputRun.iloc[treatment, 0:21] # select first 21 columns from the input file
            # print(treatment)
            results_df = run_model(
                AllParam,
                treatmentVar,
                mode_="Normal",
                Plotting=Plotting,
                numDays=duration,
                path=results_path
            )
            results_df['treatmentID'] = treatment + 1
            df_list.append(results_df)         
                            
            if Plotting:
                # storing values for respiration plot
                labels.append(results_df["treatment"][1]) # treatment label
                #eventually this could be also obtained at the end of the run somehow, not throughout:
                respSoil_mean_measure.append((inputRun.iloc[treatment, cols_measured_respSoil]).mean())
                respSubstrate_mean_measure.append((inputRun.iloc[treatment, cols_measured_respSubstrate]).mean())

        final_results_df = pd.concat(
            df_list, ignore_index=True
        )  # add all the rows to the results_df                        
        
        
        #%%--- validation plots
        #if plotting wanted, calculate also performance metrics
        if Plotting:
            #filter out modelled values for all those variables and days for which we have measured values
            #first automatically extract for which data we have measured data
            data_measured = inputRun.iloc[:, cols_data_measured]  # get measured data
            #separately get measured PE data, because these are not used to extract simulated values
            measured_PE = inputRun.iloc[:, cols_measured_PE]
            
            #now proceed with working with the measured data used to get simulated data
            data_measured_names = []
            data_measured_days = []
            data_measured_colnames = data_measured.columns.tolist() #extract column names of measured data
            
            data_measured_names, data_measured_days = (
                BayesianFunctionsPotprim.split_alphanumeric_list(data_measured_colnames)
            )           # from the column names extract variable name and day of measurement

         
            unique_variables = list(set(data_measured_names)) 

            

            #then transform the modelled data from wide format into long, to be able to filter by combinations of variable and day
            df_long = pd.melt(
                final_results_df,
                id_vars=['treatmentID','day','treatment'],  # Columns to keep as is
                value_vars=unique_variables,  # Columns to unpivot
                var_name='variable',  # Name for the new column containing former column names
                value_name='value'  # Name for the new column containing values
                )
           
            python_data_measured_days = [x - 1 for x in data_measured_days] 
            
            #pairs of variables and days for which we have measurements   
            key = pd.DataFrame({'data_measured_names': data_measured_names, 
                          'data_measured_days': python_data_measured_days
                          })
           
            #filter the modelled data by the measured variables and days by inner join             
            merged_df = pd.merge(
            df_long, 
            key,
            left_on=['variable', 'day'],
            right_on=['data_measured_names','data_measured_days'],
            how='inner'
            )

            #average by day and treatment
            averages = merged_df.groupby(['day','treatment',"variable"])['value'].mean().reset_index() #averages for each variable across all replicates for each treatment
            
            #%%------Respiration plot###########################################################################
            #filter out respiration data
            respiration = averages[averages["variable"].isin(["respSoil", "respSubstrate"])]
            #only for respiration data, calculate the average across days (for the respiration plot)
            resp_avg = respiration.groupby(['treatment',"variable"])['value'].mean().reset_index()
            #convert KEYLINK units to normal (data) units:
            # resp_avg["value"] = resp_avg["value"]  / 0.8 * 24
            
            #convert the long format back to width to have separate columns for respSoil and respSubstrate
            resp_avg = resp_avg.pivot(
                index='treatment',        # Column(s) to use as the index
                columns='variable', # Column whose unique values will become column names
                values='value'     # Column whose values will fill the new DataFrame
            ) 
            #store the measured respiration data into a dataframe
            resp_measured = pd.DataFrame({'treatment': labels, 
                                           'respSoil_mean_measure': respSoil_mean_measure,
                                           'respSubstrate_mean_measure': respSubstrate_mean_measure,                                         
                                            })
            #calculate averages for each treatment
            resp_measured_avg = resp_measured.groupby('treatment').mean()
            
            #merge modelled and measured respiration data
            resp_df = pd.merge(
            resp_avg, 
            resp_measured_avg,
            on=['treatment'],
            how='left'
            )
            pd.set_option('display.max_columns', None)
            figures_path = os.path.join(results_path, "figures")
            
            #make a directory for the figures
            try:
                os.makedirs(figures_path)
            except FileExistsError:
                # directory already exists
                pass
           
            #draw the respiration plot
            MainFunctionsPotprim.drawRespPlot(
                resp_df.index.tolist(), #extract the treatment labels from the index of the dataframe
                resp_df['respSoil'], #mean of modelled
                resp_df['respSoil_mean_measure'],
                resp_df["respSubstrate"],#mean of modelled
                resp_df["respSubstrate_mean_measure"],
                figures_path,
            )
            #%%------1:1 plots and model performance metrics###################
            #first prepare the modelled data in the same format as the input data (with same columns names)
            # add 1 to all "day" values, so it matches the naming in input day
            merged_df['day'] = merged_df['day'] + 1
            #pivot all modelled values wider 
            # Step 1: Pivot
            
            wide_df = merged_df.pivot(index=['treatmentID', "treatment"], columns=['variable', 'day'], values='value')
            
            # Step 2: Flatten MultiIndex columns
            wide_df.columns = [f'{var}{day}' for var, day in wide_df.columns]
            
            # Step 3: Reset index and rename dataframe
            data_modelled = wide_df.reset_index()               
                       
            #calculate priming effects from modelled values 
            # Loop through each timepoint column
            for tp in timepoints:
                col_name = f'respSoil{tp}'
                pe_col_name = f'PE{tp}'
            
                # Calculate baseline: average of control rows
                control_avg = data_modelled[col_name].iloc[:numControls].mean()
            
                # Subtract baseline from all rows to get priming effect
                data_modelled[pe_col_name] = data_modelled[col_name] - control_avg
    
            
            
            # Add this dataframe to the original DataFrame of measured values
            data_measured = pd.concat([data_measured, measured_PE], axis=1)
    
            
            #plot the 1:1 plots
            columns_to_plot = data_measured.columns
            #prepare a colour paletter (this could be customized for Jilkova2024 when needed)
            custom_palette = {
                'control': 'grey',
                'leachates': 'green',
                'exudates': 'orange',                
                'exudates+leachates': 'brown'
            }
            metrics_list = []
            # --- 3. Plot in 2x2 grids ---
            for i in range(0, len(columns_to_plot), 4):

                subset = columns_to_plot[i:i+4]
                
                fig, axes = plt.subplots(2, 2, figsize=(10, 10))
                axes = axes.flatten()
                

                
                # List to collect all unique treatments for the legend
                all_treatments = list(custom_palette.keys())
                
                # Create custom legend handles
               
                legend_handles = [Line2D([0], [0], marker='o', color='w', 
                                      markerfacecolor=custom_palette[treatment], 
                                      markersize=8) for treatment in all_treatments]
                
                for j, col in enumerate(subset):#loop over all four variables for this 2x2 graph
                    MainFunctionsPotprim.plot_1to1(axes[j], data_measured[col], data_modelled[col], data_modelled['treatment'], col, custom_palette)
                    # Calculate summary/diagnostic metrics for this variable
                    rmse, bias, ef = MainFunctionsPotprim.calculate_metrics(data_measured[col], data_modelled[col])    

                    # Add metrics text to each subplot
                    metrics_text = f'RMSE = {rmse:.3f}\nBias = {bias:.3f}\nEF = {ef:.3f}'
                    axes[j].text(0.05, 0.95, metrics_text, transform=axes[j].transAxes, 
                             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
                    metrics_list.append({'Variable': col, 'RMSE': rmse, 'Bias': bias, 'EF': ef})
                # Hide unused axes if fewer than 4
                for k in range(len(subset), 4):
                    fig.delaxes(axes[k])
                
                fig.tight_layout(rect=[0, 0.05, 1, 1])  # Leave space at the bottom for the legend
                # Add a single legend for the entire figure
                # Add a single legend for the entire figure
                fig.legend(handles=legend_handles, labels=all_treatments, 
                           loc='lower center', bbox_to_anchor=(0.5, 0), 
                           ncol=len(custom_palette), frameon=True, title='Treatment')
                
                filename = f'1to1_plots_{i//4 + 1}.png'
                filename2 = os.path.join(figures_path, filename)
                plt.savefig(filename2, dpi=300)
                plt.show()
                plt.close(fig)  # Close the figure to free memory
            # Save metrics
            metrics_df = pd.DataFrame(metrics_list)
            mean_ef = metrics_df['EF'].mean()
            metrics_df = pd.concat([metrics_df, pd.DataFrame([{'Variable': 'Mean EF', 'RMSE': np.nan, 'Bias': np.nan, 'EF': mean_ef}])], ignore_index=True)
             

           #save metrics from this run in the output folder
            csv_file = os.path.join(results_path, 'model_performance_metrics.csv')
            if os.path.exists(csv_file):
                metrics_df.to_csv(csv_file, mode='a', header=False, index=False)
            else:
                metrics_df.to_csv(csv_file, index=False) 
        else:
            mean_ef = []
            metrics_df = []
        return mean_ef, metrics_df, final_results_df #return mean EF as overall performance metric
        # end of normal run function
        
#%% Normal run with validation #####################################################################
if mode_ == "Normal":
        
    if dataset_ == "Jilkova2022":
        path_normal = "Normal_run_input_2022.csv"
        duration = 155 #number of days of incubation
        cols_measured_respSoil = slice(29, 37) #which columns contain measured soil derived respiration
        cols_measured_respSubstrate = slice(37, 45)  #which columns contain measured substrate derived respiration
        cols_data_measured = slice(21, 45)
        cols_measured_PE = slice(69, 77)
        timepoints = [1, 15, 29, 43, 71, 99, 127, 155] #timepoints for which to calculate the PE effect
        
    if dataset_ == "Jilkova2024":
        path_normal = "Normal_run_input_2024.csv"
        duration = 161 #number of days of incubation
        cols_measured_respSoil = slice(33, 53) #which columns contain measured soil derived respiration
        cols_measured_respSubstrate = slice(66, 77)  #which columns contain measured substrate derived respiration
        cols_data_measured = slice(21, 65) #columns with all measured data
   
                             


#run it now
    mean_ef, metrics_df, final_results_df = normal_run(path_normal, 
                   duration, 
                   cols_measured_respSoil, 
                   cols_measured_respSubstrate,
                   cols_data_measured,
                   cols_measured_PE,
                   timepoints,
                   AllParam, 
                   Plotting,
                   results_path
                   )
    
    print(metrics_df)     
    print('''
          0.75 < EF < 1  very good model performance
          0.5 < EF < 0.65 satisfactory to good
          EF < 0 poor performance (model predictions worse than the mean of observations)
          ''')
    #make output directory    
                      
    try:
            os.makedirs("./output/data")
    except FileExistsError:
        # if directory already exists
        pass
    
    #save all modelled data
    final_results_df.to_csv(
        "./output/data/Normal.csv",
        index=False,
        float_format="%.5f",
    )
    
#%% Sensitivity ###########################
if mode_ == "Sensitivity":
    t1 = time.perf_counter()
    if dataset_ == "Jilkova2022":
        path_sensitivity = "Sensitivity_run_input_2022.csv" #use simplified input with just 4 treatments (no replicates)
        path_initial_values = "datalistCalibrationParam_start_overallJuly28.json"
        # cols_data_measured = slice(21, 45) #range of the columns to be considered
        # cols_data_measured_errors = slice(45, 69) #range of the columns to be considered
        # #these extra ones are needed for validation step
        duration = 155 #number of days of incubation
        # cols_measured_respSoil = slice(29, 37) #which columns contain measured soil derived respiration
        # cols_measured_respSubstrate = slice(37, 45)  #which columns contain measured substrate derived respiration
        # cols_measured_PE = slice(69, 77)
        # cols_measured_PE_errors = slice(77, 85)
        # timepoints = [1, 15, 29, 43, 71, 99, 127, 155] #timepoints for which to calculate the PE effect
    
    if dataset_ == "Jilkova2024":
        path_sensitivity = "Bayesian_run_input_2024.csv"
        # cols_data_measured = slice(21, 65) #range of the columns to be considered
        # cols_data_measured_errors = slice(65, 109) #range of the columns to be considered
        # #these extra ones are needed for validation step
        duration = 161 #number of days of incubation
        # cols_measured_respSoil = slice(33, 53) #which columns contain measured soil derived respiration
        # cols_measured_respSubstrate = slice(66, 77)  #which columns contain measured substrate derived respiration

    # load parameter values fromthe calibration input file (from which the overall calibration started)    
    #load 
    with open(path_initial_values) as inputCalibrationParamfile:
        (numParams,
         CalibParamInit, #this will be the values at which other parameters are held constant
         CalParameters,
         CalParameterValues,
         MaximumOption,
         MinimalOption,
         keys,
         ) = BayesianFunctionsPotprim.read_parameter_data(inputCalibrationParamfile)
        
    # paramsToTestDict = CalibParamInit
    paramsToTestNames = keys
    
    # AllParam.update(paramsToTestDict) #update the calibrated parameters with the initial values from calibration
    Hypercube = pd.read_csv("Hypercube_sampled.csv", header=0, skiprows=0)
    Hypercube.columns = keys
    
    #save the hypercube sample
    Hypercube.to_csv(
        os.path.join("Hypercube_sampled_head.csv"),
        index=False,
        float_format="%.5f",
    )
    # calculate medians of hypercube for each parameter
    Median_dict = Hypercube.median().to_dict()  
    AllParam.update(Median_dict)
    
    # Treatments
    inputRun = pd.read_csv(path_sensitivity, header=0, skiprows=0)
    numTreatments = len(inputRun)
    df_list = []
    run_info_list = []
    # paramsToTestNames = ['bact_rhiz_rel', 'fungi_rhiz_rel']
    # paramsToTestNames = ['bact_rhiz_rel']
    for param in paramsToTestNames: #for each parameter
        # if param == 'DOM_EC':
        #     break
        i=0 #count which of the 15 variants is used
        for paramValue in Hypercube.loc[:,param]: #loop over all 15 parameter values from the hypercube 
            
            # # calculate by how much to change the parameter value, using a relative parameter change
            # delta = (
            #     paramsToTestDict[param] * paramChange / 100
            # )  # I change 1 parameter value
            # # caculate new value of parameter
            # value = paramsToTestDict[param] + delta
            # change the value directly in the parameter set then used by run_model

            AllParam[param] = paramValue

            for treatment in range(numTreatments):
                treatmentVar = inputRun.iloc[treatment, 0:21]

                results_df = run_model(
                    AllParam,
                    treatmentVar,
                    mode_,
                    Plotting,
                    numDays=duration,
                    path=None
                )

                
                # add to the simulated values information about the parameter, its change and value
                # temp_df_list = [param, paramChange, value] + temp_df_list

                length = len(results_df)

                info_df = pd.DataFrame(
                    {
                        "param": np.full(length, param),
                        "param_variant": np.full(length, i),
                        "value": np.full(length, paramValue),
                    }
                )

                results_df = info_df.join(results_df)
                df_list.append(
                    results_df
                )  # append doesn't work for dataframes, so the dataframes have to be appended to a list to later use concat
                # save parameter set for each run
                All_param_series = pd.Series(AllParam)
                info_series = pd.Series(
                    [param, i, paramValue], index=["param", "param_variant", "value"]
                )
                # Concatenate the three Series
                run_info = pd.concat([info_series, treatmentVar, All_param_series])
                run_info_list.append(run_info)
            #after each variant tried for certain parameter, keep track of the no of variant    
            i = i+1 

        # after all changes tried for certain parameter, reset its value to original median value
        AllParam[param] = Median_dict[param]

    final_results_df = pd.concat(
        df_list, ignore_index=True
    )  # add all the rows to the results_df
    run_info_df = pd.DataFrame(run_info_list)

    try:
        os.makedirs("./output/data")
    except FileExistsError:
        # directory already exists
        pass

    final_results_df.to_csv(
        "./output/data/Sensitivity.csv",
        index=False,
        float_format="%.5f",
    )

    run_info_df.to_csv(
        "./output/data/Sensitivity_runs.csv",
        index=False,
        float_format="%.5f",
    )

    t2 = time.perf_counter()

    print(f'Sensitivity ran for {time.strftime("%H:%M:%S", time.gmtime(t2 - t1))}\n')

#%% Bayesian optimization ###########################
if mode_ == "Bayesian":
    if dataset_ == "Jilkova2022":
        path_bayesian = "Bayesian_run_input_2022.csv"
        cols_data_measured = slice(21, 45) #range of the columns to be considered
        cols_data_measured_errors = slice(45, 69) #range of the columns to be considered
        #these extra ones are needed for validation step
        duration = 155 #number of days of incubation
        cols_measured_respSoil = slice(29, 37) #which columns contain measured soil derived respiration
        cols_measured_respSubstrate = slice(37, 45)  #which columns contain measured substrate derived respiration
        cols_measured_PE = slice(69, 77)
        cols_measured_PE_errors = slice(77, 85)
        timepoints = [1, 15, 29, 43, 71, 99, 127, 155] #timepoints for which to calculate the PE effect
    
    if dataset_ == "Jilkova2024":
        path_bayesian = "Bayesian_run_input_2024.csv"
        cols_data_measured = slice(21, 65) #range of the columns to be considered
        cols_data_measured_errors = slice(65, 109) #range of the columns to be considered
        #these extra ones are needed for validation step
        duration = 161 #number of days of incubation
        cols_measured_respSoil = slice(33, 53) #which columns contain measured soil derived respiration
        cols_measured_respSubstrate = slice(66, 77)  #which columns contain measured substrate derived respiration

        

        
    """
    key bayesian principle: the likelihood of a run is the sum of the likelihood of the parameters
    and how good the results fit.

    Bayesian calibration will run the model for a given number of fields (numFields)
    starting with initial parameters as read from json
    output from the runs for all fields is compared to measured data from each field,
    and the combined likelihood of the result and the parameters is calculated
    then a random step is taken changing each of the parameters to calibrate,
    but staying between the min and max value (for now, flat distribution)
    For this new set of parameters all fields are run again,
    then again the combined likelihood calculated. The run is accepted or not (depending on random value),
    if accepted saved into 'posterior', if not simply ignored

    There are 2 input json files: one with the measured data, field data, species data etc.
    The second one with the min/max/initial parameter values of the selected parameters that need to be calibrated

    steps in code below:
       1) calculate the variance of the parameter space: VarianceParameterSpace
       2) calculate the likelihood of the initial parameter set: loglikelihoodParam0
       3) run the model with the initial parameters and store yield and SOC in data_Simulated
       4) calculate the likelihood of each run for each field and store in LikelihoodSimulated
       5) calculate the likelihood of the entire run over all fields in LogLikelihoodSim0
       6) set first run as best fit
     start loop over NumberOfTries
       7) find new candidate parameter values to try
       8) calculate the likelihood of this candidate parameter set: loglikelihoodParam1
       9) run the model with the candidate parameters and store yield and SOC in data_Simulated
      10) calculate the likelihood of each run for each field and store in LikelihoodSimulated
      11) calculate the likelihood of the entire run over all fields in LogLikelihoodSim1
      12) compare the likelihood of this run to the previous and accept the candidate parameters into posterior or not
      13) set parameters as best fit if they are better than current best fit
    """

    # def run_model(inputData, results_path, num_treatments, data=None, parallel=False):

    #     start = time.perf_counter()

    #     for treatment in range(num_treatments):

    #         treatment, result = run_model_bayesian(inputData, results_path, treatment, data)
    #         data_Simulated[treatment] = result

    #     end = time.perf_counter()
    #     print(f'model ran for {time.strftime("%H:%M:%S", time.gmtime(end - start))}')

    # if __name__ != '__main__':
    #     exit(0)

    converged = False
    
    #empty whole output_Bayesian folder 
    for filename in os.listdir(sharable_path):
        file_path = os.path.join(sharable_path, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)  # Delete file or symbolic link
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)  # Delete subdirectory
        

    # Set a folder for logs
    logs_path, run_name = MainFunctionsPotprim.create_log_folder(mode_)

    # Initialize parser
    parser = argparse.ArgumentParser(description="Run Bayesian optimization")
    parser.add_argument(
        "-p",
        "--parallel",
        default=False,
        action="store_true",
        help="Run in parallel mode",
    )
    parser.add_argument(
        "-f",
        "--fields",
        default=0,
        type=int,
        help="Run maximum this number of fields (default: all fields)",
    )

    parser.add_argument(
        "-t",
        "--tries",
        default=3,
        type=int,
        help="Run this number of tries (default: 10000)",
    )  # was 10000
    parser.add_argument(
        "-d", "--debug", default=False, action="store_true", help="Print debug info"
    )

    # Read arguments from command line
    args = parser.parse_args()
    parallel = args.parallel
    debug = args.debug

    #%%--- Set number of tries
    # number of parameter sets to try, including the start, set very high for calibration (10000)
    NumberOfTries = args.tries
    NumberOfTries = 10000
    print("Number of Tries", NumberOfTries)
    t1 = time.perf_counter()

    # results_path = get_results_path()
    # print("Directory ", results_path, " created")

    # numTreatments = len(calibrationData_df['Treatments'])

    # read the Parameter data        
    with open("datalistCalibrationParam.json") as inputCalibrationParamfile:
        (
            numParams,
            CalibParamInit,
            CalParameters,
            CalParameterValues,
            MaximumOption,
            MinimalOption,
            keys,
        ) = BayesianFunctionsPotprim.read_parameter_data(inputCalibrationParamfile)

    # copy initial parameters in output Bayesian
    shutil.copy("datalistCalibrationParam.json",
            os.path.join(sharable_path, "datalistCalibrationParam.json"))
    
    #open the "normal"input 
    with open("datalistInput.json", "r", encoding="utf-8") as inputfileParam:
        AllParam = json.load(inputfileParam)

    # put the initial parametervalues in the correct list so overwrite some parameters        
    AllParam.update(CalibParamInit)

    # # read the measured data (towards which to calibrate) and the treatment definitions
    inputBayesianRun = pd.read_csv(path_bayesian, header=0, skiprows=0)
    numTreatments = len(inputBayesianRun)
    #how many rows contain control treatments (these should come first in the input file), needed for PE calculation
    numControls = inputBayesianRun['treatment'].str.contains('control', case=False, na=False).sum()
    
    # put the variables defining the treatments into 1 list

    # "put the measured data and their errors in separate dataframes
    data_measured = pd.DataFrame()
    data_measured_errors = pd.DataFrame()

     
    data_measured = inputBayesianRun.iloc[:, cols_data_measured]    
    data_measured_errors = inputBayesianRun.iloc[:, cols_data_measured_errors]
    #for PE make it into dictionaries for easier manipulation
    measured_PE = inputBayesianRun.iloc[:, cols_measured_PE].to_dict(orient = "records")  
    measured_PE_errors = inputBayesianRun.iloc[:, cols_measured_PE_errors].to_dict(orient = "records")



    # data_measured_errors.columns
    # inputBayesianRun({"sample"})
    # obtain nonnumeric and numeric part of variable name separately

    data_measured_names = []
    data_measured_days = []
    data_measured_colnames = data_measured.columns.tolist()
    data_measured_names, data_measured_days = (
        BayesianFunctionsPotprim.split_alphanumeric_list(data_measured_colnames)
    )

    # create list of calibrated parameters + only values starting with original
    CalibratedParametersValues = CalParameterValues
    CalibratedParameters = CalParameters

    # create list of all parameter sets tried out
    priorChain = np.zeros([NumberOfTries, numParams])  # list of tries per parameter
    priorChain[0, :] = list(
        CalParameters.values()
    )  # prior is every parameter set you try out, put in first row
    logLseries = []  # create empty list

    # create list of all parameter sets accepted, and 0 when not accepted so linnenrs are equal to priorChain
    posteriorChain = np.zeros([NumberOfTries, numParams])
    # start the chain, will hold all accepted parameter sets, and 0 if not accepted
    posteriorChain[0, :] = list(CalibratedParameters.values())

    #start a count of accepted parameter sets (prior length)
    NumOfAccepted = 0

    # create list of lists for simulations for data on different days
    data_Simulated = dict.fromkeys(
        data_measured, 0
    )  # make a dictionary from column names of measured data
    data_Simulated.update({"sim likelihood": 0})  # add one more element with likelihood, first set to zero
    data_Simulated = [
        # {"resp1": 0, "resp_sub1": 0, "sim likelihood": 0}
        copy.deepcopy(
            data_Simulated
        )  # make deep copy otherwise all point to the same values and later saving will not work
        for treatment in range(numTreatments)  # create a list of dictionaries
    ]

    """
    #%%--- actual start of calibration
    """

    """
    1) calculate the variance of the parameter space
    """
    # variance of the parameter space, is needed to define the step size for each parameter
    VarianceParameterSpace = np.diag(((0.005 * (MaximumOption - MinimalOption)) ** 2))

    """2) calculate the likelihood of the parameters chosen (for a flat distribution this will always be constant or 0)
    """
    # pdf=probability density function, the likelihood of the parameter set
    # loglikelihood_param = np.sum(
    #     np.log(
    #         stats.uniform.pdf(CalibratedParametersValues, MinimalOption, MaximumOption)
    #     )
    # )
    
    loglikelihood_param = np.sum(
        np.log(
            stats.uniform.pdf(CalibratedParametersValues, MinimalOption, MaximumOption - MinimalOption)
        )
    )

    print(loglikelihood_param)
    """
    #%% --- 3) Simulated Data in a similar frame as the measured values, 1 run is over all treatments
    """
    treatmentVar = ()
    results_df = pd.DataFrame()
    
    #make a list of likelihoods to know how many are included in the calculation
    likelihood_list = []
    colnames_list = []
    treatment_list = []
    
    #for each row in input file
    for treatment in range(numTreatments):
        # use input data for the respective treatment
        treatmentVar = inputBayesianRun.iloc[treatment, 0:21]

        results_df = run_model(
            AllParam,
            treatmentVar,
            mode_="Bayesian",
            Plotting=False,
            numDays=duration,
            path=None,
        )

        # we need to couple the output of the right day to the measured output
        # this gets the whole list of variables uploaded as input measured data
        for d in range(len(data_measured_colnames)):  # for each measured variable
            # print(d)
            # print("colname", data_measured_colnames[d], "variable", data_measured_names[d], "day", data_measured_days[d]-1)
            data_Simulated[treatment][data_measured_colnames[d]] = results_df.at[
                data_measured_days[d] - 1, data_measured_names[d]
            ]


        """
        4) calculate the likelihood of each parameter set for each treatment and store in sim likelihood from the differences between measured and simulated and error
        """

        
        for e in range(len(data_measured_colnames)):  
            # if treatment == 3:
                # print(data_measured_colnames[e])
            # for each measured variable calculate loglikelihood
            measurement = data_measured.iat[treatment, e]
            #make sure that measurement is not NA
            if pd.notna(measurement):
                likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                    data_Simulated[treatment][data_measured_colnames[e]],
                    measurement,
                    data_measured_errors.iat[treatment, e],
                )
                data_Simulated[treatment]["sim likelihood"] += likelyhood  # and add it up for all measured variables for the given treatment
                likelihood_list.append(likelyhood)
                colnames_list.append(data_measured_colnames[e])
                treatment_list.append(treatment)
            else: 
                likelyhood = pd.NA
        

        #print("line 878 first likelihood before adding up MB and FB", data_Simulated[treatment]["sim likelihood"])
        #then add likelihood from microbial biomass stability = ratio of initial and final biomass,
        # Access the initial value
        MBini = CalibratedParameters['MBini']
        # print("MBini", MBini)
        MB_simulated = results_df.at[duration - 1, "MB"] #MB on last day
        dMB = MB_simulated/MBini #ratio of final to initial MB       
        
        
        likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
            dMB,
            1, #target value = we want MB to be stable
            0.2, #guesstimated error 20%
        )
        #print("dMB and its likelihood", dMB, likelyhood)
        data_Simulated[treatment]["sim likelihood"] += likelyhood  # and add it up
        likelihood_list.append(likelyhood)
        colnames_list.append("dMB")
        treatment_list.append(treatment)
        
        #same for F:B ratio stability     
        FBini = CalibratedParameters['FBini'] #initial value
        FB_simulated = results_df.at[duration - 1, "FB"] #FB on last day        
        dFB = FB_simulated/FBini #ratio / should be 1        
        likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
            dFB,
            1, #target value, we want FB to be stable
            0.2,#error, guesstimated error 20%
        )
        #print("FBini, FB_simulated, dFB and its likelihood", FBini, FB_simulated, dFB, likelyhood)
        data_Simulated[treatment]["sim likelihood"] += likelyhood  # and add it up
        likelihood_list.append(likelyhood)
        colnames_list.append("dFB")
        treatment_list.append(treatment)


    """
    5) calculate the likelihood of the entire run over all treatments in log_likelihood_sim0
    """
    # first add up likelihoods across all treatments (from comparison of measured and simulated)
    likelihood_simulated = 0
    for treatment in range(numTreatments):  # add up likelihood across treatment
        likelihood_simulated += data_Simulated[treatment]["sim likelihood"]
        # and reset it to zero for the following parameter set trials
        data_Simulated[treatment]["sim likelihood"] = 0
    
    # print("first likelihood_sim0 before adding priming", likelihood_simulated/len(data_Simulated))
    
    #calculate priming effects and their likelihoods
    # Step 1: Calculate averages for control (first 5 entries)
    averages = defaultdict(float)
    counts = defaultdict(int)
    
    # Find all columns starting with respSoil (keys)
    resp_keys = [key for key in data_Simulated[0] if key.startswith('respSoil')]
    
    # Accumulate values from the first 5 entries
    for entry in data_Simulated[:numControls]:
        for key in resp_keys:
            averages[key] += entry.get(key, 0)
            counts[key] += 1
    
    # Calculate final average
    for key in averages:
        averages[key] /= counts[key]
    
    # Step 2: Add PE values to each entry
    for entry in data_Simulated:
        for key in resp_keys:
            pe_key = f"PE{key[8:]}"  # Extract number from 'respSoilX'
            entry[pe_key] = entry.get(key, 0) - averages[key]
   
    #then calculate likelihood connected to priming
    
    #now loop over both modelled and measured data and calculate likelihood of PE
    for idx, (sim, meas, err) in enumerate(zip(data_Simulated, measured_PE, measured_PE_errors)):
        pe_list = []
        for key in meas.keys(): #loop over all measured PE timepoints, e.g.
            err_key = f"{key}_error" #name of the column
            pe_sim = sim.get(key) #get simulated value
            pe_measured = meas.get(key) #get measured value
            pe_measured_error = err.get(err_key) #get error
            likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                pe_sim,
                pe_measured, #target value 
                pe_measured_error, #guesstimated error 0.002 (5% and then averaged)
            )
            # print("simkey, pe_sim, pe_measured, likelyhood", sim_key, pe_sim, pe_measured, likelyhood)
            likelihood_simulated +=  likelyhood  # and add it up
            likelihood_list.append(likelyhood)
            colnames_list.append(key)
            treatment_list.append(idx)
            pe_list.append(pe_sim) #store simulated values in a list
        #calculate total PE for each treatment/run
        pe_tot = sum(pe_list)/len(pe_list) * 24 * duration #priming effect in microgramsC per g of soil over the whole incubation
        sim["PE_tot"] = pe_tot  # Add it to the dictionary
        print("PE", pe_tot)
    
    #make dataframe of likelihoods of individual measured parameters
    likelihoods_overview = pd.DataFrame({
    'treatment': treatment_list,
    'parameter': colnames_list,
    'likelihood': likelihood_list
    })

   
    
    # then use this sum to calculate average likelihood of this parameter set over all treatments and save in log_likelihood_sim0
    log_likelihood_sim0 = likelihood_simulated / len(
        data_measured
    )  # divide by number of treatments
    print("1088 first likelihood_sim0 after including PE", log_likelihood_sim0)
    # print ('priming PE', PE_tot)
    logLseries.append(log_likelihood_sim0)
    """
    6) save best fit, "BestFitParam" is the parameter set giving the highest likelihood (best fit = maximum probability)
    """
    BestFitParam = CalibratedParametersValues  # the initial values of parameters are my best try at first step
    # psetMAP is max fit point, save parameter and likelihood of best run
    log_likelihood_best_fit_param = loglikelihood_param + log_likelihood_sim0

    print("start saving results")

    BayesianFunctionsPotprim.save_result([[data_Simulated]], results_path, "SimdataAll")
    # save_result([[data_Simulated]], results_path, "SimdataBestFit")
    # save_result([[CalibratedParametersValues]], results_path, "calibratedParameters")
    # save_result([[[log_likelihood_sim0]]], results_path, "logLikelyhood")

    print("resultspath", results_path)

    """
    #%% --- loop over number ot tries
    """
    for c in range(0, NumberOfTries):  # For each trial parameter set
        print("Parameter set try:", c+1)
        # 7) find new parameter values to try
        #CalibratedParametersValues,removed replaced 24/6 to test candidateValue
        # if c==0:
        #     candidateValue=CalibratedParametersValues
        candidateparameters, candidateValue, AllParam = (
            BayesianFunctionsPotprim.find_new_parameters(
                CalibratedParametersValues, 
                VarianceParameterSpace,
                MinimalOption,
                MaximumOption,
                keys,
                AllParam,
            )
        )
        print(candidateValue)

        # 8) calculate the likelihood of these new parameters, assuming a uniform distribution
        # pdf=probability density function, the likelihood of the parameter set
        # test = stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption - MinimalOption)
        # test2 = np.log(test)
        loglikelihood_param1 = np.sum(
            np.log(stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption - MinimalOption))
        )

        # a good set has no 0 likelyhood so product is a value but can be negative
        LikelyhoodTest = np.prod(
            stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption - MinimalOption)
        )
        
        # print('loglikelihood_param1',loglikelihood_param1)
        if LikelyhoodTest != 0:
            print("line 537 entered the next parameter set try yaay")
            # if the parameter you want to try is in the range between min and max
            for treatment in range(numTreatments):  # so we run for each treatment
                # print(treatment)
                # 9) run the model for each treatment with the new parameters
                treatmentVar = inputBayesianRun.iloc[
                    treatment, 0:21
                ]  # to be moved & use iloc
                results_df = run_model(
                    AllParam, treatmentVar, mode_, False, numDays=duration, path=None
                )
                # print(results_df)
                # we need to couple the output of the right day to the measured output
                for d in range(len(data_measured_colnames)):  # for each measured variable
                    # print(d)
                    # print("colname", data_measured_colnames[d], "variable", data_measured_names[d], "day", data_measured_days[d]-1)
                    data_Simulated[treatment][data_measured_colnames[d]] = results_df.at[
                        data_measured_days[d] - 1, data_measured_names[d]
                    ]

                                
                # 10) calculate the likelihood of each treatment run for given parameter set
                # we need to add for the treatment the likelyhood of all measurements added
                # if treatment == 2:
                #     print("line 552 safety break ")
                #     break  # safety for now

                for e in range(len(data_measured_colnames)):
                    measurement = data_measured.iat[treatment, e]
                    #make sure that measurement is not NA
                    if pd.notna(measurement):
                        likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                            data_Simulated[treatment][data_measured_colnames[e]],
                            measurement,
                            data_measured_errors.iat[treatment, e],
                        )
                        data_Simulated[treatment]["sim likelihood"] += likelyhood

                #then add likelihood from microbial biomass stability = ratio of initial and final biomass,
                # Access the initial value
                MBini = CalibratedParameters['MBini']
                # print("MBini", MBini)
                MB_simulated = results_df.at[duration - 1, "MB"] #MB on last day
                dMB = MB_simulated/MBini #ratio of final to initial MB           
                
                likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                    dMB,
                    1, #target value = we want MB to be stable
                    0.2, #guesstimated error 10%
                )
                data_Simulated[treatment]["sim likelihood"] += likelyhood  # and add it up
                
                #same for F:B ratio stability     
                FBini = CalibratedParameters['FBini'] #initial value
                FB_simulated = results_df.at[duration - 1, "FB"] #FB on last day        
                dFB = FB_simulated/FBini #ratio / should be 1        
                likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                    dFB,
                    1, #target value, we want FB to be stable
                    0.2,#error, guesstimated
                )
                data_Simulated[treatment]["sim likelihood"] += likelyhood  # and add it up
                
                    # if treatment == 1:
                    #     print(
                    #         "treatment",
                    #         treatment,
                    #         "e",
                    #         e,
                    #         "variable",
                    #         data_measured_colnames[e],
                    #         "simulated",
                    #         data_Simulated[treatment][data_measured_colnames[e]],
                    #         "measured",
                    #         measurement,
                    #         "error",
                    #         data_measured_errors.iat[treatment, e],
                    #         "likelihood",
                    #         likelyhood,
                    #         "overall likelihood",
                    #         data_Simulated[treatment]["sim likelihood"],
                    #     )


            """
            11) calculate the average likelihood of all the treatment runs in LogLikelihoodSim1
            the average over all runs within one try of parameters
            """
            # DiffMeasureSimulated = []  # new empty for every try
            likelihood_simulated = 0  # empty for every try
            for treatment in range(numTreatments):
                # DiffMeasureSimulated.append(data_Simulated[treatment]["sim likelihood"])
                # sum up likelihood across treatments
                likelihood_simulated += data_Simulated[treatment]["sim likelihood"]
                # empty likelihood for next parameter set tries
                data_Simulated[treatment]["sim likelihood"] = 0
            print("line1099 log likelihood_sim1 before adding priming", likelihood_simulated/ len(data_Simulated))
            
            #then calculate likelihood connected to priming
            # Step 1: Calculate averages for control (first 5 entries)
            averages = defaultdict(float)
            counts = defaultdict(int)
            
            # Find all columns starting with respSoil (keys)
            resp_keys = [key for key in data_Simulated[0] if key.startswith('respSoil')]
            
            # Accumulate values from the first 5 entries 
            for entry in data_Simulated[:numControls]:
                for key in resp_keys:
                    averages[key] += entry.get(key, 0)
                    counts[key] += 1
            
            # Calculate final average
            for key in averages:
                averages[key] /= counts[key]
            
            # Step 2: Add PE values to each entry
            for entry in data_Simulated:
                for key in resp_keys:
                    pe_key = f"PE{key[8:]}"  # Extract number from 'respSoilX'
                    entry[pe_key] = entry.get(key, 0) - averages[key]
           
            #then calculate likelihood connected to priming
   
            #now loop over both modelled and measured data and calculate likelihood of PE
            for sim, meas, err in zip(data_Simulated, measured_PE, measured_PE_errors):
                pe_list = []
                for key in meas.keys(): #loop over all measued PE timepoints, e.g.
                    err_key = f"{key}_error" #name of the column
                    pe_sim = sim.get(key) #get simulated value
                    pe_measured = meas.get(key) #get measured value
                    pe_measured_error = err.get(err_key) #get error
                    likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                        pe_sim,
                        pe_measured, #target value 
                        pe_measured_error, #guesstimated error 0.002 (5% and then averaged)
                    )
                    # print("simkey, pe_sim, pe_measured, likelyhood", sim_key, pe_sim, pe_measured, likelyhood)
                    likelihood_simulated +=  likelyhood  # and add it up
                    pe_list.append(pe_sim) #store simulated values in a list
                #calculate total PE for each treatment/run
                pe_tot = sum(pe_list)/len(pe_list) * 24 * duration #priming effect in microgramsC per g of soil over the whole incubation
                sim["PE_tot"] = pe_tot  # Add it to the dictionary
                print("PE", pe_tot)
    
            # divide by number of treatments to obtain average
            log_likelihood_sim1 = likelihood_simulated / len(data_Simulated)
            # print (DiffMeasureSimulated)
            print("line 1325 likelihood_simulated after adding priming", log_likelihood_sim1)
            
            print("finished 11)")

            """
            12) compare the likelihood of this try to the previous and accept into posterior or not
            form the ratio of this step to previous and accept/reject from this (log a/b = log a-log b)
            if the new run fits better it is always accepted (logAlpha>0), 
            if it is worse it is sometimes accepted depending on the random
            """
            print("log_likelihood_sim1", log_likelihood_sim1)
            print("log_likelihood_sim0", log_likelihood_sim0)

            # if treatment == 2:
            #     sys.exit()
            log_likelihood_diff = log_likelihood_sim1 - log_likelihood_sim0
            logalpha= log_likelihood_sim1 - log_likelihood_sim0
            #alpha = np.exp(
            #    log_likelihood_diff
            #)  # if new is better this is bigger than 1

            #alpha = log_likelihood_diff
            random =ra.random()
            lograndom =np.log(random)  # choose random value between 0 and 1
            logLseries.append(log_likelihood_sim1)

            # print('random value', lograndom)
            # print('logalpha', logalpha)
            BayesianFunctionsPotprim.save_result(
                [[data_Simulated]], results_path, "SimdataAll"
            )

            print("random:", random, "logalpha", log_likelihood_diff, "lograndom", lograndom)

            if lograndom < logalpha:
                CalibratedParametersValues = candidateValue
                loglikelihood_param = loglikelihood_param1
                log_likelihood_sim0 = (
                    log_likelihood_sim1  # if accepted move to this point
                )
                posteriorChain[c, :] = (
                    CalibratedParametersValues  # add step to the chain
                )
                NumOfAccepted = NumOfAccepted + 1 #count number of accepted parameter sets (length of posterior)

                """
                13) set parameters as best fit if they are better than current best fit
                test if we have a new best fit
                """

                if (
                    loglikelihood_param + log_likelihood_sim0
                ) > log_likelihood_best_fit_param:
                    log_likelihood_best_fit_param = (
                        loglikelihood_param + log_likelihood_sim0
                    )
                    BestFitParam = (
                        CalibratedParametersValues  # update most likely parameter set
                    )
                    

                    BayesianFunctionsPotprim.save_result(
                        [[data_Simulated]], results_path, "SimdataBestFit"
                    )
                BayesianFunctionsPotprim.save_result(
                    [[CalibratedParametersValues]],
                    sharable_path,
                    "calibratedParameters",
                )
                BayesianFunctionsPotprim.save_result(
                    [[[log_likelihood_sim0]]], sharable_path, "logLikelihood"
                )
                # BayesianFunctionsPotprim.save_result(
                #     [[[log_likelihood_sim0]]], logs_path, "logLikelihood"
                # )
                
                """
                14) test if we have enough runs: avg and stdev are table for each column of posterior
                """
                parameters = pd.read_csv(
                    os.path.join(sharable_path, "calibratedParameters.csv")
                )
               
                #print how many parameter sets were accepted from how many tries
                # numaccepted= len(parameters)
                # print(numaccepted, "parameter sets accepted from ",  c, "number of Tries")
                
                if BayesianFunctionsPotprim.check_dataframe_significant_change(
                    parameters, alpha=0.5, num_identical_results=500
                ):
                    print("Hurraaayyy!!! converged")

                    converged = True

                    break

        # the prior chain saves all tries, also the ones that are not 'saved' in the posterior chain
        priorChain[c, :] = candidateValue
        #print number of accepted versus tried parameters
        print(NumOfAccepted, "parameter sets accepted out of", c+1, "tries")
        #save all tested parameters (prior chain) in a csv also / includes all parameter tries after first try
        BayesianFunctionsPotprim.save_result(
            [[candidateValue]],
            sharable_path,
            "AllTestedParameters",
        )

    """
    end of loop
    """
    #%% --- final export of accepted parameters
    t2 = time.perf_counter()

    # save all accepted parameters set as json >> this is then used for Validation
    AcceptedParamsJsonName = "AcceptedParams_" + run_name
    BayesianFunctionsPotprim.save_json(
        "calibratedParameters.csv",
        "logLikelihood.csv",
        keys,
        sharable_path,
        sharable_path,
        AcceptedParamsJsonName,
    )
    
    #%%--- Histograms of tried and accepted parameters
    #make 2 figures output folders or check if they exist
    sharable_figures = os.path.join(sharable_path, "figures") #output_Bayesian/figures folder
    # logs_figures = os.path.join(logs_path, "figures") #logs/run/figures folder
    
    try:
        os.makedirs(sharable_figures)
    except FileExistsError:
        # directory already exists
        pass
    
    # try:
    #     os.makedirs(logs_figures)
    # except FileExistsError:
    #     # directory already exists
    #     pass
    
    # Load the CSV file into a DataFrame
    df = pd.read_csv("./output_Bayesian/calibratedParameters.csv", header=None)
    df_all = pd.read_csv("./output_Bayesian/AllTestedParameters.csv", header=None)
    # inputCalibrationParamfile = open("datalistCalibrationParam.json")
    # (
    #     numParams,
    #     CalibParamInit,
    #     CalParameters,
    #     CalParameterValues,
    #     MaximumOption,
    #     MinimalOption,
    #     keys,
    # ) = BayesianFunctionsPotprim.read_parameter_data(inputCalibrationParamfile)

    df.columns = keys
    df_all.columns = keys
    
    # Iterate through each column in the DataFrame
    for i, column in enumerate(df.columns):
        # Create a figure for the histograms

        fig, axes = plt.subplots(2, 3, sharey="row", sharex ="all", figsize=(15, 10))
        
        # Flatten the 2D array of axes to a 1D array
        axes = axes.flatten()

        # Get the values of the column for both accepted and all tried parameters
        values = df[column].dropna()  # Drop NaN values if any
        values_all = df_all[column].dropna()  # Drop NaN values if any
        n = len(values)
        
        #Prior
        # Plot histogram for all tried values
        axes[0].hist(values_all, bins=30, color="blue", alpha=0.7)
        axes[0].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[0].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[0].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[0].set_title(f"Prior all Values – {column}")

        # Plot histogram for the last 500 values
        last500 = values_all[-500:]
        axes[1].hist(last500, bins=30, color="blue", alpha=0.7)
        axes[1].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[1].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[1].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[1].set_title(f"Prior last 500 values – {column}")

        # Plot histogram for the previous 500 values
        previous500 = values_all[-1000:-500]
        axes[2].hist(previous500, bins=30, color="blue", alpha=0.7)
        axes[2].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[2].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[2].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[2].set_title(f"Prior previous 500 values – {column}")
        
        #Posterior
        # Plot histogram for all accepted values
        axes[3].hist(values, bins=30, color="blue", alpha=0.7)
        axes[3].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[3].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[3].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[3].set_title(f"Posterior all Values – {column}")

        # Plot histogram for the last 500 values
        last500 = values[-500:]
        axes[4].hist(last500, bins=30, color="blue", alpha=0.7)
        axes[4].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[4].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[4].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[4].set_title(f"Posterior last 500 values – {column}")

        # Plot histogram for the previous 500 values
        previous500 = values[-1000:-500]
        axes[5].hist(previous500, bins=30, color="blue", alpha=0.7)
        axes[5].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[5].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[5].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[5].set_title(f"Posterior previous 500 values – {column}")

        plt.tight_layout()

        # plt.savefig(os.path.join(logs_figures, "hist_" + column + ".png"))
        plt.savefig(os.path.join(sharable_figures, "hist_" + column + ".png")) #save also into output_Bayesian
        plt.close()
        
    #%%Validation part of Bayesian mode ###########################
    #run normal run using best parameter set
    #%% ---  step 1: select the best parameter set
    # Load the "AcceptedParams" json from output_Bayesian/, its name contains also a date, so search for file starting AcceptedParams
    file_path = os.path.join(sharable_path, "AcceptedParams_*")
    files = glob.glob(file_path)
    
    calibParam = [] #empty calibParam after last calibration
    
    
    if len(files) == 1: #safety check, there should be just one file like this
        with open(files[0], 'r') as f:
            calibParam = json.load(f)
        # print(calibParam)
    else:
        raise FileNotFoundError("Expected exactly one file starting with 'AcceptedParams', found: {}".format(len(files)))
    
    
    # Select parameter sets
    selected_sets = []
    selected_likelihoods = []
    bestlikelihood = [] # also empty this
    
    # Select 1 set of parameters with the highest likelihood
    calibParam = {k: calibParam[k] for k in sorted(calibParam)} #sort the calibrated parameter sets by likelihood
    bestlikelihood = list(calibParam.keys())[0] #likelihood of the best parameter set
    selected_sets.append(calibParam[bestlikelihood]) #add to the list, this can be later expanded if multiple parameter sets are chosen for validation
    selected_likelihoods.append(bestlikelihood) #add to the list, this can be later expanded if multiple parameter sets are chosen for validation
    

    # print("all likelihoods", list(calibParam.keys()))
    print("best parameter set likelihood", bestlikelihood)
    print("best parameter set", selected_sets)
    
    # Save the set of parameters that will be used
    with open(os.path.join(sharable_path, "BestParamSetValidation.json"), "w") as json_file:
        json.dump(selected_sets, json_file, indent=4)
    
    #%% --- step 2: combine it with the "fixed" parameters which are not calibrated for
    with open("fixedParameters.json", "r") as f2:
        fixedParam = json.load(
            f2
        )  # only the parameters that are fixed, ie not calibrated
        
    for index, parset in enumerate(selected_sets): #ready for multiple parameter sets
        # Combine the calibrated parameters and the fixed parameters into one variable
        AllParam = {**parset, **fixedParam}
        
        
       #%% --- step 3: run normal run with best parameter set
        Plotting = True # needs to be set to True cuz for Bayesian it is automatically switched to False
        mean_ef, metrics_df, final_results_df = normal_run(path_bayesian, 
                       duration, 
                       cols_measured_respSoil, 
                       cols_measured_respSubstrate,
                       cols_data_measured,
                       cols_measured_PE,
                       timepoints,
                       AllParam,
                       Plotting, 
                       sharable_path
                       )


    #%% --- Final reports of calibration + validation
    # save metadata of the calibration
    content = f"""        
        {NumOfAccepted} parameter sets accepted out of {c+1} tries
        converged = {converged}
        best likelihood = {bestlikelihood}
        mean EF = {mean_ef}
        Calibration ran for {time.strftime('%H:%M:%S', time.gmtime(t2 - t1))}
        {metrics_df}
        """
    print(content)
    with open(os.path.join(sharable_path, "info.txt"), "w") as output_file:
        output_file.write(content)
    
    #copy the whole output_Bayesian folder to logs folder
    shutil.copytree(sharable_path, logs_path, dirs_exist_ok=True)
    
    print(f"Output folder copied to {logs_path}")  
    
    # save metadata in one csv file
    file_exists = os.path.isfile("./logs/logs_Bayesian.csv")
    file_is_empty = file_exists and os.path.getsize("./logs/logs_Bayesian.csv") == 0

    log_likelihood_csv = pd.read_csv(
        os.path.join(sharable_path, "logLikelihood.csv"), header=None
    )
    bestLikelihood = log_likelihood_csv[0].max()
    numberAccepted = log_likelihood_csv.shape[0]

    with open("./logs/logs_Bayesian.csv", "a", newline="") as f:
        writer = csv.writer(f)

        # Write the header
        if not file_exists or file_is_empty:
            writer.writerow(
                [
                    "Run",
                    "Tries",
                    "Converged",
                    "Best likelihood",
                    "Number of accepted",
                    "Mean EF",
                    "Date",
                ]
            )

        writer.writerow(
            [
                run_name,
                NumberOfTries,
                converged,
                bestLikelihood,
                numberAccepted,
                mean_ef,
                datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            ]
        )
#%% Latin Hypercube mode ###########################
if mode_ == "Hypercube":
    #manually set the path of the folder containing the accepted parameter csv and calibration input file
    # logs_path = os.path.join("logs/250809_Bayesian_overall_saved_manually")
    # for i in range(5):
    #     i = i+1
    #     print(i)
    
    #version for kfold validation
    # i = 1
    # foldername =f"set{i}"
    # logs_path = os.path.join("logs/250809_Kfolds_v5_saved_manually",foldername,"output_Bayesian")
    
    #version for scenarios
    logs_path = "logs/251030_Negpriming_availDOMtobulk"
    
   
    # #uncomment if in need to make the hypercube sample again
    # #%% --- step 0.5 make hypercube sample and save them as csv and dictionary    
    # #load csv of accepted params as dataframe
    
    # #path for overall calibration
    # file_path = os.path.join(logs_path, "calibratedParameters.csv") #was created by code stored away in MainFunctionsPotprim in function compile_Bayesians

    
    # # Load the CSV files into a DataFrame    
    # acceptedParams_df = pd.read_csv(file_path, header=None)
    
    # # Number of samples
    # n_samples = 15
    # n_params = acceptedParams_df.shape[1]
    # n_grid = acceptedParams_df.shape[0]  # number of values per parameter
    
    # # Latin Hypercube Sampling in [0,1]
    # sampler = qmc.LatinHypercube(d=n_params)
    # lhs_sample = sampler.random(n=n_samples)
    
    # # Map to grid indices
    # indices = (lhs_sample * n_grid).astype(int)
    # indices = np.clip(indices, 0, n_grid - 1)
    
    # # Sample from DataFrame using the LHS indices
    # sampled_df = pd.DataFrame({
    #     i: acceptedParams_df[i].values[indices[:, i]] for i in range(n_params)
    # })
    
    
    # #save the hypercube sample
    # sampled_df.to_csv(
    #     os.path.join(logs_path, "Hypercube_sampled.csv"),
    #     index=False,
    #     float_format="%.5f",
    # )
    #load hypercube sample from csv
    #version for kfold validation
    # path_hypercube = os.path.join(logs_path, "Hypercube_sampled.csv")
    #version for scenarios, just use the overall hypercube
    # path_hypercube = "Hypercube_sampled.csv"
    path_hypercube = "Negativepriming_availDOMtobulk.csv"

    sampled_df = pd.read_csv(path_hypercube, header=0, skiprows=0)
    #convert the dataframe to a list of dictionaries required further
    
    #to get the names of the parameters that were calibrated for
    # read the calibration Parameter data, obtain keys      
    #version fo kfold
    # with open(os.path.join(logs_path, "datalistCalibrationParam.json")) as inputCalibrationParamfile:
    #version for scenarios
    with open( "datalistCalibrationParam_start_overallJuly28.json") as inputCalibrationParamfile:
        (
            numParams,
            CalibParamInit,
            CalParameters,
            CalParameterValues,
            MaximumOption,
            MinimalOption,
            keys,
        ) = BayesianFunctionsPotprim.read_parameter_data(inputCalibrationParamfile)
        
    
        
    # Convert the dataframe with hypercube sample to list of dictionaries
    dict_list = sampled_df.to_dict(orient='records')
    
    # Manually map keys to each row
    selected_sets = [dict(zip(keys, row)) for row in sampled_df.values]
 
    
    #%% --- step 2: run for all parameter sets from hypercube sample
    #load input
    #version for cross-validation:
    # filename = f"Normal_run_input_2022_subset{i}.csv"
    # path_normal = os.path.join("input files crossvalidation 2022", "validation", filename)
    #version for scenario anylsis
    path_normal = "Normal_run_input_2022baseline_only.csv"
    inputRun = pd.read_csv(path_normal, header=0, skiprows=0)    
    numTreatments = len(inputRun)
    duration = 155 #number of days of incubation
    Plotting = False # switch off for now
    
    df_list1 = [] #list of dataframes, created again for each hypercube sample
    df_list2 = [] #list of all dataframes that are finally joined together
    
    #for each parameter set (this loop happens 15 times)
    for index, parset in enumerate(selected_sets): #ready for multiple parameter sets
       
        # Combine the calibrated parameters and the fixed parameters into one variable
        AllParam.update(parset) #this overwrites part of the input parameters which are calibrated
        print("parameter set no.", index)
        print(AllParam)
        #for each treatment (row in Normal_input) (happens ca 20 times)
        for treatment in range(numTreatments):
            treatmentVar = inputRun.iloc[treatment, 0:21] # select first 21 columns from the input file
            print("treatment", treatment)
            print("main row 153treatmentVar", treatmentVar)
            results_df = run_model(
                AllParam,
                treatmentVar,
                mode_="Normal",
                Plotting=Plotting,
                numDays=duration,
                path=logs_path
            )
            #create treatment ID starting from 1
            results_df['treatmentID'] = treatment + 1
            #create parameter set ID starting from 1
            results_df['setID'] = index + 1
            #store the results in a list of dataframes
            df_list1.append(results_df)         
                            

       #create dataset for all treatments from one parameter set
        final_results_df = pd.concat(
            df_list1, ignore_index=True
        ) 
        #empty list to store the results of the next selected set of hypercube
        df_list1 = [] 
        #store this dataframe in a second list of dataframes
        df_list2.append(final_results_df)
        
   
    #overall output dataset including simulated data from all parameter sets in hypercube sample
    overall_results_df = pd.concat(
        df_list2, ignore_index=True
    ) 
    #%% --- step 3: save output data into Simdata_Hypercube
    #save as csv file
    
    #make output directory (data)
    data_path = os.path.join(logs_path, "data")
    
    try:
            os.makedirs(data_path)
    except FileExistsError:
        # if directory already exists
        pass
    #save all modelled data
    
    overall_results_df.to_csv(
        os.path.join(logs_path, "data/Simdata_Hypercube_validation.csv"),
        index=False,
        float_format="%.5f",
    )
    
# # %%Old Validation run ###########################
# if mode_ == "Validation":
#     # Clear files
#     csv_files = [
#         "selectedLikelihoods.csv",
#         "rmse.csv",
#         "ef.csv",
#     ]

#     for file in csv_files:
#         file_path = os.path.join(sharable_path, file)
#         if os.path.exists(file_path):
#             with open(file_path, "w") as file:
#                 file.write("")  # Clear the contents of the file

#     logs_path, run_name = MainFunctionsPotprim.create_log_folder(mode_)
#     # logs_figures = os.path.join(logs_path, "figures")

#     # try:
#     #     os.makedirs(logs_figures)
#     # except FileExistsError:
#     #     # directory already exists
#     #     pass

#     # Input values
#     inputRun = pd.read_csv("Validation_run_input.csv", header=0, skiprows=0)
#     numTreatments = len(inputRun)

#     # Clear file with mean respirations if it already exists
#     file_path = os.path.join(results_path, "selectedSets_MeanRespiration.csv")
#     if os.path.exists(file_path):
#         with open(file_path, "w") as file:
#             file.write("")  # Clear the contents of the file

#     # Load AcceptedParams
#     pattern = os.path.join(sharable_path, "AcceptedParams_*")
#     matching_files = glob.glob(pattern)

#     # Check if any matching files were found
#     if matching_files:
#         # Open the first matching file (or handle multiple files as needed)
#         with open(matching_files[0], "r") as f1:
#             calibParam = json.load(f1)  # Load the JSON data from the file

#         # Extract the filename without the path
#         filename = os.path.basename(matching_files[0])

#         # Extract the part after "AcceptedParams_"
#         Bayesian_version = "_".join(filename.split("_")[1:])

#     else:
#         print("No AcceptedParams file found.")
#         sys.exit()

#     # save input AcceptedParams to logs
#     with open(os.path.join(logs_path, "AcceptedParams_input.json"), "w") as output_file:
#         json.dump(calibParam, output_file, indent=4)

#     # Select parameter sets
#     selected_sets = []
#     selected_likelihoods = []

#     ##### Latin Hypercube ############
#     # # Extract likelihoods (keys) and parameters sets (values) from uploaded json
#     # likelihoods = []
#     # parameter_sets = []

#     # for likelihood, parameters in calibParam.items():
#     #     likelihoods.append(float(likelihood))
#     #     parameter_sets.append(parameters)

#     # # Calculate weights of the likelihoods
#     # total_likelihood = sum(likelihoods)
#     # weights = [likelihood / total_likelihood for likelihood in likelihoods]

#     # # Actual Latin Hypercube
#     # n_samples = 3  # Adjust based on your needs
#     # sampler = qmc.LatinHypercube(d=len(parameter_sets))
#     # sample = sampler.random(n=n_samples)

#     # for i in range(n_samples):
#     #     index = np.random.choice(len(parameter_sets), p=weights)
#     #     selected_sets.append(parameter_sets[index])
#     #     selected_likelihoods.append(index)

#     # Select 1 set of parameters with the highest likelihood
#     calibParam = {k: calibParam[k] for k in sorted(calibParam)}
#     likelihood = list(calibParam.keys())[0]
#     selected_sets.append(calibParam[likelihood])
#     selected_likelihoods.append(likelihood)

#     # Save the set of parameters that will be used
#     with open(os.path.join(logs_path, "setParamValidation.json"), "w") as json_file:
#         json.dump(selected_sets, json_file, indent=4)

#     # Save the likelihoods
#     file_exists = os.path.isfile(os.path.join(logs_path, "selectedLikelihoods.csv"))
#     file_is_empty = (
#         file_exists
#         and os.path.getsize(os.path.join(logs_path, "selectedLikelihoods.csv")) == 0
#     )
#     with open(
#         os.path.join(logs_path, "selectedLikelihoods.csv"), "w", newline=""
#     ) as file:
#         csv_writer = csv.writer(file)

#         # Write the header
#         if not file_exists or file_is_empty:
#             csv_writer.writerow(["Set", "Likelihood"])

#         for index, likelihood in enumerate(selected_likelihoods):
#             # Write the key and values to the CSV file
#             csv_writer.writerow([index + 1, likelihood])

#     # Merge calibrated parameters with the fixed ones – this should happen inside the for loop in the future
#     with open("fixedParameters.json", "r") as f2:
#         fixedParam = json.load(
#             f2
#         )  # only the parameters that are fixed, ie not calibrated

#     # create lists for respiration plot
#     respSoil_mean_model = []
#     respSubstrate_mean_model = []
#     respSoil_mean_measure = []
#     respSubstrate_mean_measure = []
#     labels = []

#     for index, parset in enumerate(selected_sets):
#         # Combine the calibrated parameters and the fixed parameters into one variable
#         AllParam = {**parset, **fixedParam}

#         # Empty these variables so every plot shows only the values of the specific set
#         # However, it won't be here later on, when we rewrite the Plotting for the mean of the results over sets
#         # This is only provisional
#         respSoil_mean_model = []
#         respSubstrate_mean_model = []
#         respSoil_mean_measure = []
#         respSubstrate_mean_measure = []
#         labels = []

#         final_results_df = pd.DataFrame()  # empty, so every set has its own file
#         df_list = []

#         for treatment in range(numTreatments):
#             treatmentVar = inputRun.iloc[treatment, 0:21]

#             results_df = run_model(
#                 AllParam,
#                 treatmentVar,
#                 mode_="Normal",
#                 Plotting=Plotting,
#                 numDays=161,
#                 path=logs_figures,
#             )
#             df_list.append(results_df)

#             # storing values for respiration plot
#             if Plotting:
#                 labels.append(results_df["treatment"][1])
#                 respSoil_mean_model.append((results_df["respSoil"].mean()) / 0.8 * 24)
#                 respSubstrate_mean_model.append(
#                     (results_df["respSubstrate"].mean()) / 0.8 * 24
#                 )
#                 respSoil_mean_measure.append((inputRun.iloc[treatment, 33:53]).mean())
#                 respSubstrate_mean_measure.append((inputRun.iloc[treatment, 66:77]).mean())
            
            
            
              
                

#         final_results_df = pd.concat(
#             df_list, ignore_index=True
#         )  # add all the rows to the results_df

#         try:
#             os.makedirs("./output/data")
#         except FileExistsError:
#             # directory already exists
#             pass

#         final_results_df.to_csv(
#             os.path.join("./output/data", "Validation_" + str(index + 1) + ".csv"),
#             index=False,
#             float_format="%.5f",
#         )

#         # Save mean respiration in treatments, sets under each other
#         file_path = os.path.join(results_path, "selectedSets_MeanRespiration.csv")
#         file_exists = os.path.isfile(file_path)
#         file_is_empty = file_exists and os.path.getsize(file_path) == 0

#         with open(
#             file_path,
#             mode="a",
#             newline="",
#         ) as csvfile:
#             csv_writer = csv.writer(csvfile)

#             # Write the header
#             if not file_exists or file_is_empty:
#                 csv_writer.writerow(["Set", "Treatment", "respSoil", "respSubstrate"])

#             for label, value1, value2 in zip(
#                 labels, respSoil_mean_model, respSubstrate_mean_model
#             ):
#                 # Write the key and values to the CSV file
#                 csv_writer.writerow([index + 1, label, value1, value2])

#         ###### Plotting is now inside the for loop over selected sets
#         # Later, we should put the Plotting outside the loop and draw it using mean respirations over sets
#         if Plotting:
#             name = "respPlot_Validation_" + str(index + 1) + ".png"
#             MainFunctionsPotprim.drawRespPlot(
#                 labels,
#                 respSoil_mean_model,
#                 respSoil_mean_measure,
#                 respSubstrate_mean_model,
#                 respSubstrate_mean_measure,
#                 name,
#                 logs_figures,
#             )

#     ######## Calculate RMSE ########################
#     # It's now calculated from the last set, needs to be changed for the mean of everything !!!!!!!!!!!

#     # Derive respiration from simulated values (modelled in Validation mode)
#     obs_days_soil = [
#         0,
#         2,
#         6,
#         13,
#         21,
#         23,
#         27,
#         34,
#         49,
#         51,
#         55,
#         62,
#         91,
#         93,
#         97,
#         104,
#         147,
#         149,
#         153,
#         160,
#     ]  # list of days in which the respiration was measured for soil; note that it is 1 smaller than in the input file as in the output file, it starts with 0
#     obs_days_sub = [
#         0,
#         2,
#         6,
#         13,
#         49,
#         51,
#         55,
#         62,
#         147,
#         149,
#         153,
#         160,
#     ]  # list of days in which the respiration was measured for substrate; also starts with 0
#     respSoil_sim = []
#     respSub_sim = []

#     # Iterate over the DataFrame rows
#     for index, row in final_results_df.iterrows():
#         if row["day"] in obs_days_soil:
#             respSoil_sim.append(row["respSoil"])

#         if row["day"] in obs_days_sub:
#             respSub_sim.append(row["respSubstrate"])

#     # Derive respiration from observed values
#     respSoil_obs = []
#     respSub_obs = []

#     for index, row in inputRun.iterrows():
#         for column_name, column_value in row.items():
#             if (
#                 column_name.startswith("resp")
#                 and "sub" not in column_name
#                 and "error" not in column_name
#             ):
#                 if "obs" in column_name:
#                     # Append to respSub_obs
#                     respSub_obs.append(column_value)
#                 else:
#                     # Append to respSoil_obs
#                     respSoil_obs.append(column_value)

#     ######## Calculate RMSE ##################################################
#     RMSE_Soil = calculateRMSE(respSoil_obs, respSoil_sim, "respSoil", logs_path)
#     RMSE_Substrate = calculateRMSE(
#         respSub_obs, respSub_sim, "respSubstrate", logs_path
#     )

#     ######## Calculate EF (Nash-Sutcliffe Efficiency) ########################
#     EF_Soil = calculateEF(respSoil_obs, respSoil_sim, "respSoil", logs_path)
#     EF_Substrate = calculateEF(respSub_obs, respSub_sim, "respSubstrate", logs_path)

#     ######## Calculate Bias ##################################################
#     Bias_Soil = calculateBias(respSoil_obs, respSoil_sim, "respSoil", logs_path)
#     Bias_Substrate = calculateBias(
#         respSub_obs, respSub_sim, "respSubstrate", logs_path
#     )

#     # save metadata in one csv file
#     file_exists = os.path.isfile("./logs/logs_Validation.csv")
#     file_is_empty = file_exists and os.path.getsize("./logs/logs_Validation.csv") == 0

#     # Save overall log file
#     with open("./logs/logs_Validation.csv", "a", newline="") as f:
#         writer = csv.writer(f)

#         # Write the header
#         if not file_exists or file_is_empty:
#             writer.writerow(
#                 [
#                     "Run",
#                     "RMSE_Soil",
#                     "EF_Soil",
#                     "Bias_Soil",
#                     "RMSE_Substrate",
#                     "EF_Substrate",
#                     "Bias_Substrate",
#                     "Bayesian version",
#                     "Date",
#                 ]
#             )

#         writer.writerow(
#             [
#                 run_name,
#                 RMSE_Soil,
#                 EF_Soil,
#                 Bias_Soil,
#                 RMSE_Substrate,
#                 EF_Substrate,
#                 Bias_Substrate,
#                 Bayesian_version,
#                 datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
#             ]
#         )

#     ################ In case we ever need to calculate RMSE for each treatment separately ######################
#     # # Derive respiration from simulated values (modelled in Validation mode)
#     # respSoil_sim = {}
#     # respSub_sim = {}
#     # obs_days_soil = [
#     #     0,
#     #     2,
#     #     6,
#     #     13,
#     #     21,
#     #     23,
#     #     27,
#     #     34,
#     #     49,
#     #     51,
#     #     55,
#     #     62,
#     #     91,
#     #     93,
#     #     97,
#     #     104,
#     #     147,
#     #     149,
#     #     153,
#     #     160,
#     # ]  # list of days in which the respiration was measured for soil; note that it is 1 smaller than in the input file as in the output file, it starts with 0
#     # obs_days_sub = [
#     #     0,
#     #     2,
#     #     6,
#     #     13,
#     #     49,
#     #     51,
#     #     55,
#     #     62,
#     #     147,
#     #     149,
#     #     153,
#     #     160,
#     # ]  # list of days in which the respiration was measured for substrate; also starts with 0

#     # # Iterate over the DataFrame rows
#     # for index, row in final_results_df.iterrows():
#     #     treatment = row["treatment"]

#     #     # Append respSoil values
#     #     if treatment not in respSoil_sim:
#     #         respSoil_sim[treatment] = (
#     #             []
#     #         )  # Create a new list if the treatment is not in the dictionary
#     #     if row["day"] in obs_days_soil:
#     #         respSoil_sim[treatment].append(row["respSoil"])

#     #     # Append respSubstrate values
#     #     if treatment not in respSub_sim:
#     #         respSub_sim[treatment] = (
#     #             []
#     #         )  # Create a new list if the treatment is not in the dictionary
#     #     if row["day"] in obs_days_sub:
#     #         respSub_sim[treatment].append(row["respSubstrate"])

#     # # Derive respiration from observed values
#     # respSoil_obs = {}
#     # respSub_obs = {}

#     # for index, row in inputRun.iterrows():
#     #     treatment = row["treatment"]

#     #     if treatment not in respSoil_obs:
#     #         respSoil_obs[treatment] = []

#     #     if treatment not in respSub_obs:
#     #         respSub_obs[treatment] = []

#     #     for column_name, column_value in row.items():
#     #         if column_name.startswith("resp") and "error" not in column_name:
#     #             if "sub" in column_name:
#     #                 # Append to respSub_obs
#     #                 respSub_obs[treatment].append(column_value)
#     #             else:
#     #                 # Append to respSoil_obs
#     #                 respSoil_obs[treatment].append(column_value)

#     # # Now we will actually calculate the RMSE, yaaay!
#     # rmse_soil = {}
#     # rmse_sub = {}

#     # for (key1, value1), (key2, value2) in zip(
#     #     respSoil_obs.items(), respSoil_sim.items()
#     # ):
#     #     if key1 == key2:  # Ensure the keys match
#     #         actual = np.array(value1)
#     #         predicted = np.array(value2)

#     #         rmse_soil[key1] = np.sqrt(((predicted - actual) ** 2).mean())

#     # for (key1, value1), (key2, value2) in zip(respSub_obs.items(), respSub_sim.items()):
#     #     if key1 == key2:  # Ensure the keys match
#     #         actual = np.array(value1)
#     #         predicted = np.array(value2)

#     #         rmse_sub[key1] = np.sqrt(((predicted - actual) ** 2).mean())

#     # # save as a one csv file
#     # try:
#     #     os.makedirs("./output_Bayesian")
#     # except FileExistsError:
#     #     # directory already exists
#     #     pass

#     # with open(os.path.join(sharable_path, "rmse.csv"), mode="w", newline="") as csvfile:
#     #     csv_writer = csv.writer(csvfile)

#     #     # Write the header
#     #     csv_writer.writerow(["Treatment", "respSoil", "respSub"])

#     #     # Iterate over the keys in the dictionaries
#     #     for key in rmse_soil.keys():
#     #         value1 = rmse_soil[key]
#     #         value2 = rmse_sub[key]
#     #         # Write the key and values to the CSV file
#     #         csv_writer.writerow([key, value1, value2])
