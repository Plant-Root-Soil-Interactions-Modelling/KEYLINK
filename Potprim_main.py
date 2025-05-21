# stuff needed for Bayesian mode
import argparse
# import concurrent.futures
import json

# import math
import os
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import random as ra
from scipy import stats
import BayesianFunctionsPotprim #needed in Bayesian mode
import MainFunctionsPotprim #needed in other modes, plotting, validation etc.
import sys
import csv
from datetime import datetime
import glob
import shutil
from matplotlib.lines import Line2D

# needed by all modes / normal, sensitivity and bayesian mode
from potPrimingMAOMfunction import *
from typing import Literal, get_args

# needed for sensitivity and bayesian
import copy

############## Modes #################################
# set allowed values for mode
modes = Literal["Normal", "Sensitivity", "Bayesian", "Validation"]
options = get_args(modes)

#%% set the mode to Normal, Sensitivity or Bayesian
mode_ = "Normal"

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
inputfileParam = open(
    "datalistInput.json"
)  # input parameters (all) is always in same filenam
AllParam = json.load(inputfileParam)


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
               AllParam, 
               Plotting,
               results_path
               ):
    
        # load input data - Treatments
        inputRun = pd.read_csv(path_normal, header=0, skiprows=0)    
        numTreatments = len(inputRun)

        # create lists for respiration plot
        respSoil_mean_measure = []
        respSubstrate_mean_measure = []
        labels = []
        
        #%%--- run the model
        #run the model for all treatments/rows in treatment input file

        for treatment in range(numTreatments):
            treatmentVar = inputRun.iloc[treatment, 0:21] # select first 21 columns from the input file

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
        #%%--- validation plots
        if Plotting:
            #filter out modelled values for all those variables and days for which we have measured values
            #first automatically extract for which data we have measured data
            data_measured = inputRun.iloc[:, cols_data_measured]  # get measured data
            data_measured_names = []
            data_measured_days = []
            data_measured_colnames = data_measured.columns.tolist() #extract column names of measured data
            
            data_measured_names, data_measured_days = (
                BayesianFunctionsPotprim.split_alphanumeric_list(data_measured_colnames)
            )           # from the column names extract variable name and day of measurement

         
            unique_variables = list(set(data_measured_names)) 
            # print(unique_variables)
            
            # print(final_results_df)
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
            # print(key)              
            #filter the modelled data by the measured variables and days by inner join             
            merged_df = pd.merge(
            df_long, 
            key,
            left_on=['variable', 'day'],
            right_on=['data_measured_names','data_measured_days'],
            how='inner'
            )
            # print(merged_df)
 
            # print(merged_df)
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
            
            # print(data_modelled)
            # print(data_modelled.columns) # has the modelled data but also treatmentID and treatment
            # print(data_measured) 
            # print(data_measured.columns) #has only the actual columns with data, 26 columns
            
            #first plot the non-respiration 1:1 plots
            # -----select columns NOT starting with 'resp' ---
            columns_to_plot = data_measured.columns
            # columns_to_plot = [col for col in data_measured.columns if not col.startswith("resp")]
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
                
            return mean_ef, metrics_df #return mean EF as overall performance metric
        
#%% Normal run with validation #####################################################################
if mode_ == "Normal":
        
    if dataset_ == "Jilkova2022":
        path_normal = "Normal_run_input_2022.csv"
        duration = 155 #number of days of incubation
        cols_measured_respSoil = slice(29, 37) #which columns contain measured soil derived respiration
        cols_measured_respSubstrate = slice(37, 45)  #which columns contain measured substrate derived respiration
        cols_data_measured = slice(21, 45)
        
    if dataset_ == "Jilkova2024":
        path_normal = "Normal_run_input_2024.csv"
        duration = 161 #number of days of incubation
        cols_measured_respSoil = slice(33, 53) #which columns contain measured soil derived respiration
        cols_measured_respSubstrate = slice(66, 77)  #which columns contain measured substrate derived respiration
        cols_data_measured = slice(21, 65) #columns with all measured data
   
                             

# end of normal run function
#run it now
    mean_ef, metrics_df = normal_run(path_normal, 
                   duration, 
                   cols_measured_respSoil, 
                   cols_measured_respSubstrate,
                   cols_data_measured,
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
    
#%% Sensitivity ###########################
if mode_ == "Sensitivity":
    t1 = time.perf_counter()
    # load parameter values from AllParam
    paramsToTestValues = (
        AllParam["bact_DOM_rel"],
        AllParam["DOM_EC"],
        AllParam["kpriming"],
        AllParam["KS"],
        AllParam["KSfungi"],
        AllParam["KSbact"],
        AllParam["kPOM_MAOM"],
        AllParam["kMAOMs_MAOMp"],
        AllParam["MAOMpmaxrate"],
        AllParam["MAOMsmaxrate"],
        AllParam["MAOMratioSP"],
        AllParam["maxEffectBactMAOM"],
        AllParam["maxEffectSA_MAOM"],
        AllParam["maxEffectN_MAOM"],
        AllParam["MM_N_MAOM"],
        AllParam["MM_Bact_MAOM"],
        AllParam["MM_SA_MAOM"],
        AllParam["MM_DOM_MAOM"],
        AllParam["Priming_max"],
    )

    paramsToTestNames = (
        "bact_DOM_rel",
        "DOM_EC",
        "kpriming",
        "KS",
        "KSfungi",
        "KSbact",
        "kPOM_MAOM",
        "kMAOMs_MAOMp",
        "MAOMpmaxrate",
        "MAOMsmaxrate",
        "MAOMratioSP",
        "maxEffectBactMAOM",
        "maxEffectSA_MAOM",
        "maxEffectN_MAOM",
        "MM_N_MAOM",
        "MM_Bact_MAOM",
        "MM_SA_MAOM",
        "MM_DOM_MAOM",
        "Priming_max",
    )

    paramsToTestDict = dict(zip(paramsToTestNames, paramsToTestValues))
    origValues = copy.deepcopy(
        paramsToTestDict
    )  # need deepcopy to not have a pointer but really full copy of values
    paramChanges = np.array([-50, 0, 100])  # % changes to try for each parameter
    numParams = len(paramsToTestValues)
    numValues = len(paramChanges)
    # numRuns_total= len(paramChanges) * len(paramsToTestValues)  # number of sensitivity runs

    # Treatments
    inputRun = pd.read_csv("Bayesian_run_input.csv", header=0, skiprows=0)
    numTreatments = len(inputRun)
    run_info_list = []

    for param in paramsToTestNames:
        # if param == 'DOM_EC':
        #     break
        for paramChange in paramChanges:
            # calculate by how much to change the parameter value, using a relative parameter change
            delta = (
                paramsToTestDict[param] * paramChange / 100
            )  # I change 1 parameter value
            # caculate new value of parameter
            value = paramsToTestDict[param] + delta
            # change the value directly in the parameter set then used by run_model
            AllParam[param] = value

            # for i in range(len(DOMinput_treatments)):
            # numruns = numruns + 1
            # I want to use thevalues from the dict, for sensitivity, so i put all of the values back in the variable (not the fastest way)
            # needs to be changed ifyou change the parameters to test
            # bact_DOM_rel = paramsToTestDict["bact_DOM_rel"]
            # DOM_EC = paramsToTestDict["DOM_EC"]
            # kpriming = paramsToTestDict["kpriming"]
            # KS = paramsToTestDict["KS"]
            # KSfungi = paramsToTestDict["KSfungi"]
            # KSbact = paramsToTestDict["KSbact"]
            # kPOM_MAOM = paramsToTestDict["kPOM_MAOM"]
            # kMAOMs_MAOMp = paramsToTestDict["kMAOMs_MAOMp"]
            # MAOMpmaxrate = paramsToTestDict["MAOMpmaxrate"]
            # MAOMsmaxrate = paramsToTestDict["MAOMsmaxrate"]
            # MAOMratioSP = paramsToTestDict["MAOMratioSP"]
            # maxEffectBactMAOM = paramsToTestDict["maxEffectBactMAOM"]
            # maxEffectSA_MAOM = paramsToTestDict["maxEffectSA_MAOM"]
            # maxEffectN_MAOM = paramsToTestDict["maxEffectN_MAOM"]
            # MM_N_MAOM = paramsToTestDict["MM_N_MAOM"]
            # MM_Bact_MAOM = paramsToTestDict["MM_Bact_MAOM"]
            # MM_SA_MAOM = paramsToTestDict["MM_SA_MAOM"]
            # MM_DOM_MAOM = paramsToTestDict["MM_DOM_MAOM"]
            # Priming_max = paramsToTestDict["Priming_max"]

            for treatment in range(numTreatments):
                treatmentVar = inputRun.iloc[treatment, 0:21]

                results_df = run_model(
                    AllParam,
                    treatmentVar,
                    mode_="Sensitivity",
                    Plotting=Plotting,
                    numDays=161,
                    path=None,
                )
                # add to the simulated values information about the parameter, its change and value
                # temp_df_list = [param, paramChange, value] + temp_df_list
                # print("temp_df_list", temp_df_list)
                length = len(results_df)

                info_df = pd.DataFrame(
                    {
                        "param": np.full(length, param),
                        "paramChange": np.full(length, paramChange),
                        "value": np.full(length, value),
                    }
                )

                results_df = info_df.join(results_df)
                df_list.append(
                    results_df
                )  # append doesn't work for dataframes, so the dataframes have to be appended to a list to later use concat
                # save parameter set for each run
                All_param_series = pd.Series(AllParam)
                info_series = pd.Series(
                    [param, paramChange, value], index=["param", "paramChange", "value"]
                )
                # Concatenate the three Series

                run_info = pd.concat([info_series, treatmentVar, All_param_series])
                run_info_list.append(run_info)

            # after all runs with one parameter set, reset parameters to original, before next parameter value change
            # paramsToTestDict = copy.deepcopy(origValues)
        # after all changes tried for certain parameter, reset its value to original value
        AllParam[param] = paramsToTestDict[param]

    final_results_df = pd.concat(
        df_list, ignore_index=True
    )  # add all the rows to the results_df
    run_info_df = pd.concat(run_info_list, ignore_index=True)

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
    #     print("Start running model")
    #     start = time.perf_counter()

    #     for treatment in range(num_treatments):
    #         print('treatment =', treatment)
    #         treatment, result = run_model_bayesian(inputData, results_path, treatment, data)
    #         data_Simulated[treatment] = result

    #     end = time.perf_counter()
    #     print(f'model ran for {time.strftime("%H:%M:%S", time.gmtime(end - start))}')

    # if __name__ != '__main__':
    #     exit(0)

    converged = False

    # clear the csv files so it won't append after the existing values from the run before
    csv_files = [
        "calibratedParameters.csv",
        "AllTestedParameters.csv",
        "logLikelihood.csv",
    ]

    csv_files2 = ["SimdataAll.csv", "SimdataBestFit.csv"]

    for file in csv_files:
        file_path = os.path.join(sharable_path, file)
        if os.path.exists(file_path):
            with open(file_path, "w") as file:
                file.write("")  # Clear the contents of the file

    for file in csv_files2:
        file_path = os.path.join(results_path, file)
        if os.path.exists(file_path):
            with open(file_path, "w") as file:
                file.write("")  # Clear the contents of the file


    # Delete previous AcceptedParams file in output_Bayesian folder
    pattern = os.path.join(sharable_path, "AcceptedParams_*")
    matching_files = glob.glob(pattern)

    # Check if any matching files were found
    try:
        os.remove(matching_files[0])
    except:
        pass

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
    #%%--- Set number of tries
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

    # number of parameter sets to try, including the start, set very high for calibration (10000)
    NumberOfTries = args.tries
    NumberOfTries = 10000
    print("Number of Tries", NumberOfTries)
    t1 = time.perf_counter()

    # results_path = get_results_path()
    # print("Directory ", results_path, " created")

    # numTreatments = len(calibrationData_df['Treatments'])

    # read the Parameter data to  and the fixed parameter values and put them together
    inputCalibrationParamfile = open("datalistCalibrationParam.json")
    (
        numParams,
        CalibParamInit,
        CalParameters,
        CalParameterValues,
        MaximumOption,
        MinimalOption,
        keys,
    ) = BayesianFunctionsPotprim.read_parameter_data(inputCalibrationParamfile)

    # Save initial parameters to output Bayesian
    with open("datalistCalibrationParam.json", "r") as inputCalibrationParamfile:
        data = json.load(inputCalibrationParamfile)

    with open(
        os.path.join(sharable_path, "datalistCalibrationParam.json"), "w"
    ) as output_file:
        json.dump(data, output_file, indent=4)

    # put the initial parametervalues in the correct list so overwrite some parameters
    # you can start your 'walk' from another point then the old parameter value
    AllParam.update(CalibParamInit)

    # # read the measured data (towards which to calibrate) and the treatment definitions
    inputBayesianRun = pd.read_csv(path_bayesian, header=0, skiprows=0)
    numTreatments = len(inputBayesianRun)

    # put the variables defining the treatments into 1 list

    # "put the measured data and their errors in separate dataframes
    data_measured = pd.DataFrame()
    data_measured_errors = pd.DataFrame()
     
    data_measured = inputBayesianRun.iloc[:, cols_data_measured]    
    data_measured_errors = inputBayesianRun.iloc[:, cols_data_measured_errors]


    # data_measured_errors.columns
    # inputBayesianRun({"sample"})
    # obtain nonnumeric and numeric part of variable name separately

    data_measured_names = []
    data_measured_days = []
    data_measured_colnames = data_measured.columns.tolist()
    # print(data_measured_colnames)
    #data currently used for calibration Jilkova 2022:
        # 'POM155', 
        # 'MAOM155',
        # 'POM_sub155', 
        # 'MAOM_sub155',
        # 'bact155', 
        # 'fungi155',
        # 'bact_sub155',
        # 'fungi_sub155', 
        # 'respSoil1', 'respSoil15', 'respSoil29', 'respSoil43', 'respSoil71', 'respSoil99', 'respSoil127', 'respSoil155', 
        # 'respSubstrate1', 'respSubstrate15', 'respSubstrate29', 'respSubstrate43', 'respSubstrate71', 'respSubstrate99', 'respSubstrate127', 'respSubstrate155'
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
    # print(parameterlist_df)
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
    loglikelihood_param = np.sum(
        np.log(
            stats.uniform.pdf(CalibratedParametersValues, MinimalOption, MaximumOption)
        )
    )
    
    # loglikelihood_param = np.sum(
    #     np.log(
    #         stats.uniform.pdf(CalibratedParametersValues, MinimalOption, MaximumOption - MinimalOption)
    #     )
    # )

    print(loglikelihood_param)
    """
    #%% --- 3) Simulated Data in a similar frame as the measured values, 1 run is over all treatments
    """
    treatmentVar = ()
    results_df = pd.DataFrame()

    for treatment in range(numTreatments):
        # use input data for the respective treatment
        treatmentVar = inputBayesianRun.iloc[treatment, 0:21]

        # to be  corrected for nr of columns needed
        # print(
        #     "treatmentID",
        #     treatmentVar["treatmentID"],
        #     "treatment",
        #     treatmentVar["treatment"],
        # )

        results_df = run_model(
            AllParam,
            treatmentVar,
            mode_="Bayesian",
            Plotting=False,
            numDays=161,
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

        # we need to add for the treatment the likelyhood of all measurements added, data_measured is df so other indexing
            
    
        # print(treatment)
        # if treatment == 2:
        #     print("line 447 safety break ")
        #     break  # safety for now

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
            else: 
                likelyhood = pd.NA
                
            # if treatment == 3:
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

            # if treatment == 3:
            #     sys.exit()
            # if treatment == 3:
            #     break
            # sys.exit('safety stop to run only for first 2 treatments')

    """
    5) calculate the likelihood of the entire run over all treatments in log_likelihood_sim0
    """
    # first add up likelihoods across all treatments
    likelihood_simulated = 0
    for treatment in range(numTreatments):  # add up likelihood across treatment
        likelihood_simulated += data_Simulated[treatment]["sim likelihood"]
        # and reset it to zero for the following parameter set trials
        data_Simulated[treatment]["sim likelihood"] = 0

    # then use this sum to calculate average likelihood of this parameter set over all treatments and save in log_likelihood_sim0
    log_likelihood_sim0 = likelihood_simulated / len(
        data_measured
    )  # divide by number of treatments
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
    loop over number ot tries
    """
    for c in range(0, NumberOfTries):  # For each trial parameter set
        print("Parameter set try:", c+1)
        # 7) find new parameter values to try

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
        test = stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption)
        test2 = np.log(test)
        loglikelihood_param1 = np.sum(
            np.log(stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption))
        )
        # a good set has no 0 likelyhood so product is a value but can be negative
        LikelyhoodTest = np.product(
            stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption)
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
                    AllParam, treatmentVar, mode_, False, numDays=161, path=None
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

            # log_likelihood_sim1 = sum(DiffMeasureSimulated) / len(data_Simulated)
            # divide by number of treatments to obtain average
            log_likelihood_sim1 = likelihood_simulated / len(data_Simulated)
            # print (DiffMeasureSimulated)

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

    try:
        os.remove(os.path.join(sharable_path, "AcceptedParams_*"))
    except Exception as e:
        pass

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
        # fig, axes = plt.subplots(1, 3, figsize=(15, 5))
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
    
    if len(files) == 1: #safety check, there should be just one file like this
        with open(files[0], 'r') as f:
            calibParam = json.load(f)
        # print(calibParam)
    else:
        raise FileNotFoundError("Expected exactly one file starting with 'AcceptedParams', found: {}".format(len(files)))
    
    
    # Select parameter sets
    selected_sets = []
    selected_likelihoods = []
    
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
        mean_ef, metrics_df = normal_run(path_bayesian, 
                       duration, 
                       cols_measured_respSoil, 
                       cols_measured_respSubstrate,
                       cols_data_measured,
                       AllParam,
                       Plotting, 
                       sharable_path
                       )
    #copy the whole output_Bayesian folder to logs folder
    shutil.copytree(sharable_path, logs_path, dirs_exist_ok=True)
    
    print(f"Output folder copied to {logs_path}")  

    #%% --- Final reports of calibration + validation
    # save metadata of the calibration
    content = f"""\        
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
