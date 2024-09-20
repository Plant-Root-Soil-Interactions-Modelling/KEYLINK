#stuff needed for Bayesian mode
import argparse
# import concurrent.futures
import json
# import math
import os
import time

import numpy as np
import pandas as pd
from numpy import random as ra
from scipy import stats
import BayesianFunctionsPotprim
import sys

#needed by all modes / normal, sensitivity and bayesian mode
from potPrimingMAOMfunction import *
from typing import Literal, get_args

#needed for sensitivity and bayesian
import copy


the_current_path = os.path.abspath(os.getcwd())
results_path = ".\output"

############## Modes #################################
#set allowed values for mode
modes = Literal["Normal", "Sensitivity", "Bayesian", "Jilkova2022"]
options = get_args(modes)

#set the mode to Normal, Sensitivity or Bayesian
mode_ = "Bayesian" #'Jilkova2022'
#check if mode was set correctly, if not stop the run
assert mode_ in options, f'"{mode_}" is not in "{options}"'

#if the Normal or Jilkova2022 mode was chosen, you can decide to turn on the Plotting
Plotting = False

############## Creating the Respiration Plot #########
def drawRespPlot (labels, respSoil_mean_model, respSoil_mean_measure, respSubstrate_mean_model, respSubstrate_mean_measure):
    # create plot
    plt.figure(figsize=(10,12))
    x = np.arange(len(labels))  # label locations
    width = 0.3  # width of the bars

    # create first subplot
    plt.subplot(2, 1, 1)
    plt.bar(x - width / 2, respSoil_mean_model, width, label="Modeled", color="gray")
    plt.bar(
        x + width / 2, respSoil_mean_measure, width, label="Measured", color="black"
    )
    plt.title("Soil derived")
    plt.ylabel("Respiration [µg C-CO2/g soil/h]")
    plt.xticks(x, labels, rotation=90, ha='center')
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1), shadow=True)

    # create second subplot
    plt.subplot(2, 1, 2)
    plt.bar(x - width / 2, respSubstrate_mean_model, width, label="Modeled", color="gray")
    plt.bar(
        x + width / 2, respSubstrate_mean_measure, width, label="Measured", color="black"
    )
    plt.title("Substrate derived")
    plt.ylabel("Respiration [µg C-CO2/g soil/h]")
    plt.xticks(x, labels, rotation=90, ha='center')
    

    plt.tight_layout()

    plt.savefig("./output/figures/respPlot.png")
    plt.close()
        

############## Read data #############################
# read the fixed parameter list
inputfileParam = open(
    "datalistInput.json"
)  # input parameters (all) is always in same filenam
AllParam = json.load(inputfileParam)



# Create variables
df_list = []
temp_df = []  # temporary df to store returned dataframe


############## Normal run ###########################
if mode_ == "Normal":
    # Treatments
    inputRun = pd.read_csv(
        "Normal_run_input.csv", header=0, skiprows=0
    )  
    
    numTreatments = len(inputRun)

    # create lists for respiration plot
    respSoil_mean_model = []
    respSubstrate_mean_model = []
    respSoil_mean_measure = []
    respSubstrate_mean_measure = []
    labels = []
    
    for treatment in range(numTreatments):
        treatmentVar = inputRun.iloc[treatment, 0:17]
        
        temp_df = run_model(
            AllParam, treatmentVar, mode_ = "Normal", Plotting=Plotting
        )
        df_list.append(temp_df)

        # storing values for respiration plot
        if Plotting:
            labels.append(temp_df["treatment"][1])
            respSoil_mean_model.append((temp_df["respSoil"].mean()) / 0.8 * 24)
            respSubstrate_mean_model.append((temp_df["respSubstrate"].mean()) / 0.8 * 24)
            respSoil_mean_measure.append((inputRun.iloc[treatment, 17:37]).mean())
            respSubstrate_mean_measure.append((inputRun.iloc[treatment, 37:48]).mean())

    results_df = pd.concat(
        df_list, ignore_index=True
    )  # add all the rows to the results_df
    
    try:
        os.makedirs("./output/data")
    except FileExistsError:
        # directory already exists
        pass
    
    results_df.to_csv(
        "./output/data/Normal.csv",
        index=False,
        float_format="%.2f",
    )

    if Plotting:
        drawRespPlot(labels, respSoil_mean_model, respSoil_mean_measure, respSubstrate_mean_model, respSubstrate_mean_measure)

############## Jilkova2022 run ###########################
if mode_ == "Jilkova2022":
    # Treatments
    inputRun = pd.read_csv(
        "Normal_run_input_2022.csv", header=0, skiprows=0
    )  
    
    numTreatments = len(inputRun)

    # create lists for respiration plot
    respSoil_mean_model = []
    respSubstrate_mean_model = []
    respSoil_mean_measure = []
    respSubstrate_mean_measure = []
    labels = []
    
    for treatment in range(numTreatments):
        treatmentVar = inputRun.iloc[treatment, 0:17]
        
        temp_df = run_model(
            AllParam, treatmentVar, mode_ = "Normal", Plotting=Plotting
        )
        df_list.append(temp_df)

        # storing values for respiration plot
        if Plotting:
            labels.append(temp_df["treatment"][1])
            respSoil_mean_model.append((temp_df["respSoil"].mean()) / 0.8 * 24)
            respSubstrate_mean_model.append((temp_df["respSubstrate"].mean()) / 0.8 * 24)

            if temp_df["treatment"][1] == "control":
                respSoil_mean_measure.append(0.3199)
                respSubstrate_mean_measure.append(0)
            if temp_df["treatment"][1] == "leachates":
                respSoil_mean_measure.append(0.3980)
                respSubstrate_mean_measure.append(0.0833)
            if temp_df["treatment"][1] == "exudates":
                respSoil_mean_measure.append(0.3371)
                respSubstrate_mean_measure.append(0.1028)
                

    results_df = pd.concat(
        df_list, ignore_index=True
    )  # add all the rows to the results_df
    
    try:
        os.makedirs("./output/data")
    except FileExistsError:
        # directory already exists
        pass
    
    results_df.to_csv(
        "./output/data/Normal2022.csv",
        index=False,
        float_format="%.2f",
    )

    if Plotting:
        drawRespPlot(labels, respSoil_mean_model, respSoil_mean_measure, respSubstrate_mean_model, respSubstrate_mean_measure)
        
############## Sensitivity ###########################
if mode_ == 'Sensitivity':
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
    inputRun = pd.read_csv(
        "Bayesian_run_input.csv", header=0, skiprows=0
    )  
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
            #change the value directly in the parameter set then used by run_model
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
                treatmentVar = inputRun.iloc[treatment, 0:17]
                
                temp_df = run_model(
                    AllParam, treatmentVar, mode_ = "Sensitivity", Plotting=False
                )
                #add to the simulated values information about the parameter, its change and value
                # temp_df_list = [param, paramChange, value] + temp_df_list    
                # print("temp_df_list", temp_df_list)
                length = len(temp_df)
                
                info_df = pd.DataFrame({
                    'param': np.full(length, param),
                    'paramChange': np.full(length, paramChange),
                    'value': np.full(length, value)
                })
                
                temp_df = info_df.join(temp_df)
                df_list.append(
                    temp_df
                )  # append doesn't work for dataframes, so the dataframes have to be appended to a list to later use concat
                #save parameter set for each run
                All_param_series = pd.Series(AllParam)
                info_series = pd.Series([param, paramChange, value], index=['param', 'paramChange', 'value'])
                # Concatenate the three Series
                
                run_info = pd.concat([info_series, treatmentVar, All_param_series])
                run_info_list.append(run_info)
                
            # after all runs with one parameter set, reset parameters to original, before next parameter value change
            # paramsToTestDict = copy.deepcopy(origValues)
        #after all changes tried for certain parameter, reset its value to original value
        AllParam[param] = paramsToTestDict[param] 
     
    results_df = pd.concat(
        df_list, ignore_index=True
    )  # add all the rows to the results_df
    run_info_df = pd.concat(
        run_info_list, ignore_index=True
    ) 
    
    try:
        os.makedirs("./output/data")
    except FileExistsError:
        # directory already exists
        pass
    
    results_df.to_csv(
        ".\output\data\Sensitivity.csv",
        index=False,
        float_format="%.2f",
    )
    
    run_info_df.to_csv(
        ".\output\data\Sensitivity_runs.csv",
        index=False,
        float_format="%.2f",
    )

    t2 = time.perf_counter()
    
    print(f'Sensitivity ran for {time.strftime("%H:%M:%S", time.gmtime(t2 - t1))}\n')
    
############## Bayesian optimization ###########################
if mode_ == 'Bayesian':
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
    
    
    # Initialize parser
    parser = argparse.ArgumentParser(description="Run Bayesian optimization")
    parser.add_argument(
        "-p", "--parallel", default=False, action="store_true", help="Run in parallel mode"
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
        default=10000,
        type=int,
        help="Run this number of tries (default: 10000)",
    ) #was 10000
    parser.add_argument(
        "-d", "--debug", default=False, action="store_true", help="Print debug info"
    )
    
    # Read arguments from command line
    args = parser.parse_args()
    parallel = args.parallel
    debug = args.debug
    
    # number of parameter sets to try, including the start, set very high for calibration (10000)
    NumberOfTries = args.tries 
    print('Number of Tries', NumberOfTries)
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
    
    # put the initial parametervalues in the correct list so overwrite some parameters
    # you can start your 'walk' from another point then the old parameter value
    AllParam.update(CalibParamInit)
    
    # # read the measured data (towards which to calibrate) and the treatment definitions
    inputBayesianRun = pd.read_csv("Bayesian_run_input.csv", header=0, skiprows=0)
    numTreatments = len(inputBayesianRun)
    
    # put the variables defining the treatments into 1 list
    
    
    # "put the measured data and their errors in separate dataframes
    data_measured = pd.DataFrame()
    data_measured_errors = pd.DataFrame()
    
    data_measured = inputBayesianRun.iloc[:, 17:49]    
    data_measured_errors = inputBayesianRun.iloc[:, 49:81]
    # data_measured_errors.columns
    # inputBayesianRun({"sample"})
    #obtain nonnumeric and numeric part of variable name separately
    
    data_measured_names = []
    data_measured_days = []
    data_measured_colnames = data_measured.columns.tolist()
    data_measured_names, data_measured_days = BayesianFunctionsPotprim.split_alphanumeric_list(data_measured_colnames)
    
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
    
    # create list of lists for simulations for data on different days
    data_Simulated = dict.fromkeys(data_measured, 0) #make a dictionary from column names of measured data
    data_Simulated.update({"sim likelihood": 0}) #add one more element with likelihood
    data_Simulated = [
        # {"resp1": 0, "resp_sub1": 0, "sim likelihood": 0}
        copy.deepcopy(data_Simulated) #make deep copy otherwise all point to the same values and later saving will not work
        for treatment in range(numTreatments) #create a list of dictionaries
    ]
 
    
    """
     actual start of calibration
    """
    
    # 1) calculate the variance of the parameter space
    # variance of the parameter space, is needed to define the step size for each parameter
    VarianceParameterSpace = np.diag(((0.005 * (MaximumOption - MinimalOption)) ** 2))
    
    # 2) calculate the likelihood of the parameters chosen (for a flat distribution this will always be constant or 0)
    # pdf=probability density function, the likelihood of the parameter set
    loglikelihood_param = np.sum(
        np.log(stats.uniform.pdf(CalibratedParametersValues, MinimalOption, MaximumOption))
    )
    # 3) Simulated Data in a similar frame as the measured values, 1 run is over all treatments
    treatmentVar = ()
    results_df = pd.DataFrame()
    
    # data_Simulated=pd.DataFrame()
    
    for treatment in range(numTreatments):
        #use input data for the respective treatment
        treatmentVar = inputBayesianRun.iloc[
            treatment, 0:17
        ] 
        print('treatmentVar',treatmentVar)
        # to be  corrected for nr of columns needed
        # print(
        #     "treatmentID",
        #     treatmentVar["treatmentID"],
        #     "treatment",
        #     treatmentVar["treatment"],
        # )

        results_df = run_model(
            AllParam, treatmentVar, mode_ = "Bayesian", Plotting=False
        )
        # df_list.append(
        #     temp_df_list
        # )  # append doesn't work for dataframes, so the lists have to be appended to later use concat
        
        # we need to couple the output of the right day to the measured output
        #this gets the whole list of variables uploaded as input measured data
        for d in range(len(data_measured_colnames)): #for each measured variable
            # print(d)
            # print("colname", data_measured_colnames[d], "variable", data_measured_names[d], "day", data_measured_days[d]-1)
            data_Simulated[treatment][data_measured_colnames[d]] = results_df.at[data_measured_days[d]-1,data_measured_names[d]]
            
            #old version
            # data_Simulated[treatment]["resp1"] = results_df.at[0,'resp']
            # data_Simulated[treatment]["resp_sub1"] = results_df.at[0, 'resp_sub']
        # we need to add for the treatment the likelyhood of all measurements added, data_measured is df so other indexing
        
        # print("datasim resp1 treatment 1", data_Simulated[treatment]["resp1"])
        # print("datameasured", data_measured["resp1"][treatment])
        
        if treatment == 3: 
            print("line 447 safety break ")
            break #safety for now
            
        for e in range(len(data_measured_colnames)):
        
            likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                data_Simulated[treatment][data_measured_colnames[e]],
                data_measured.iat[treatment, e],
                data_measured_errors.iat[treatment, e],
            )
            data_Simulated[treatment]["sim likelihood"] += likelyhood
            
            print("treatment", treatment,
                  "e", e, 
                  "variable", data_measured_colnames[e],
                  'simulated', data_Simulated[treatment][data_measured_colnames[e]],
                  'measured', data_measured.iat[treatment, e],
                  'error', data_measured_errors.iat[treatment, e],
                  'likelihood', likelyhood,    
                  "overall likelihood", data_Simulated[treatment]["sim likelihood"]                  
                  )
            if treatment == 3: 
                break
                # sys.exit('safety stop to run only for first 2 treatments')
            
        #old version only for one variable resp1
        # likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
        #     data_Simulated[treatment]["resp1"],
        #     data_measured["resp1"][treatment],
        #     data_measured_errors["resp1_error"][treatment],
        # )
        # data_Simulated[treatment]["sim likelihood"] += likelyhood
    
    
      
    
    # 4) calculate the likelihood of each run for each treatment from the differences between measured and simulated and error
    # create empty list(logLi) and store all the differences between measured and simulated per day
    likelihood_simulated = 0
    for treatment in range(
        numTreatments
    ):  # add up likelihood across treatment
        likelihood_simulated += data_Simulated[treatment]["sim likelihood"]
        #and reset it to zero for the following parameter set trials
        data_Simulated[treatment]["sim likelihood"] = 0
    
    
    # 5) calculate average total likelihood of this parameter set over all treatments (Log0)
    log_likelihood_sim0 = likelihood_simulated / len(data_measured)
    logLseries.append(log_likelihood_sim0)
    
    #  6) save best fit, "BestFitParam" is the parameter set giving the highest likelihood (best fit = maximum probability)
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
        # print("c", c)
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
        test=stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption)
        test2=np.log(test)
        loglikelihood_param1 = np.sum(
            np.log(stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption))
        )
        # a good set has no 0 likelyhood so product is a value but can be negative
        LikelyhoodTest=np.product(stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption))
        # print('loglikelihood_param1',loglikelihood_param1)
        if (
            LikelyhoodTest != 0
        ):  
            print('line 537 entered the next parameter set try yaay')
            # if the parameter you want to try is in the range between min and max
            for treatment in range(
                numTreatments
            ):  # so we run for each treatment
                # print(treatment)
                # 9) run the model for each treatment with the new parameters
                treatmentVar = inputBayesianRun.iloc[treatment, 0:17]  # to be moved & use iloc
                results_df = run_model(AllParam, treatmentVar, mode_, False)
                # print(results_df)
                # we need to couple the output of the right day to the measured output
                data_Simulated[treatment]["resp1"] = results_df.at[0,'resp']
                data_Simulated[treatment]["resp_sub1"] = results_df.at[0, 'resp_sub']
                
                # 10) calculate the likelihood of each treatment run for given parameter set
                # we need to add for the treatment the likelyhood of all measurements added 
                if treatment == 2: 
                    print("line 552 safety break ")
                    break #safety for now
                
                for e in range(len(data_measured_colnames)):
                
                    likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                        data_Simulated[treatment][data_measured_colnames[e]],
                        data_measured.iat[treatment, e],
                        data_measured_errors.iat[treatment, e],
                    )
                    data_Simulated[treatment]["sim likelihood"] += likelyhood
                    
                
                
                #old version
                # likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                #     data_Simulated[treatment]["resp1"],
                #     data_measured["resp1"][treatment],
                #     data_measured["resp1_error"][treatment],
                # )
                # data_Simulated[treatment]["sim likelihood"] = likelyhood
    
            
    

            """
            11) calculate the average likelihood of all the treatment runs in LogLikelihoodSim1
            the average over all runs within one try of parameters
            """
            # DiffMeasureSimulated = []  # new empty for every try
            likelihood_simulated = 0 #empty for every try
            for treatment in range(numTreatments):
                # DiffMeasureSimulated.append(data_Simulated[treatment]["sim likelihood"])
                #sum up likelihood across treatments
                likelihood_simulated += data_Simulated[treatment]["sim likelihood"]
                #empty likelihood for next parameter set tries
                data_Simulated[treatment]["sim likelihood"] = 0
            
            # log_likelihood_sim1 = sum(DiffMeasureSimulated) / len(data_Simulated)
            #divide by number of treatments to obtain average
            log_likelihood_sim1 = likelihood_simulated / len(data_Simulated)
            # print (DiffMeasureSimulated)
    
            """
            12) compare the likelihood of this try to the previous and accept into posterior or not
            form the ratio of this step to previous and accept/reject from this (log a/b = log a-log b)
            if the new run fits better it is always accepted (logAlpha>0), 
            if it is worse it is sometimes accepted depending on the random
            """
            print('log_likelihood_sim1', log_likelihood_sim1)
            alpha = log_likelihood_sim1 / log_likelihood_sim0 # if new is better this is bigger than 1
            random = ra.random()  # choose random value between 0 and 1
            logLseries.append(log_likelihood_sim1)
    
            # print('random value', lograndom)
            # print('logalpha', logalpha)
            BayesianFunctionsPotprim.save_result(
                [[data_Simulated]], results_path, "SimdataAll"
            )
    
            if random < alpha:
                CalibratedParametersValues = candidateValue
                loglikelihood_param = loglikelihood_param1
                log_likelihood_sim0 = log_likelihood_sim1  # if accepted move to this point
                posteriorChain[c, :] = CalibratedParametersValues  # add step to the chain
    
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
                    [[CalibratedParametersValues]], results_path, "calibratedParameters"
                )
                BayesianFunctionsPotprim.save_result(
                    [[[log_likelihood_sim0]]], results_path, "logLikelihood   "
                )
    
                """
                14) test if we have enough runs: avg and stdev are table for each column of posterior
                """
    
                parameters = pd.read_csv(
                    os.path.join(results_path, "calibratedParameters.csv")
                )
                if BayesianFunctionsPotprim.check_dataframe_significant_change(
                    parameters, alpha=0.5, num_identical_results=50
                ):
                    print("Hurraaayyy!!! converged")
                    break
    
        # the prior chain saves all tries, also the ones that are not 'saved' in the posterior chain
        priorChain[c, :] = candidateValue
    
    """
    end of loop
    """
    
    t2 = time.perf_counter()
    
    print(f'Calibration ran for {time.strftime("%H:%M:%S", time.gmtime(t2 - t1))}\n')
    # df = pd.DataFrame(priorChain, columns=['a0Photo_Eff', 'SoilRootCond', 'minimalStomatalResistance', 'a9DistriRoot',
    #                                        'a10DistriFruit', 'a91DistriHyphae', 'ratioExudRoot', 'f_MIC_DOM',
    #                                        'f_FOM_fPOM', 'f_fPOM_oPOM', 'MAOMsmaxrate', 'MAOMpmaxrate', 'AgRate',
    #                                        'MIC_gmax', 'Mic_RespRate', 'MIC_TR', 'HyphalExploration',
    #                                        'HyphaeTurnoverrate', 'Hyphae_fGRSP', 'SoilHyphaeCond'])
    df = pd.DataFrame(priorChain, columns=list(CalibratedParameters.keys()))
    
    # bayesian_plots(df=df, path=file_name, columns=5, save_to_file=True)
