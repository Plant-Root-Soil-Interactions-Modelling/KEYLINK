# stuff needed for Bayesian mode
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
from scipy.stats import qmc
import BayesianFunctionsPotprim
import sys
import csv

# needed by all modes / normal, sensitivity and bayesian mode
from potPrimingMAOMfunction import *
from typing import Literal, get_args

# needed for sensitivity and bayesian
import copy


the_current_path = os.path.abspath(os.getcwd())
results_path = "./output"
sharable_path = "./output_Bayesian"

try:
    os.makedirs("./output_Bayesian")
except FileExistsError:
    # directory already exists
    pass

############## Modes #################################
# set allowed values for mode
modes = Literal["Normal", "Sensitivity", "Bayesian", "Jilkova2022", "Validation"]
options = get_args(modes)

# set the mode to Normal, Sensitivity or Bayesian
mode_ = "Bayesian"
# check if mode was set correctly, if not stop the run
assert mode_ in options, f'"{mode_}" is not in "{options}"'

# if the Normal or Jilkova2022 mode was chosen, you can decide to turn on the Plotting
Plotting = True

# safety condition / do not allow plotting with Sensitivity or Bayesian mode
if mode_ == "Sensitivity" or mode_ == "Bayesian":
    Plotting = False


############## Creating the Respiration Plot #########
def drawRespPlot(
    labels,
    respSoil_mean_model,
    respSoil_mean_measure,
    respSubstrate_mean_model,
    respSubstrate_mean_measure,
    name,
):
    # create plot
    plt.figure(figsize=(10, 12))
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
    plt.xticks(x, labels, rotation=90, ha="center")
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1), shadow=True)

    # create second subplot
    plt.subplot(2, 1, 2)
    plt.bar(
        x - width / 2, respSubstrate_mean_model, width, label="Modeled", color="gray"
    )
    plt.bar(
        x + width / 2,
        respSubstrate_mean_measure,
        width,
        label="Measured",
        color="black",
    )
    plt.title("Substrate derived")
    plt.ylabel("Respiration [µg C-CO2/g soil/h]")
    plt.xticks(x, labels, rotation=90, ha="center")

    plt.tight_layout()

    if "Validation" in name:
        try:
            os.makedirs("./output_Bayesian/figures")
        except FileExistsError:
            # directory already exists
            pass

        plt.savefig(os.path.join("./output_Bayesian/figures/", name))

    else:
        plt.savefig(os.path.join("./output/figures/", name))

    plt.close()


############## RMSE ##################################
def calculateRMSE(actual, predicted, variable, filepath):
    actual = np.array(actual)
    predicted = np.array(predicted)

    rmse = np.sqrt(((predicted - actual) ** 2).mean())

    # save as a one csv file
    try:
        os.makedirs(filepath)
    except FileExistsError:
        # directory already exists
        pass

    file_exists = os.path.isfile(os.path.join(filepath, "rmse.csv"))
    file_is_empty = (
        file_exists and os.path.getsize(os.path.join(filepath, "rmse.csv")) == 0
    )

    with open(os.path.join(filepath, "rmse.csv"), newline="", mode="a") as file:
        csv_writer = csv.writer(file)

        # Write the header
        if not file_exists or file_is_empty:
            csv_writer.writerow(["Variable", "RMSE"])

        # Write the key and values to the CSV file
        csv_writer.writerow([variable, rmse])


############## EF (Nash-Sutcliffe Efficiency) ########
def calculateEF(actual, predicted, variable, filepath):
    actual = np.array(actual)
    predicted = np.array(predicted)

    ef = 1 - (
        np.sum((actual - predicted) ** 2) / np.sum((actual - np.mean(actual)) ** 2)
    )

    # save as a one csv file
    try:
        os.makedirs(filepath)
    except FileExistsError:
        # directory already exists
        pass

    file_exists = os.path.isfile(os.path.join(filepath, "ef.csv"))
    file_is_empty = (
        file_exists and os.path.getsize(os.path.join(filepath, "ef.csv")) == 0
    )

    with open(os.path.join(filepath, "ef.csv"), newline="", mode="a") as file:
        csv_writer = csv.writer(file)

        # Write the header
        if not file_exists or file_is_empty:
            csv_writer.writerow(["Variable", "EF"])

        # Write the key and values to the CSV file
        csv_writer.writerow([variable, ef])


############## Read data #############################
# read the fixed parameter list
inputfileParam = open(
    "datalistInput.json"
)  # input parameters (all) is always in same filenam
AllParam = json.load(inputfileParam)


# Create variables
df_list = []
results_df = []  # temporary df to store returned dataframe


############## Normal run ###########################
if mode_ == "Normal":
    # Treatments
    inputRun = pd.read_csv("Normal_run_input.csv", header=0, skiprows=0)

    numTreatments = len(inputRun)

    # create lists for respiration plot
    respSoil_mean_model = []
    respSubstrate_mean_model = []
    respSoil_mean_measure = []
    respSubstrate_mean_measure = []
    labels = []

    for treatment in range(numTreatments):
        treatmentVar = inputRun.iloc[treatment, 0:17]

        results_df = run_model(
            AllParam, treatmentVar, mode_="Normal", Plotting=Plotting, numDays=161
        )
        df_list.append(results_df)

        # storing values for respiration plot
        if Plotting:
            labels.append(results_df["treatment"][1])
            respSoil_mean_model.append((results_df["respSoil"].mean()) / 0.8 * 24)
            respSubstrate_mean_model.append(
                (results_df["respSubstrate"].mean()) / 0.8 * 24
            )
            respSoil_mean_measure.append((inputRun.iloc[treatment, 17:37]).mean())
            respSubstrate_mean_measure.append((inputRun.iloc[treatment, 50:61]).mean())

    final_results_df = pd.concat(
        df_list, ignore_index=True
    )  # add all the rows to the results_df

    try:
        os.makedirs("./output/data")
    except FileExistsError:
        # directory already exists
        pass

    final_results_df.to_csv(
        "./output/data/Normal.csv",
        index=False,
        float_format="%.5f",
    )

    if Plotting:
        drawRespPlot(
            labels,
            respSoil_mean_model,
            respSoil_mean_measure,
            respSubstrate_mean_model,
            respSubstrate_mean_measure,
            "respPlot_Normal.png",
        )

############## Jilkova2022 run ###########################
if mode_ == "Jilkova2022":
    # Treatments
    inputRun = pd.read_csv("Normal_run_input_2022.csv", header=0, skiprows=0)

    numTreatments = len(inputRun)

    # create lists for respiration plot
    respSoil_mean_model = []
    respSubstrate_mean_model = []
    respSoil_mean_measure = []
    respSubstrate_mean_measure = []
    labels = []

    for treatment in range(numTreatments):
        treatmentVar = inputRun.iloc[treatment, 0:17]

        results_df = run_model(
            AllParam, treatmentVar, mode_="Normal", Plotting=Plotting, numDays=155
        )
        df_list.append(results_df)

        # storing values for respiration plot
        if Plotting:
            labels.append(results_df["treatment"][1])
            respSoil_mean_model.append((results_df["respSoil"].mean()) / 0.8 * 24)
            respSubstrate_mean_model.append(
                (results_df["respSubstrate"].mean()) / 0.8 * 24
            )

            if results_df["treatment"][1] == "control":
                respSoil_mean_measure.append(0.3199)
                respSubstrate_mean_measure.append(0)
            if results_df["treatment"][1] == "leachates":
                respSoil_mean_measure.append(0.3980)
                respSubstrate_mean_measure.append(0.0833)
            if results_df["treatment"][1] == "exudates":
                respSoil_mean_measure.append(0.3371)
                respSubstrate_mean_measure.append(0.1028)

    final_results_df = pd.concat(
        df_list, ignore_index=True
    )  # add all the rows to the results_df

    try:
        os.makedirs("./output/data")
    except FileExistsError:
        # directory already exists
        pass

    final_results_df.to_csv(
        "./output/data/Normal2022.csv",
        index=False,
        float_format="%.5f",
    )

    if Plotting:
        drawRespPlot(
            labels,
            respSoil_mean_model,
            respSoil_mean_measure,
            respSubstrate_mean_model,
            respSubstrate_mean_measure,
            "respPlot_Jilkova2022.png",
        )

############## Sensitivity ###########################
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
                treatmentVar = inputRun.iloc[treatment, 0:17]

                results_df = run_model(
                    AllParam,
                    treatmentVar,
                    mode_="Sensitivity",
                    Plotting=Plotting,
                    numDays=161,
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

############## Bayesian optimization ###########################
if mode_ == "Bayesian":
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

    # clear the csv files so it won't append after the existing values from the run before
    csv_files = [
        "calibratedParameters.csv",
        "logLikelihood.csv",
        "SimdataAll.csv",
        "SimdataBestFit.csv",
        "BestFitParams.csv",
    ]

    for file in csv_files:
        file_path = os.path.join(results_path, file)
        if os.path.exists(file_path):
            with open(file_path, "w") as file:
                file.write("")  # Clear the contents of the file

    with open(os.path.join(results_path, "BestFitParams.json"), "w") as file:
        file.write("")

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
        default=100000,
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
    # print(parameterlist_df)

    # create list of lists for simulations for data on different days
    data_Simulated = dict.fromkeys(
        data_measured, 0
    )  # make a dictionary from column names of measured data
    data_Simulated.update({"sim likelihood": 0})  # add one more element with likelihood
    data_Simulated = [
        # {"resp1": 0, "resp_sub1": 0, "sim likelihood": 0}
        copy.deepcopy(
            data_Simulated
        )  # make deep copy otherwise all point to the same values and later saving will not work
        for treatment in range(numTreatments)  # create a list of dictionaries
    ]

    """
    actual start of calibration
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
    """
    # 3) Simulated Data in a similar frame as the measured values, 1 run is over all treatments
    """
    treatmentVar = ()
    results_df = pd.DataFrame()

    # data_Simulated=pd.DataFrame()

    for treatment in range(numTreatments):
        # use input data for the respective treatment
        treatmentVar = inputBayesianRun.iloc[treatment, 0:17]
        print("treatmentVar", treatmentVar)
        # to be  corrected for nr of columns needed
        # print(
        #     "treatmentID",
        #     treatmentVar["treatmentID"],
        #     "treatment",
        #     treatmentVar["treatment"],
        # )

        results_df = run_model(
            AllParam, treatmentVar, mode_="Bayesian", Plotting=False, numDays=161
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

        # if treatment == 3:
        #     print("line 447 safety break ")
        #     break  # safety for now

        """
        4) calculate the likelihood of each parameter set for each treatment and store in sim likelihood from the differences between measured and simulated and error
        """
        for e in range(
            len(data_measured_colnames)
        ):  # for each measured variable calculate loglikelihood

            likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                data_Simulated[treatment][data_measured_colnames[e]],
                data_measured.iat[treatment, e],
                data_measured_errors.iat[treatment, e],
            )
            data_Simulated[treatment][
                "sim likelihood"
            ] += likelyhood  # and add it up for all measured variables for the given treatment

            print(
                "treatment",
                treatment,
                "e",
                e,
                "variable",
                data_measured_colnames[e],
                "simulated",
                data_Simulated[treatment][data_measured_colnames[e]],
                "measured",
                data_measured.iat[treatment, e],
                "error",
                data_measured_errors.iat[treatment, e],
                "likelihood",
                likelyhood,
                "overall likelihood",
                data_Simulated[treatment]["sim likelihood"],
            )
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
        print(c)
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
                    treatment, 0:17
                ]  # to be moved & use iloc
                results_df = run_model(
                    AllParam, treatmentVar, mode_, False, numDays=161
                )
                # print(results_df)
                # we need to couple the output of the right day to the measured output
                data_Simulated[treatment]["resp1"] = results_df.at[0, "resp"]
                data_Simulated[treatment]["resp_sub1"] = results_df.at[0, "resp_sub"]

                # 10) calculate the likelihood of each treatment run for given parameter set
                # we need to add for the treatment the likelyhood of all measurements added
                # if treatment == 2:
                #     print("line 552 safety break ")
                #     break  # safety for now

                for e in range(len(data_measured_colnames)):

                    likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                        data_Simulated[treatment][data_measured_colnames[e]],
                        data_measured.iat[treatment, e],
                        data_measured_errors.iat[treatment, e],
                    )
                    data_Simulated[treatment]["sim likelihood"] += likelyhood

                # old version
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
            alpha = (
                log_likelihood_sim1 / log_likelihood_sim0
            )  # if new is better this is bigger than 1

            random = ra.random()  # choose random value between 0 and 1
            logLseries.append(log_likelihood_sim1)

            # print('random value', lograndom)
            # print('logalpha', logalpha)
            BayesianFunctionsPotprim.save_result(
                [[data_Simulated]], results_path, "SimdataAll"
            )

            print("random:", random, "alpha:", alpha)

            if random < alpha:
                CalibratedParametersValues = candidateValue
                loglikelihood_param = loglikelihood_param1
                log_likelihood_sim0 = (
                    log_likelihood_sim1  # if accepted move to this point
                )
                posteriorChain[c, :] = (
                    CalibratedParametersValues  # add step to the chain
                )

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

                """
                14) test if we have enough runs: avg and stdev are table for each column of posterior
                """
                parameters = pd.read_csv(
                    os.path.join(sharable_path, "calibratedParameters.csv")
                )

                if BayesianFunctionsPotprim.check_dataframe_significant_change(
                    parameters, alpha=0.5, num_identical_results=500
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

    # save all accepted parameters set as json
    BayesianFunctionsPotprim.save_json(
        "calibratedParameters.csv",
        "logLikelihood.csv",
        keys,
        sharable_path,
        "BestFitParams",
    )

    ################## Histograms of accepted parameters
    # if mode_ == "Histogram":
    # Load the CSV file into a DataFrame
    df = pd.read_csv("./output_Bayesian/calibratedParameters.csv", header=None)

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

    # Iterate through each column in the DataFrame
    for i, column in enumerate(df.columns):
        # Create a figure for the histograms
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # Get the values of the column
        values = df[column].dropna()  # Drop NaN values if any
        n = len(values)

        # Plot histogram for all values
        axes[0].hist(values, bins=30, color="blue", alpha=0.7)
        axes[0].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[0].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[0].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[0].set_title(f"All Values – {column}")

        # Plot histogram for the last 500 values
        last500 = values[-500:]
        axes[1].hist(last500, bins=30, color="blue", alpha=0.7)
        axes[1].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[1].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[1].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[1].set_title(f"Last 500 values – {column}")

        # Plot histogram for the previous 500 values
        previous500 = values[-1000:-500]
        axes[2].hist(previous500, bins=30, color="blue", alpha=0.7)
        axes[2].axvline(MinimalOption[i], color="black", linestyle="--", label="max")
        axes[2].axvline(MaximumOption[i], color="black", linestyle="--", label="min")
        axes[2].axvline(
            CalParameterValues[i], color="red", linestyle="--", label="initial"
        )
        axes[2].set_title(f"Previous 500 values – {column}")

        plt.tight_layout()

        try:
            os.makedirs("./output_Bayesian/figures")
        except FileExistsError:
            # directory already exists
            pass
        plt.savefig(
            os.path.join("./output_Bayesian/figures/", "hist_" + column + ".png")
        )
        plt.close()


############## Validation run ###########################
if mode_ == "Validation":
    # Clear files
    csv_files = [
        "selectedLikelihoods.csv",
        "rmse.csv",
        "ef.csv",
    ]

    for file in csv_files:
        file_path = os.path.join(sharable_path, file)
        if os.path.exists(file_path):
            with open(file_path, "w") as file:
                file.write("")  # Clear the contents of the file

    # Input values
    inputRun = pd.read_csv("Validation_run_input.csv", header=0, skiprows=0)
    numTreatments = len(inputRun)

    # Clear file with mean respirations if it already exists
    file_path = os.path.join(results_path, "selectedSets_MeanRespiration.csv")
    if os.path.exists(file_path):
        with open(file_path, "w") as file:
            file.write("")  # Clear the contents of the file

    # Calibrated Parameters
    with open("./output_Bayesian/BestFitParams.json", "r") as f1:
        calibParam = json.load(
            f1
        )  # only the parameters that were accepted, all of them

    ##### Latin Hypercube ############
    # Extract likelihoods (keys) and parameters sets (values) from uploaded json
    likelihoods = []
    parameter_sets = []

    for likelihood, parameters in calibParam.items():
        likelihoods.append(float(likelihood))
        parameter_sets.append(parameters)

    # Calculate weights of the likelihoods
    total_likelihood = sum(likelihoods)
    weights = [likelihood / total_likelihood for likelihood in likelihoods]

    # Actual Latin Hypercube
    n_samples = 3  # Adjust based on your needs
    sampler = qmc.LatinHypercube(d=len(parameter_sets))
    sample = sampler.random(n=n_samples)

    selected_sets = []
    selected_likelihoods = []
    for i in range(n_samples):
        index = np.random.choice(len(parameter_sets), p=weights)
        selected_sets.append(parameter_sets[index])
        selected_likelihoods.append(index)

    # Save the set of parameters that will be used
    with open(os.path.join(sharable_path, "setParamValidation.json"), "w") as json_file:
        json.dump(selected_sets, json_file, indent=4)

    # Save the likelihoods
    file_exists = os.path.isfile(os.path.join(sharable_path, "selectedLikelihoods.csv"))
    file_is_empty = (
        file_exists
        and os.path.getsize(os.path.join(sharable_path, "selectedLikelihoods.csv")) == 0
    )
    with open(
        os.path.join(sharable_path, "selectedLikelihoods.csv"), "w", newline=""
    ) as file:
        csv_writer = csv.writer(file)

        # Write the header
        if not file_exists or file_is_empty:
            csv_writer.writerow(["Set", "Variable"])

        for index, likelihood in enumerate(selected_likelihoods):
            # Write the key and values to the CSV file
            csv_writer.writerow([index + 1, likelihood])

    ############# That's all for Latin Hypercube ########################

    # # Select 1 set of parameters with the highest likelihood
    # calibParam = {k: calibParam[k] for k in sorted(calibParam)}
    # likelihood = list(calibParam.keys())[-1]
    # setCalibParam = calibParam[likelihood]

    # Merge calibrated parameters with the fixed ones – this should happen inside the for loop in the future
    with open("fixedParameters.json", "r") as f2:
        fixedParam = json.load(
            f2
        )  # only the parameters that are fixed, ie not calibrated

    # create lists for respiration plot
    respSoil_mean_model = []
    respSubstrate_mean_model = []
    respSoil_mean_measure = []
    respSubstrate_mean_measure = []
    labels = []

    for index, set in enumerate(selected_sets):
        # Combine the calibrated parameters and the fixed parameters into one variable
        AllParam = {**set, **fixedParam}

        # Empty these variables so every plot shows only the values of the specific set
        # However, it won't be here later on, when we rewrite the Plotting for the mean of the results over sets
        # This is only provisional
        respSoil_mean_model = []
        respSubstrate_mean_model = []
        respSoil_mean_measure = []
        respSubstrate_mean_measure = []
        labels = []

        final_results_df = pd.DataFrame()  # empty, so every set has its own file
        df_list = []

        for treatment in range(numTreatments):
            treatmentVar = inputRun.iloc[treatment, 0:17]

            results_df = run_model(
                AllParam, treatmentVar, mode_="Normal", Plotting=Plotting, numDays=161
            )
            df_list.append(results_df)

            # storing values for respiration plot
            if Plotting:
                respSoil_mean_measure.append((inputRun.iloc[treatment, 17:37]).mean())
                respSubstrate_mean_measure.append(
                    (inputRun.iloc[treatment, 50:61]).mean()
                )

            # This is now needed for the mean respiration output with all the treatments and sets
            respSoil_mean_model.append((results_df["respSoil"].mean()) / 0.8 * 24)
            respSubstrate_mean_model.append(
                (results_df["respSubstrate"].mean()) / 0.8 * 24
            )
            labels.append(results_df["treatment"][1])

        final_results_df = pd.concat(
            df_list, ignore_index=True
        )  # add all the rows to the results_df

        try:
            os.makedirs("./output/data")
        except FileExistsError:
            # directory already exists
            pass

        final_results_df.to_csv(
            os.path.join("./output/data", "Validation_" + str(index + 1) + ".csv"),
            index=False,
            float_format="%.5f",
        )

        # Save mean respiration in treatments, sets under each other
        file_path = os.path.join(results_path, "selectedSets_MeanRespiration.csv")
        file_exists = os.path.isfile(file_path)
        file_is_empty = file_exists and os.path.getsize(file_path) == 0

        with open(
            file_path,
            mode="a",
            newline="",
        ) as csvfile:
            csv_writer = csv.writer(csvfile)

            # Write the header
            if not file_exists or file_is_empty:
                csv_writer.writerow(["Set", "Treatment", "respSoil", "respSubstrate"])

            for label, value1, value2 in zip(
                labels, respSoil_mean_model, respSubstrate_mean_model
            ):
                # Write the key and values to the CSV file
                csv_writer.writerow([index + 1, label, value1, value2])

        ###### Plotting is now inside the for loop over selected sets
        # Later, we should put the Plotting outside the loop and draw it using mean respirations over sets
        if Plotting:
            name = "respPlot_Validation_" + str(index + 1) + ".png"
            drawRespPlot(
                labels,
                respSoil_mean_model,
                respSoil_mean_measure,
                respSubstrate_mean_model,
                respSubstrate_mean_measure,
                name,
            )

    ######## Calculate RMSE ########################
    # It's now calculated from the last set, needs to be changed for the mean of everything !!!!!!!!!!!

    # Derive respiration from simulated values (modelled in Validation mode)
    obs_days_soil = [
        0,
        2,
        6,
        13,
        21,
        23,
        27,
        34,
        49,
        51,
        55,
        62,
        91,
        93,
        97,
        104,
        147,
        149,
        153,
        160,
    ]  # list of days in which the respiration was measured for soil; note that it is 1 smaller than in the input file as in the output file, it starts with 0
    obs_days_sub = [
        0,
        2,
        6,
        13,
        49,
        51,
        55,
        62,
        147,
        149,
        153,
        160,
    ]  # list of days in which the respiration was measured for substrate; also starts with 0
    respSoil_sim = []
    respSub_sim = []

    # Iterate over the DataFrame rows
    for index, row in final_results_df.iterrows():
        if row["day"] in obs_days_soil:
            respSoil_sim.append(row["respSoil"])

        if row["day"] in obs_days_sub:
            respSub_sim.append(row["respSubstrate"])

    # Derive respiration from observed values
    respSoil_obs = []
    respSub_obs = []

    for index, row in inputRun.iterrows():
        for column_name, column_value in row.items():
            if (
                column_name.startswith("resp")
                and "sub" not in column_name
                and "error" not in column_name
            ):
                if "obs" in column_name:
                    # Append to respSub_obs
                    respSub_obs.append(column_value)
                else:
                    # Append to respSoil_obs
                    respSoil_obs.append(column_value)

    # Now we will actually calculate the RMSE, yaaay!
    calculateRMSE(respSoil_obs, respSoil_sim, "respSoil", sharable_path)
    calculateRMSE(respSub_obs, respSub_sim, "respSubstrate", sharable_path)

    ######## Calculate EF (Nash-Sutcliffe Efficiency) ########################
    calculateEF(respSoil_obs, respSoil_sim, "respSoil", sharable_path)
    calculateEF(respSub_obs, respSub_sim, "respSubstrate", sharable_path)

    ################ In case we ever need to calculate RMSE for each treatment separately ######################
    # # Derive respiration from simulated values (modelled in Validation mode)
    # respSoil_sim = {}
    # respSub_sim = {}
    # obs_days_soil = [
    #     0,
    #     2,
    #     6,
    #     13,
    #     21,
    #     23,
    #     27,
    #     34,
    #     49,
    #     51,
    #     55,
    #     62,
    #     91,
    #     93,
    #     97,
    #     104,
    #     147,
    #     149,
    #     153,
    #     160,
    # ]  # list of days in which the respiration was measured for soil; note that it is 1 smaller than in the input file as in the output file, it starts with 0
    # obs_days_sub = [
    #     0,
    #     2,
    #     6,
    #     13,
    #     49,
    #     51,
    #     55,
    #     62,
    #     147,
    #     149,
    #     153,
    #     160,
    # ]  # list of days in which the respiration was measured for substrate; also starts with 0

    # # Iterate over the DataFrame rows
    # for index, row in final_results_df.iterrows():
    #     treatment = row["treatment"]

    #     # Append respSoil values
    #     if treatment not in respSoil_sim:
    #         respSoil_sim[treatment] = (
    #             []
    #         )  # Create a new list if the treatment is not in the dictionary
    #     if row["day"] in obs_days_soil:
    #         respSoil_sim[treatment].append(row["respSoil"])

    #     # Append respSubstrate values
    #     if treatment not in respSub_sim:
    #         respSub_sim[treatment] = (
    #             []
    #         )  # Create a new list if the treatment is not in the dictionary
    #     if row["day"] in obs_days_sub:
    #         respSub_sim[treatment].append(row["respSubstrate"])

    # # Derive respiration from observed values
    # respSoil_obs = {}
    # respSub_obs = {}

    # for index, row in inputRun.iterrows():
    #     treatment = row["treatment"]

    #     if treatment not in respSoil_obs:
    #         respSoil_obs[treatment] = []

    #     if treatment not in respSub_obs:
    #         respSub_obs[treatment] = []

    #     for column_name, column_value in row.items():
    #         if column_name.startswith("resp") and "error" not in column_name:
    #             if "sub" in column_name:
    #                 # Append to respSub_obs
    #                 respSub_obs[treatment].append(column_value)
    #             else:
    #                 # Append to respSoil_obs
    #                 respSoil_obs[treatment].append(column_value)

    # # Now we will actually calculate the RMSE, yaaay!
    # rmse_soil = {}
    # rmse_sub = {}

    # for (key1, value1), (key2, value2) in zip(
    #     respSoil_obs.items(), respSoil_sim.items()
    # ):
    #     if key1 == key2:  # Ensure the keys match
    #         actual = np.array(value1)
    #         predicted = np.array(value2)

    #         rmse_soil[key1] = np.sqrt(((predicted - actual) ** 2).mean())

    # for (key1, value1), (key2, value2) in zip(respSub_obs.items(), respSub_sim.items()):
    #     if key1 == key2:  # Ensure the keys match
    #         actual = np.array(value1)
    #         predicted = np.array(value2)

    #         rmse_sub[key1] = np.sqrt(((predicted - actual) ** 2).mean())

    # # save as a one csv file
    # try:
    #     os.makedirs("./output_Bayesian")
    # except FileExistsError:
    #     # directory already exists
    #     pass

    # with open(os.path.join(sharable_path, "rmse.csv"), mode="w", newline="") as csvfile:
    #     csv_writer = csv.writer(csvfile)

    #     # Write the header
    #     csv_writer.writerow(["Treatment", "respSoil", "respSub"])

    #     # Iterate over the keys in the dictionaries
    #     for key in rmse_soil.keys():
    #         value1 = rmse_soil[key]
    #         value2 = rmse_sub[key]
    #         # Write the key and values to the CSV file
    #         csv_writer.writerow([key, value1, value2])
