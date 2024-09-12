import argparse
import concurrent.futures
import json
import math
import os
import time

import numpy as np
import pandas as pd
from numpy import random as ra
from scipy import stats
from potPrimingMAOMfunction import *

import BayesianFunctionsPotprim


the_current_path = os.path.abspath(os.getcwd())
results_path = ".\output"

############## Modes #################################
Sensitivity = True
Bayesian = False

############## Read data #############################
# read the fixed parameter list
inputfileParam = open(
    "datalistInput.json"
)  # input parameters (all) is always in same filenam
AllParam = json.load(inputfileParam)

# Parameteres from AllParam (for Sensitivity)
bact_DOM_rel = AllParam["bact_DOM_rel"]
DOM_EC = AllParam["DOM_EC"]
kpriming = AllParam["kpriming"]
KS = AllParam["KS"]
KSfungi = AllParam["KSfungi"]
KSbact = AllParam["KSbact"]
kMAOMs_MAOMp = AllParam["kMAOMs_MAOMp"]
kPOM_MAOM = AllParam["kPOM_MAOM"]
MAOMpmaxrate = AllParam["MAOMpmaxrate"]
MAOMsmaxrate = AllParam["MAOMsmaxrate"]
MAOMmaxrate = AllParam["MAOMmaxrate"]
MAOMratioSP = AllParam["MAOMratioSP"]
maxEffectBactMAOM = AllParam["maxEffectBactMAOM"]
maxEffectSA_MAOM = AllParam["maxEffectSA_MAOM"]
maxEffectN_MAOM = AllParam["maxEffectN_MAOM"]
MM_N_MAOM = AllParam["MM_N_MAOM"]
MM_Bact_MAOM = AllParam["MM_Bact_MAOM"]
MM_SA_MAOM = AllParam["MM_SA_MAOM"]
MM_DOM_MAOM = AllParam["MM_DOM_MAOM"]
Priming_max = AllParam["Priming_max"]

# Create variables
df_list = []
temp_df_list = []  # temporary df_list to store returned dataframe

############## Sensitivity ###########################
if Sensitivity:
    paramsToTestValues = (
        bact_DOM_rel,
        DOM_EC,
        kpriming,
        KS,
        KSfungi,
        KSbact,
        kPOM_MAOM,
        MAOMpmaxrate,
        MAOMsmaxrate,
        MAOMmaxrate,
        MAOMratioSP,
        maxEffectBactMAOM,
        maxEffectSA_MAOM,
        maxEffectN_MAOM,
        MM_N_MAOM,
        MM_Bact_MAOM,
        MM_SA_MAOM,
        MM_DOM_MAOM,
        Priming_max,
    )

    paramsToTestNames = (
        "bact_DOM_rel",
        "DOM_EC",
        "kpriming",
        "KS",
        "KSfungi",
        "KSbact",
        "kPOM_MAOM",
        "MAOMpmaxrate",
        "MAOMsmaxrate",
        "MAOMmaxrate",
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
    paramChanges = np.array([-50, 0, 100])  # % changes to try for each parameter
    numParams = len(paramsToTestValues)
    numValues = len(paramChanges)
    # numRuns_total= len(paramChanges) * len(paramsToTestValues)  # number of sensitivity runs
    # todo, check how MAOM, MAOMs and MAOMp are calculate throughout the run

    # Treatments
    inputRun = pd.read_csv(
        "Bayesian_run_input.csv", header=0, skiprows=0
    )  # maybe we need to change the file??
    numTreatments = len(inputRun)

    for param in paramsToTestNames:
        for paramChange in paramChanges:
            # calculate by how much to change the parameter value, using a relative parameter change
            delta = (
                paramsToTestDict[param] * paramChange / 100
            )  # I change 1 parameter value
            # caculate new value of parameter
            value = paramsToTestDict[param] + delta
            paramsToTestDict[param] = value

            # for i in range(len(DOMinput_treatments)):
            # numruns = numruns + 1
            # I want to use thevalues from the dict, for sensitivity, so i put all of the values back in the variable (not the fastest way)
            # needs to be changed ifyou change the parameters to test
            bact_DOM_rel = paramsToTestDict["bact_DOM_rel"]
            DOM_EC = paramsToTestDict["DOM_EC"]
            kpriming = paramsToTestDict["kpriming"]
            KS = paramsToTestDict["KS"]
            KSfungi = paramsToTestDict["KSfungi"]
            KSbact = paramsToTestDict["KSbact"]
            kPOM_MAOM = paramsToTestDict["kPOM_MAOM"]
            MAOMpmaxrate = paramsToTestDict["MAOMpmaxrate"]
            MAOMsmaxrate = paramsToTestDict["MAOMsmaxrate"]
            MAOMmaxrate = paramsToTestDict["MAOMmaxrate"]
            MAOMratioSP = paramsToTestDict["MAOMratioSP"]
            maxEffectBactMAOM = paramsToTestDict["maxEffectBactMAOM"]
            maxEffectSA_MAOM = paramsToTestDict["maxEffectSA_MAOM"]
            maxEffectN_MAOM = paramsToTestDict["maxEffectN_MAOM"]
            MM_N_MAOM = paramsToTestDict["MM_N_MAOM"]
            MM_Bact_MAOM = paramsToTestDict["MM_Bact_MAOM"]
            MM_SA_MAOM = paramsToTestDict["MM_SA_MAOM"]
            MM_DOM_MAOM = paramsToTestDict["MM_DOM_MAOM"]

            for treatment in range(numTreatments):
                treatmentVar = inputRun.iloc[treatment, 0:17]

                temp_df_list, Bayesian, Sensitivity = run_model(
                    AllParam, treatmentVar, Bayesian=False, Sensitivity=False
                )
                df_list.append(
                    temp_df_list
                )  # append doesn't work for dataframes, so the lists have to be appended to later use concat

    # in the end, save data output
    try:
        os.makedirs("./output/data")
    except FileExistsError:
        # directory already exists
        pass

    results_df = pd.concat(
        df_list, ignore_index=True
    )  # add all the rows to the results_df
    results_df.to_csv(
        ".\output\data\Sensitivity.csv",
        index=False,
        float_format="%.2f",
    )


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
)
parser.add_argument(
    "-d", "--debug", default=False, action="store_true", help="Print debug info"
)

# Read arguments from command line
args = parser.parse_args()
parallel = args.parallel
debug = args.debug

# number of parameter sets to try, including the start, set very high for calibration (10000)
NumberOfTries = args.tries

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


# "put the measured data in another dataframe per treatment number"
data_measured = pd.DataFrame()
data_measured = inputBayesianRun.iloc[:, 17:81]

# inputBayesianRun({"sample"})

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
data_Simulated = [
    {"resp1": 0, "resp_sub1": 0, "sim likelihood": 0}
    for treatment in range(numTreatments)
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
    treatmentVar = inputBayesianRun.iloc[
        treatment, 0:17
    ]  # to be  corrected for nr of columns needed
    # print(
    #     "treatmentID",
    #     treatmentVar["treatmentID"],
    #     "treatment",
    #     treatmentVar["treatment"],
    # )
    temp_df_list, Bayesian, Sensitivity = run_model(
        AllParam, treatmentVar, Bayesian=False, Sensitivity=False
    )
    df_list.append(
        temp_df_list
    )  # append doesn't work for dataframes, so the lists have to be appended to later use concat

    if Sensitivity:
        temp_df_list.to_csv(
            ".\output\data\Sensitivity.csv", index=False, float_format="%.2f"
        )
        break  # to have only one set of Sensitivity data (for one treatment)

    # we need to couple the output of the right day to the measured output
    data_Simulated[treatment]["resp1"] = temp_df_list["resp"]
    data_Simulated[treatment]["resp_sub1"] = temp_df_list["resp_sub"]
    # we need to add for the treatment the likelyhood of all measurements added, data_measured is df so other indexing

    # print("datasim resp1 treatment 1", data_Simulated[treatment]["resp1"])
    # print("datameasured", data_measured["resp1"][treatment])

    if Bayesian:
        likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
            data_Simulated[treatment]["resp1"],
            data_measured["resp1"][treatment],
            data_measured["resp1_error"][treatment],
        )
        data_Simulated[treatment]["sim likelihood"] += likelyhood


results_df = pd.concat(df_list, ignore_index=True)  # add all the rows to the results_df

# in the end, save data output
try:
    os.makedirs("./output/data")
except FileExistsError:
    # directory already exists
    pass

if (Sensitivity is False) and (Bayesian is False):
    # after running the outermost loop (for three different treatments)
    # merge the three dataframes to create a data output containing all three treatments
    # dfAll=pd.concat(outDataframes)
    # dfAll.to_csv(".\output\data\Output.csv", index=False)
    results_df.to_csv(".\output\data\Output.csv", index=False, float_format="%.2f")


print(results_df)


# 4) calculate the likelihood of each run for each field from the differences between measured and simulated and error
# create empty list(logLi) and store all the differences between measured and simulated per day
likelihood_simulated = 0
for treatment in range(
    numTreatments
):  # For each datapoint (should be the same for both datas) calculate the differance
    likelihood_simulated += data_Simulated[treatment]["sim likelihood"]

# 5) calculate total likelihood of this parameter set over all fields (Log0)
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
for c in range(1, NumberOfTries):  # For each trial Run

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

    # 8) calculate the likelihood of these new parameters, assuming a uniform distribution
    # pdf=probability density function, the likelihood of the parameter set
    loglikelihood_param1 = np.sum(
        np.log(stats.uniform.pdf(candidateValue, MinimalOption, MaximumOption))
    )

    if (
        loglikelihood_param1 > 0
    ):  # if the parameter you want to try is in the range between min and max
        for treatment in range(
            numTreatments
        ):  # so we run for each datapoint measured (=field)

            # 9) run the model for each treatment with the new parameters
            treatmentVar = inputBayesianRun.iloc[treatment]  # to be moved & use iloc
            results_df = run_model(AllParam, treatmentVar, numTreatments, False, False)
            print(results_df)
            # we need to couple the output of the right day to the measured output
            data_Simulated[treatment]["resp1"] = results_df["resp"][0]
            # we need to add for the treatment the likelyhood of all measurements added
            likelyhood = BayesianFunctionsPotprim.calc_sim_likelyhood(
                data_Simulated[treatment]["resp1"],
                data_measured["resp1"][treatment],
                data_measured["resp1_error"][treatment],
            )
            data_Simulated[treatment]["sim likelihood"] += likelyhood

        # 10) calculate the likelihood of each run for each field and store in

        DiffMeasureSimulated = []  # new empty for every try
        for treatment in range(numTreatments):
            DiffMeasureSimulated.append(data_Simulated[treatment]["sim likelihood"])

        """
        11) calculate the likelihood of the entire run over all fields in LogLikelihoodSim1
        the average over all runs within one try of parameters
        """

        log_likelihood_sim1 = sum(DiffMeasureSimulated) / len(data_Simulated)
        # print (DiffMeasureSimulated)

        """
        12) compare the likelihood of this run to the previous and accept into posterior or not
        form the ratio of this step to previous and accept/reject from this (log a/b = log a-log b)
        if the new run fits better it is always accepted (logAlpha>0), 
        if it is worse it is sometimes accepted depending on the random
        """
        alpha = log_likelihood_sim1 / log_likelihood_sim0
        random = ra.random()  # choose random value and take the log
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
            if BayesianFunctions.check_dataframe_significant_change(
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
