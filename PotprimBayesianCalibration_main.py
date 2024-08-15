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

from macsime import BayesianFunctions
from macsime.plots import bayesian_plots
from macsime.utilities import get_results_path, save_result

the_current_path = os.path.abspath(os.getcwd())


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


def run_model(inputData, results_path, num_treatments, data=None, parallel=False):
    print("Start running model")
    start = time.perf_counter()
    
    for treatment in range(num_treatments):
        print('treatment =', treatment)
           treatment, result = potPrimingMAOM.run_model_bayesian(inputData, results_path, treatment, data)
            data_Simulated[treatment] = result

    end = time.perf_counter()
    print(f'model ran for {time.strftime("%H:%M:%S", time.gmtime(end - start))}')


if __name__ != '__main__':
    exit(0)

# Initialize parser
parser = argparse.ArgumentParser(description="Run Bayesian optimization")
parser.add_argument("-p", "--parallel", default=False, action='store_true', help="Run in parallel mode")
parser.add_argument("-f", "--fields", default=0, type=int,
                    help="Run maximum this number of fields (default: all fields)")
parser.add_argument("-t", "--tries", default=10000, type=int,
                    help="Run this number of tries (default: 10000)")
parser.add_argument("-d", "--debug", default=False, action='store_true', help="Print debug info")

# Read arguments from command line
args = parser.parse_args()
parallel = args.parallel
debug = args.debug

# number of parameter sets to try, including the start, set very high for calibration (10000)
NumberOfTries = args.tries

t1 = time.perf_counter()

results_path = get_results_path()
# print("Directory ", results_path, " created")

# file_name = 'aa'  # ? who wrote this?
# management_df=CornManagement_df
# plant_input_df=corn_input_df

# input should be in filename below

    # input is always in same filename, contains all field data from database
inputfile = open(os.path.join('', 'BayesianPotprim.json'))
    #   inputfile = open('corn_datalistBayesian_2fields_1double.json')
    #   inputfile = open('corn_datalistBayesianDebug1field1year.json')
    #   inputfile = open('corn_datalistBayesianDebug1field1day.json')
    #   inputfile = open('corn_datalistBayesian_noNan1field.json')
inputData = json.load(inputfile)
    
# the number of fields equals the number of different management files.
# Assume each management always has its own soil input

numFields = len(totallist_df['Fields'])
if args.fields > 0:
    numFields = min(args.fields, numFields)

# read the Parameter data to calibrate
numParams, treatments, Parameters, ParameterValues, MaximumOption, MinimalOption, keys = (
    BayesianFunctions.read_parameter_data(inputData))

# read the measured data (towards which to calibrate)
data_measured = BayesianFunctions.read_measured_data(numTreatments, inputData)

# create list of all parameter sets tried out
priorChain = np.zeros([NumberOfTries, numParams])  # list of tries per parameter
priorChain[0, :] = list(CalibratedParameters.values())  # prior is every parameter set you try out, put in first row
logLseries = []  # create empty list

# create list of all parameter sets accepted, and 0 when not accepted so linnenrs are equal to priorChain
posteriorChain = np.zeros([NumberOfTries, numParams])
# start the chain, will hold all accepted parameter sets, and 0 if not accepted
posteriorChain[0, :] = list(CalibratedParameters.values())
# print(parameterlist_df)

# create list for simulations
data_Simulated = [{'Yield': 0, 'SOC-PostHarvest': 0, 'SOC-Spring': 0} for treatment in range(numTreatments)]

'''
 actual start of calibration
'''

# 1) calculate the variance of the parameter space
# variance of the parameter space, is needed to define the step size for each parameter
VarianceParameterSpace = np.diag(((0.005 * (MaximumOption - MinimalOption)) ** 2))

# 2) calculate the likelihood of the parameters chosen (for a flat distribution this will always be constant or 0)
# pdf=probability density function, the likelihood of the parameter set
loglikelihood_param = np.sum(np.log(stats.uniform.pdf(CalibratedParametersValues, MinimalOption,
                                                      MaximumOption)))
# 3) Simulated Data in a similar frame as the measured values,
run_model(totallist_df, results_path, numFields, data_measured, parallel)
print(data_Simulated)

Csequestered = BayesianFunctions.calc_sequestration(totallist_df, numFields, data_Simulated)
print('C seq=', Csequestered)

# 4) calculate the likelihood of each run for each field from the differences between measured and simulated and error
# create empty list(logLi) and store all the differences between measured and simulated per day
likelihood_simulated = 0
for Fieldnr in range(numFields):  # For each datapoint (should be the same for both datas) calculate the differance
    likelihood_simulated += data_Simulated[Fieldnr]['sim likelihood']

# 5) calculate total likelihood of this parameter set over all fields (Log0)
log_likelihood_sim0 = likelihood_simulated / len(data_measured)
logLseries.append(log_likelihood_sim0)

#  6) save best fit, "BestFitParam" is the parameter set giving the highest likelihood (best fit = maximum probability)
BestFitParam = CalibratedParametersValues  # the initial values of parameters are my best try at first step
# psetMAP is max fit point, save parameter and likelihood of best run
log_likelihood_best_fit_param = loglikelihood_param + log_likelihood_sim0

print("start saving results")
save_result([[data_Simulated]], results_path, "SimdataAll")
save_result([[Csequestered]], results_path, "Csequestered")
save_result([[data_Simulated]], results_path, "SimdataBestFit")
save_result([[CalibratedParametersValues]], results_path, "calibratedParameters")
save_result([[[log_likelihood_sim0]]], results_path, "logLikelyhood")

print('resultspath', results_path)

'''
loop over number ot tries
'''
for c in range(1, NumberOfTries):  # For each trial Run

    # 7) find new parameter values to try
   
    candidateparameters, candidateValue, crops, soilbiota = BayesianFunctions.find_new_parameters(
        CalibratedParametersValues, VarianceParameterSpace, MinimalOption, MaximumOption, keys, soilbiota, crops)
    
    # 8) calculate the likelihood of these new parameters, assuming a uniform distribution
    # pdf=probability density function, the likelihood of the parameter set
    loglikelihood_param1 = np.sum(np.log(stats.uniform.pdf(candidateValue, MinimalOption,
                                                           MaximumOption)))

    if loglikelihood_param1 > 0:  # if the parameter you want to try is in the range between min and max
        # for Fieldnr in range(numFields):  # so we run for each datapoint measured (=field)

        # 9) run the model for each field with the new parameters
        run_model(totallist_df, results_path, numFields, data_measured, parallel)

        # 10) calculate the likelihood of each run for each field and store in

        DiffMeasureSimulated = []  # new empty for every try
        for Fieldnr in range(numFields):
            DiffMeasureSimulated.append(data_Simulated[Fieldnr]['sim likelihood'])

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
        logalpha = (loglikelihood_param1 + log_likelihood_sim1) - (loglikelihood_param + log_likelihood_sim0)
        lograndom = math.log(ra.random())  # choose random value and take the log
        logLseries.append(log_likelihood_sim1)

        # print('random value', lograndom)
        # print('logalpha', logalpha)
        save_result([[data_Simulated]], results_path, "SimdataAll")

        if lograndom < logalpha:
            CalibratedParametersValues = candidateValue
            loglikelihood_param = loglikelihood_param1
            log_likelihood_sim0 = log_likelihood_sim1  # if accepted move to this point
            posteriorChain[c, :] = CalibratedParametersValues  # add step to the chain

            """
            13) set parameters as best fit if they are better than current best fit
            test if we have a new best fit
            """

            if (loglikelihood_param + log_likelihood_sim0) > log_likelihood_best_fit_param:
                log_likelihood_best_fit_param = (loglikelihood_param + log_likelihood_sim0)
                BestFitParam = CalibratedParametersValues  # update most likely parameter set

            Csequestered = BayesianFunctions.calc_sequestration(totallist_df, numFields, data_Simulated)
            
            save_result([[Csequestered]], results_path, "Csequestered")
            save_result([[data_Simulated]], results_path, "SimdataBestFit")
            save_result([[CalibratedParametersValues]], results_path, "calibratedParameters")
            save_result([[[log_likelihood_sim0]]], results_path, "logLikelihood   ")

            """
            14) test if we have enough runs: avg and stdev are table for each column of posterior
            """

            parameters = pd.read_csv(os.path.join(results_path, 'calibratedParameters.csv'))               
            if BayesianFunctions.check_dataframe_significant_change(parameters, alpha=0.5, num_identical_results=50):
                print('Hurraaayyy!!! converged')
                break

# the prior chain saves all tries, also the ones that are not 'saved' in the posterior chain
    priorChain[c, :] = candidateValue

'''
end of loop
'''

t2 = time.perf_counter()

print(f'Calibration ran for {time.strftime("%H:%M:%S", time.gmtime(t2 - t1))}\n')
# df = pd.DataFrame(priorChain, columns=['a0Photo_Eff', 'SoilRootCond', 'minimalStomatalResistance', 'a9DistriRoot',
#                                        'a10DistriFruit', 'a91DistriHyphae', 'ratioExudRoot', 'f_MIC_DOM',
#                                        'f_FOM_fPOM', 'f_fPOM_oPOM', 'MAOMsmaxrate', 'MAOMpmaxrate', 'AgRate',
#                                        'MIC_gmax', 'Mic_RespRate', 'MIC_TR', 'HyphalExploration',
#                                        'HyphaeTurnoverrate', 'Hyphae_fGRSP', 'SoilHyphaeCond'])
df = pd.DataFrame(priorChain, columns=list(CalibratedParameters.keys()))

bayesian_plots(df=df, path=file_name, columns=5, save_to_file=True)
