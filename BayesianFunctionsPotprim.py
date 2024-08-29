# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 11:06:34 2024

@author: gdeckmyn
"""
import json
import os
import random
import sys
import csv
import numpy as np
from numpy import random as ra
from scipy import stats

from macsime.MacsimeCore import run_model


def block_print():
    sys.stdout = open(os.devnull, 'w')


def calc_sim_likelyhood(measurement, simulation, error):
    if measurement == 0:
        sim_likelihood = 0

    else:
        if error == 0:  # assumption for missing error values
            error = measurement / 5
        sim_likelihood = -0.5 * ((measurement - simulation) / error) ** 2 - np.log(error)

    return sim_likelihood





def read_measured_data(num_treatments, inputData, measurementdays):
    
    data_measured = [{} for treatment in range(num_treatments)]
    RespList=inputData['resp'].to_list()
    for treatment in range(num_treatments):
        # read the measured values of yield from json, assume stndard deviation on yield is 30%
        for day in list (measurementdays):
            data_measured[treatment]["Resp"][day] = inputData["Resp"][day]
            data_measured[treatment]["Resp_sub"][day] = inputData["Resp_sub"][day]
        

       
    return data_measured  # TODO dataframe 'outbayesian' for all  simulated data + likelyhood?


def read_treatmenInput(Inputlist_df, treatment):
    run_input = Inputlist_df['Treatments'][treatment]

    return run_input


def read_parameter_data(inputCalibParamfile):
    # open and read file with the range for the parameters to calibrate
    
    parameterForCalibList = json.load(inputCalibParamfile)
    maximum_option_all = parameterForCalibList['MaxParam']  # input is always in same filename
    CalibParamInit = parameterForCalibList['InitParam']  # MaximumOptionDic = json.load(inputfileParamMax)
    minimal_option_all = parameterForCalibList['MinParam']  # input is always in same filename
 
    # create the list where the parameter values will be updated
    calibrated_parameters=CalibParamInit
    
    
    print('initial', calibrated_parameters)  # check
    # calibrated_parameters.append(soil_biota_init)
    keys = [*calibrated_parameters]
    calibrated_parameters_values = [calibrated_parameters[x] for x in keys]

    num_params = len(calibrated_parameters)

    # create 1 array (not list or dict) for minima and maxima of all parameters to calibrate
   
    minimal_option_list = [minimal_option_all[x] for x in keys]
    minimal_option = np.array(minimal_option_list)

    maximum_option_list = [maximum_option_all[x] for x in keys]
    maximum_option = np.array(maximum_option_list)

    return (num_params, CalibParamInit,calibrated_parameters, calibrated_parameters_values, maximum_option, minimal_option,
            keys)


def find_new_parameters(calibrated_parameters_values, variance_parameter_space, minimal_option, maximum_option, keys, Allparam):
    candidate_value = ra.multivariate_normal(calibrated_parameters_values, variance_parameter_space)  # new parameter value
    for j in range(len(calibrated_parameters_values)):  # For each parameter
        ref_min = min(0, candidate_value[j] - minimal_option[j])  # reflection from minimum if you walked too far
        ref_max = max(0, candidate_value[j] - maximum_option[j])  # reflection from maximum
        if ref_min != 0 or ref_max != 0:
            print(f"For parameter {j} we have ref_min {ref_min} and ref_max {ref_max}")
            print(f"Old candidate {candidate_value[j]}")
            candidate_value[j] = candidate_value[j] - 2 * ref_min - 2 * ref_max  # new values
            print(f"New candidate {candidate_value[j]}")
            ref_min = min(0, candidate_value[j] - minimal_option[j])  # reflection from minimum if you walked too far
            ref_max = max(0, candidate_value[j] - maximum_option[j])  # reflection from maximum
            if ref_min != 0 or ref_max != 0:
                candidate_value[j] = random.gauss(calibrated_parameters_values[j], variance_parameter_space[j, j])
        assert (candidate_value[j] > minimal_option[j])
        assert (candidate_value[j] < maximum_option[j])

    # assign the new values to the parameters
    candidate_parameters = dict(zip(keys, candidate_value))  # Combine the new values to the parameters name

    Allparam.update(candidate_parameters)
    return candidate_parameters, candidate_value, Allparam

def save_result(result, path, filename="Output", Bayesian=True):
    
    if Bayesian !=  True:
   #  for treatment in result:
   #      print(f"{row[0]}: {row[1]} {row[2]}")
           Bayesian=Bayesian
    # Save as csv
    with open(os.path.join(path, filename+'.csv'), 'a', newline='') as f:
        writer = csv.writer(f)
        for treatment in result:
            writer.writerows(treatment)


def check_significant_change(values, num_identical_results, alpha=0.05):
    """
    Checks if there is a significant change in the number of chosen results and the previous values before them
    in terms of averages and standard deviations.

    Parameters:
    averages (list of float): List of average values.
    std_devs (list of float): List of standard deviation values.
    num_identical_results(int): number of results to compare
    alpha(float): desired alpha

    Returns:
    bool: True if there is a significant change, False otherwise.
    """

    if len(values) < 2 * num_identical_results:
        return True
        # raise ValueError(f"Lists must contain at least {2 * num_identical_results} elements each.")

    # Perform t-test for averages
    t_stat_avg, p_val_avg = stats.ttest_ind(values[-num_identical_results:],
                                            values[-2 * num_identical_results:-num_identical_results])

    # Perform F-test for variances (standard deviations)
    f_stat_std, p_val_std = stats.f_oneway(values[-num_identical_results:],
                                           values[-2 * num_identical_results:-num_identical_results])

    # Check if there is a significant change
    significant_change = (p_val_avg > alpha) or (p_val_std > alpha)  # may change to 'and' here
    # print('last_results_avg=', last_results_avg,'prev_results_avg', prev_results_avg)
    # print('t_stat_avg=', t_stat_avg,'p_val_avg=', p_val_avg)
    return significant_change


def check_dataframe_significant_change(df, alpha, num_identical_results):
    for column in df.columns:
        if df[column].dtype == 'object' or df[column].dtype.name == 'category':
            raise ValueError(f"Column {column} is not numeric.")
        # values bevat de waarde van 1 parameter over alle runs    
        values = df[column].dropna().values
        # as long as 1 column does not converge (gives significant difference) run continues
        if check_significant_change(values, num_identical_results, alpha):
            return False
    # when all parameters converge return true
    return True
