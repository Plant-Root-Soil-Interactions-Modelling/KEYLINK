# -*- coding: utf-8 -*-
"""
Created on Mon May 12 14:31:46 2025

@author: Olga
"""
import os
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import csv
from sklearn.metrics import mean_squared_error
import pandas as pd


# setup the logs folder
def create_log_folder(mode_):
    Logs_path = "./logs"
    date_folder = datetime.now().strftime("%y%m%d")

    # Base folder name
    base_folder_name = f"{date_folder}_{mode_}"
    logs_path = os.path.join(Logs_path, base_folder_name)
    counter = 1

    try:
        os.makedirs(logs_path)
        return logs_path, base_folder_name
    except FileExistsError:
        # Try appending _1, _2, etc. directly to base folder name
        while os.path.exists(logs_path):
            base_folder_name2 = f"{base_folder_name}_{counter}"
            logs_path = os.path.join(Logs_path, base_folder_name2)
            counter += 1
        os.makedirs(logs_path)
        return logs_path, base_folder_name2


############## Creating the Respiration Plot #########
def drawRespPlot(
    labels,
    respSoil_mean_model,
    respSoil_mean_measure,
    respSubstrate_mean_model,
    respSubstrate_mean_measure,
    figures_path
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

    # if "Validation" in name:
    #     try:
    #         os.makedirs(os.path.join(logs_path))
    #     except FileExistsError:
    #         # directory already exists
    #         pass

    #     plt.savefig(os.path.join(os.path.join(logs_path, name)))

    # else:
    plt.savefig(os.path.join(figures_path, "RespPlot.png"))

    plt.close()

#leftovers after Terka probably not useful anymore
# ############## RMSE ##################################
# def calculateRMSE(actual, predicted, variable, filepath):
#     actual = np.array(actual)
#     predicted = np.array(predicted)

#     rmse = np.sqrt(((predicted - actual) ** 2).mean())

#     # save as a one csv file
#     try:
#         os.makedirs(filepath)
#     except FileExistsError:
#         # directory already exists
#         pass

#     file_exists = os.path.isfile(os.path.join(filepath, "rmse.csv"))
#     file_is_empty = (
#         file_exists and os.path.getsize(os.path.join(filepath, "rmse.csv")) == 0
#     )

#     with open(os.path.join(filepath, "rmse.csv"), newline="", mode="a") as file:
#         csv_writer = csv.writer(file)

#         # Write the header
#         if not file_exists or file_is_empty:
#             csv_writer.writerow(["Variable", "RMSE"])

#         # Write the key and values to the CSV file
#         csv_writer.writerow([variable, rmse])

#     return rmse


# ############## EF (Nash-Sutcliffe Efficiency) ########
# def calculateEF(actual, predicted, variable, filepath):
#     actual = np.array(actual)
#     predicted = np.array(predicted)

#     ef = 1 - (
#         np.sum((actual - predicted) ** 2) / np.sum((actual - np.mean(actual)) ** 2)
#     )

#     # save as a one csv file
#     try:
#         os.makedirs(filepath)
#     except FileExistsError:
#         # directory already exists
#         pass

#     file_exists = os.path.isfile(os.path.join(filepath, "ef.csv"))
#     file_is_empty = (
#         file_exists and os.path.getsize(os.path.join(filepath, "ef.csv")) == 0
#     )

#     with open(os.path.join(filepath, "ef.csv"), newline="", mode="a") as file:
#         csv_writer = csv.writer(file)

#         # Write the header
#         if not file_exists or file_is_empty:
#             csv_writer.writerow(["Variable", "EF"])

#         # Write the key and values to the CSV file
#         csv_writer.writerow([variable, ef])

#     return ef


# ############## Bias ##################################
# def calculateBias(actual, predicted, variable, filepath):
#     actual = np.array(actual)
#     predicted = np.array(predicted)

#     bias = np.mean(actual - predicted)

#     # save as a one csv file
#     try:
#         os.makedirs(filepath)
#     except FileExistsError:
#         # directory already exists
#         pass

#     file_exists = os.path.isfile(os.path.join(filepath, "bias.csv"))
#     file_is_empty = (
#         file_exists and os.path.getsize(os.path.join(filepath, "bias.csv")) == 0
#     )

#     with open(os.path.join(filepath, "bias.csv"), newline="", mode="a") as file:
#         csv_writer = csv.writer(file)

#         # Write the header
#         if not file_exists or file_is_empty:
#             csv_writer.writerow(["Variable", "Bias"])

#         # Write the key and values to the CSV file
#         csv_writer.writerow([variable, bias])

#     return bias

############## 1:1 plot ##################################
def plot_1to1(ax, measured, modelled, treatment, label, palette):
    # Map treatment values to colors
    colors = treatment.map(palette).tolist() 
    ax.scatter(measured, modelled, color=colors, alpha=0.7)
    
    min_val = min(measured.min(), modelled.min())
    max_val = max(measured.max(), modelled.max())
    
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='1:1 line')
    ax.set_xlabel(f'Measured {label}')
    ax.set_ylabel(f'Modelled {label}')
    ax.set_title(label)
    # ax.legend()
    ax.grid(True)
    
    # Legend (one entry per treatment)
    # handles = []
    # for treatment_value, color in palette.items():
    #     handles.append(plt.Line2D([], [], marker='o', linestyle='', color=color, label=treatment_value))
    # ax.legend(handles=handles, title="Treatment")
    
def calculate_metrics(measured, modelled):
    """Calculate RMSE, bias, and model efficiency (EF) using only valid data pairs"""
    # First, ensure we're working with arrays or Series

    if isinstance(measured, pd.Series):
        measured = measured.values
    if isinstance(modelled, pd.Series):
        modelled = modelled.values
    
    # Find indices where both measured and modelled have valid values
    valid_indices = np.logical_and(~np.isnan(measured), ~np.isnan(modelled))
    
    # Extract only valid pairs
    valid_measured = measured[valid_indices]
    valid_modelled = modelled[valid_indices]
    
    # Check if we have any valid pairs
    if len(valid_measured) == 0:
        return np.nan, np.nan, np.nan
    
    # Root Mean Square Error
    rmse = np.sqrt(mean_squared_error(valid_measured, valid_modelled))
    
    # Bias (mean difference)
    bias = np.mean(valid_modelled - valid_measured)
    
    # Model Efficiency (EF)
    ss_res = np.sum((valid_measured - valid_modelled) ** 2)
    ss_tot = np.sum((valid_measured - np.mean(valid_measured)) ** 2)
    # Avoid division by zero
    ef = 1 - (ss_res / ss_tot) if ss_tot != 0 else np.nan
    
    return rmse, bias, ef

#%% --- combines all calibrations we want to combine into one csv file and calculates min and max of the accepted parameters
def compile_Bayesians():  
    #not sure why this is needed but somehow yes
    os.chdir("C:/Users/Olga/Dropbox/git/KEYLINK")
    #load all the accepted parameter sets and calculate minimum and maximum of each parameter
    #list of all calibrations that I want to string:
    # calibrations = [
    # "250707_Bayesian",
    # "250708_Bayesian",
    # "250710_Bayesian",
    # "250715_Bayesian",
    # "250717_Bayesian",
    # "250717_Bayesian_1",
    # "250719_Bayesian"
    # ]
    
    calibrations = ["250725_Bayesian_2_saved_manually",
                    "250727_Bayesian_2_saved_manually",
                    "250726_Bayesian_saved_manually"
        ]
    

    
    df_list = []
    df_all_list = []
    
    for i in calibrations:        
        file_path1 = os.path.join("./logs/", i, "calibratedParameters.csv")
        # file_path2 = os.path.join("./logs/", i, "AllTestedParameters.csv")
        # Load the CSV files into a DataFrame    
        df = pd.read_csv(file_path1, header=None)
        # df_all = pd.read_csv(file_path2, header=None)
        #append the dataframe to a list of dataframes
        df_list.append(df)
        # df_all_list.append(df_all)
        
    #create dataset from all calibrations
    acceptedParams_df = pd.concat(
        df_list, ignore_index=True
    )     
    # testedParams_df = pd.concat(
        # df_all_list, ignore_index=True
    # )     
    
    #calculate minimum and maximum values for each parameter
    # Min and max per column
    min_values = acceptedParams_df.min()
    max_values = acceptedParams_df.max()
    

    
    # Combine into a new DataFrame
    summary_df = pd.DataFrame({
        'min': min_values,
        'max': max_values
    })

    #save summary to csv
    summary_df.to_csv(
        os.path.join("acceptedParams_MinMax.csv"),
        index=False,
        float_format="%.5f",
    )
    
    #save All accepted params to csv
    acceptedParams_df.to_csv(
        os.path.join("acceptedParams.csv"),
        header=False,
        index=False,
        float_format="%.5f",
    )
    
    return

