# -*- coding: utf-8 -*-
"""
Created on Sat Nov 29 12:26:27 2025

@author: Olga
"""
#running with second option parameter set, keeping all
import runpy

def run_with_path(path_value):
    runpy.run_path("Potprim_main.py", init_globals={"path_bayesian": path_value})

paths = [
    # "./input files crossvalidation 2022/calibration/Bayesian_run_input_2022_subset1.csv",
    "./input files crossvalidation 2022/calibration/Bayesian_run_input_2022_subset2.csv",
    "./input files crossvalidation 2022/calibration/Bayesian_run_input_2022_subset3.csv",
    "./input files crossvalidation 2022/calibration/Bayesian_run_input_2022_subset4.csv",
    "./input files crossvalidation 2022/calibration/Bayesian_run_input_2022_subset5.csv"
    
]

for p in paths:
    print("************************************************************************")
    print("running for:")
    print(p)
    print("************************************************************************")
    run_with_path(p)
