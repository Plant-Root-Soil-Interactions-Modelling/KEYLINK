# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 11:12:59 2024

@author: Olga
"""
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DOM_EC = 2 # DOM energetic quality = energy stored per one gram of DOM [J/g] was 5
kpriming=0.001 #decay rate of negative exponential decay curve of decay price, was 0.3
Priming_max=0.5 #maximum decay price [J/gC] should be maximum fraction that can be primed
ExtraGrowth = 20
#what can be primed
#SOMprimable=POM+MAOMs+MAOMp  # MAOMp is primed, decision 13/8/2024
SOMprimable = 10000

#changing DOM_EC
DOM_EC_range = np.array([0.5,1,2,5,10])

for i in range(len(DOM_EC_range)):
    #DecayCost = how much energy will be spent on SOM decay, definite integral of a decay price function [J]
    DOM_E = ExtraGrowth*DOM_EC_range[i] # total energy stored in the DOM that bacteria can still assimilate [J]
         
    #SOMprimed = [gC] how much gC in POM or MAOM can be decayed with energy in DOM (DOM_E)
    SOMprimed=Priming_max*SOMprimable*(1-math.exp(-kpriming*DOM_E)) 
    print('DOM EC=', DOM_EC_range[i], 'DOM_E=', DOM_E, 'SOMprimed=', SOMprimed)

#with higher DOM_EC, more SOM can be primed



column_names=['DOM_EC','DOM_E', 	'kpriming', 'Priming_max',	'SOMprimed'	, 'SOMprimable']
results_df = pd.DataFrame(columns=column_names)

#changing Priming max
Priming_max = np.array([0.5,1,2,5,10]) #maximum decay price [J/gC]
for i in range(len(Priming_max)):
    #DecayCost = how much energy will be spent on SOM decay, definite integral of a decay price function [J]
    DOM_E = ExtraGrowth*DOM_EC # total energy stored in the DOM that bacteria can still assimilate [J]
         
    #SOMprimed = [gC] how much gC in POM or MAOM can be decayed with energy in DOM (DOM_E)
    SOMprimed=Priming_max[i]*SOMprimable*(1-math.exp(-kpriming*DOM_E)) 
    
    print('Priming_max=', Priming_max[i], 'DOM_E=', DOM_E, 'SOMprimed=', SOMprimed)
    
# with higher Priming_max, Priming is higher / weird!!


#changing third part of equation
#changing DOM_EC
DOM_EC_range = np.array([0.5,1,2,5,10])
Priming_max=0.5

for i in range(len(DOM_EC_range)):
    #DecayCost = how much energy will be spent on SOM decay, definite integral of a decay price function [J]
    DOM_E = ExtraGrowth*DOM_EC_range[i] # total energy stored in the DOM that bacteria can still assimilate [J]
                 
    #SOMprimed = [gC] how much gC in POM or MAOM can be decayed with energy in DOM (DOM_E)
    bracket = (1-math.exp(-kpriming*DOM_E)) 
    SOMprimed=Priming_max*SOMprimable*(1-math.exp(-kpriming*DOM_E)) 
    # SOMprimed=SOMprimable* (1-math.exp(-kpriming*DOM_E))/Priming_max 
    results_df.loc[len(results_df)] = [DOM_EC, DOM_E, kpriming, Priming_max, SOMprimed, SOMprimable]

    print('bracket', bracket, 'Priming_max=', Priming_max, 'DOM_E=', DOM_E, 'SOMprimed=', SOMprimed)
    
    # scatter plot
results_df.plot(kind='scatter',
        x='DOM_E',
        y='SOMprimed',
        color='red')

# set the title
plt.title('ScatterPlot')

# show the plot
plt.show()


#changing kriming max
kpriming = np.array([0.001,0.01,0.1,1]) #decay rate of negative exponential decay curve of decay price
Priming_max=0.5 #maximum decay price [J/gC] should be maximum fraction that can be primed
DOM_E = np.array([1, 10, 100])

for i in range(len(kpriming)):
    for j in range(len(DOM_E)):
    #DecayCost = how much energy will be spent on SOM decay, definite integral of a decay price function [J]
    # DOM_E = ExtraGrowth*DOM_EC # total energy stored in the DOM that bacteria can still assimilate [J]
         
        #SOMprimed = [gC] how much gC in POM or MAOM can be decayed with energy in DOM (DOM_E)
        SOMprimed=Priming_max*SOMprimable*(1-math.exp(-kpriming[i]*DOM_E[j])) 
        
        print('kpriming=', kpriming[i], 'DOM_E=', DOM_E[j], 'SOMprimed=', SOMprimed)
    
# with higher kpriming, more SOM primed in total but also less response of priming to energy available