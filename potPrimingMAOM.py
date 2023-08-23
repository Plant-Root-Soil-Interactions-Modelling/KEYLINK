# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:08:20 2023

@author: gdeckmyn
"""
import keylink_functions as mf
import matplotlib.pyplot as plt
import numpy as np

time_d=[]
outMAOM=[]
outMAOMsaturation=[]
outPOM=[]
outDOM=[]
outDOM_CN=[]
outBact_RS=[]
outPRIMING=[]
outResp=[]
outRespPriming=[]
DOM_RS=0.05  # rhizosphere DOM [gC/m3]
bact_RS=61 #bacterial biomass [gC/m3], noadd average from Jílková2022, was 2.4
CN_DOM_RSinput=80  #CN of the daily input [unitless], Jílková2022: leachates 80, exudates 6, was 50
MAOMsaturation=0.9 #proportion of MAOM capacity filled [unitless], based on MAOM content from Jílková2022, was 0.5
maxMAOM=  28208 #MAOM capacity [gC/m3] based on silt+clay from Jílková2022 and Georgiou et al. 2022, was 100
MAOMini=MAOMsaturation*maxMAOM #initial MAOM concentration [gC/m3]
MM_DOMtoMAOM=0.025  # DOM concentration for speed being half max speed (Michaelis Menten) [gC/m3]
MAOMmaxrate=0.2 # max proportion of DOM stabilized in MAOM per day [unitless]
surface_RS= 10000   # surface area of all roots/hyphae [m2/m3] ??unit correct
DOMinput=10 #DOM added in each addition [gC/m3], Jílková2022, was 0.5
numDays=150 #number of days of incubation experiment/how long to run the model, Jílková2022
SOMini=38400  # total SOM, only used for priming [gC/m3], Jílková2022, was 150
resp=0 #respiration [gC/m3] ??is unit correct
GMAX=1.24 #growth rate [gC/(gC day)], KEYLINK
DEATH=0.05 #death rate [gC/(gC day)], KEYLINK
rRESP=0.03  #respiration rate resp, [gC/(gC day)], KEYLINK
KS=0.05  # concentration of ?substrate (DOM) for half speed growth, for growing on DOM [gC/m3], related to substrate quality
MCN=0.8 #???, KEYLINK
CN_bact=4 #CN of bacteria, Jílková2022 initial is 10
pH=4.1 #Jílková2022, was 3.5
Nmin=0.00000001 # was 0.0005 [gN/m3] was 0.00000001
SOMCN=24 #CN of SOM, Jílková2022, was 20
PV=15 # volume of micropores [l/m3]
primingIntensity=5 #ratio of POM decayed for DOM decayed [gC/gC, unitless], depends on DOM quality, we know DOM CN which is something else
CN_DOM_RS=CN_DOM_RSinput # set inital DOM CN equal to input
SOM=SOMini
MAOM=MAOMini

for d in range(numDays):
      time_d.append(d)  #store days in an array for plotting
      DOM_N=DOM_RS/CN_DOM_RS
      if (d%14)==0: #where does this if end?
          DOM_RS+=DOMinput/CN_DOM_RS #??why is DOMinput divided by CN
          
      DOM_N+=DOMinput/CN_DOM_RSinput
      CN_DOM_RS=DOM_RS/ DOM_N
      DOM_RS,CN_DOM_RS, bact_RS, SOM, resp, respPriming= mf.calcRhizosphere(MAOMsaturation,maxMAOM, bact_RS, DOM_RS, GMAX, DEATH,CN_bact, CN_DOM_RS, MCN, pH, rRESP, KS, SOMCN, Nmin, SOM, PV, primingIntensity)
                                                         # ((MAOMsaturation,maxMAOM,bact_RS, DOM_RS, gmax, DEATH,CN_bact, CN_DOM_RS, pCN, pH, res, Ks, fCN, CN_SOM, Nmin, SOM,PVstruct,  primingIntensity)        # MAOM formation
      MAOMsaturation=mf.calcMAOMsaturation (maxMAOM,MM_DOMtoMAOM, MAOMsaturation, MAOMmaxrate, DEATH, PV,bact_RS, surface_RS, DOM_RS)
      MAOM=MAOMsaturation*maxMAOM   
      outMAOM.append(MAOM)
      outMAOMsaturation.append(MAOMsaturation)
      outPOM.append(SOM-MAOM)
      outDOM.append(DOM_RS)
      outBact_RS.append(bact_RS)
      outPRIMING.append((SOM-MAOM)-(SOMini-MAOMini))
      outResp.append(resp)
      outRespPriming.append(respPriming)
      
# plt.plot(outMAOM)      
# plt.plot(outDOM)
# plt.plot(outPOM)
# plt.plot(outMAOMsaturation) 
# plt.plot(outBact_RS)   
# plt.plot(outPRIMING) 
# plt.plot(outResp)
# plt.plot(outRespPriming)   

outResp2 =np.divide(outResp, 0.8*24) #change units from gC/m3/day to µg CO2-C/g soil/h
outRespPriming2=np.divide(outRespPriming,0.8*24)
outBact_RS2=np.divide(outBact_RS,0.8) #change units from gC/m3 µgC/g soil

def Dailyplot(outBact_RS, outResp, outRespPriming, outDOM, outPOM, outPRIMING, outMAOM): #plot in original KEYLINK units
    # df2 = pd.DataFrame(df)
    # x = []
    # y = []
    fig, ((p1, p2), (p3,p4),(p5,p6)) = plt.subplots(nrows=3,
                                                                                      ncols=2,
                                                                                      figsize=(10, 12))
    fig.tight_layout(pad=2.0)
    ps = (p1, p2, p3, p4, p5, p6)
    # counter = count(0, 1)
    # columns = list(df)
    ps[0].set_title("Respiration, gC m-3 day-1")
    ps[1].set_title("Microbial biomass, gC m-3")
    ps[2].set_title("DOM, gC m-3")   
    ps[3].set_title("POM, gC m-3")
    ps[4].set_title("POM change, gC m-3")
    ps[5].set_title("MAOM, gC m-3")
    
    p1.plot(time_d, outResp, label="substrate-derived")
    p1.plot(time_d, outRespPriming, label="soil-derived")
    ps[0].legend(loc=(0.1, 0.7), shadow=True) #loc='upper left',

    p2.plot(time_d, outBact_RS, label="bacterial biomass")
    p3.plot(time_d, outDOM)
    p4.plot(time_d, outPOM)
    p5.plot(time_d, outPRIMING)
    p6.plot(time_d, outMAOM)
    # plt.legend(loc=(1.01, 0), shadow=True) #loc='upper right',

def Dailyplot2(outBact_RS2,outResp2, outRespPriming2): #plot in adjusted units matching the data
    # df2 = pd.DataFrame(df)
    # x = []
    # y = []
    fig, (p1, p2) = plt.subplots(nrows=2,
                                                                                      ncols=1,
                                                                                      figsize=(10, 12))
    fig.tight_layout(pad=2.0)
    ps = (p1, p2)
    # counter = count(0, 1)
    # columns = list(df)
    ps[0].set_title("Respiration, µg C-CO2 g-1 soil h-1")
    ps[1].set_title("Microbial biomass, µgC g-1 soil")
   
    p1.plot(time_d, outResp2, label="substrate-derived")
    p1.plot(time_d, outRespPriming2, label="soil-derived")
    ps[0].legend(loc=(1.01, 0), shadow=True) #loc='upper right',

    p2.plot(time_d, outBact_RS2, label="bacterial biomass")
    # plt.legend(loc=(1.01, 0), shadow=True) #loc='upper right',
    
# Dailyplot(outBact_RS,outResp, outRespPriming)
Dailyplot(outBact_RS, outResp, outRespPriming, outDOM, outPOM, outPRIMING, outMAOM)
Dailyplot2(outBact_RS2,outResp2, outRespPriming2)