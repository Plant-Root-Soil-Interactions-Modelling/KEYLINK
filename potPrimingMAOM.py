# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:08:20 2023

@author: gdeckmyn
"""
import keylink_functions as mf
import matplotlib.pyplot as plt

outMAOM=[]
outMAOMsaturation=[]
outDOM=[]
outDOM_CN=[]
outBact_RS=[]
outPRIMING=[]
DOM_RS=0.05  # rhizosphere DOM
bact_RS=0.01
CN_DOM_RSinput=50  #CN of the daily input

MAOMsaturation=0.5
maxMAOM=  100#TODO needs to be set to literature values clay & silt
MAOMini=MAOMsaturation*maxMAOM
MM_DOMtoMAOM=0.025  # DOM concentration for speed being half max speed (Michaelis Menten)
MAOMmaxrate=0.2 # max proportion of DOM stabilized in MAOM per day
surface_RS= 10000   # surface area of all roots/hyphae (in m2)
DOMinput=0.5
numDays=1500
SOMini=150  # total SOM, only used for priming
resp=0
GMAX=0.3
DEATH=0.1
rRESP=0.1  
KS=0.05  # concentration for half speed growth, for growin on DOM
MCN=0.8
CN_bact=4 
pH=3.5
Nmin=0.0005
SOMCN=20
PV=15 # volume of micropores 
primingIntensity=0.01
CN_DOM_RS=CN_DOM_RSinput # set inital DOM CN equal to input
SOM=SOMini
MAOM=MAOMini

for d in range(numDays):
      DOM_N=DOM_RS/CN_DOM_RS
      DOM_RS+=DOMinput/CN_DOM_RS
      DOM_N+=DOMinput/CN_DOM_RSinput
      CN_DOM_RS=DOM_RS/ DOM_N
      DOM_RS,CN_DOM_RS, bact_RS, SOM, resp= mf.calcRhizosphere(MAOMsaturation,maxMAOM, bact_RS, DOM_RS, GMAX, DEATH,CN_bact, CN_DOM_RS, MCN, pH, rRESP, KS, SOMCN, Nmin, SOM, PV, primingIntensity)
                                                         # ((MAOMsaturation,maxMAOM,bact_RS, DOM_RS, gmax, DEATH,CN_bact, CN_DOM_RS, pCN, pH, res, Ks, fCN, CN_SOM, Nmin, SOM,PVstruct,  primingIntensity)        # MAOM formation
      MAOMsaturation=mf.calcMAOMsaturation (maxMAOM,MM_DOMtoMAOM, MAOMsaturation, MAOMmaxrate, DEATH, PV,bact_RS, surface_RS, DOM_RS)
      outMAOM.append(MAOMsaturation*maxMAOM)
      outMAOMsaturation.append(MAOMsaturation)
      outDOM.append(DOM_RS)
      outBact_RS.append(bact_RS)
      outPRIMING.append((SOM-MAOMsaturation*maxMAOM)-(SOMini-MAOMini))
plt.plot(outMAOM)      
plt.plot (outDOM)
plt.plot (outMAOMsaturation) 
plt.plot(outBact_RS)   
plt.plot(outPRIMING)  