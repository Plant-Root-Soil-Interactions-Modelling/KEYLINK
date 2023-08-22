# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:08:20 2023

@author: gdeckmyn
"""
import keylink_functions as mf
import matplotlib.pyplot as plt

time_d=[]
outMAOM=[]
outMAOMsaturation=[]
outDOM=[]
outDOM_CN=[]
outBact_RS=[]
outPRIMING=[]
outResp=[]
outRespPriming=[]
DOM_RS=0.05  # rhizosphere DOM
bact_RS=2.4
CN_DOM_RSinput=50  #CN of the daily input

MAOMsaturation=0.5
maxMAOM=  100#TODO needs to be set to literature values clay & silt
MAOMini=MAOMsaturation*maxMAOM
MM_DOMtoMAOM=0.025  # DOM concentration for speed being half max speed (Michaelis Menten)
MAOMmaxrate=0.2 # max proportion of DOM stabilized in MAOM per day
surface_RS= 10000   # surface area of all roots/hyphae (in m2)
DOMinput=0.5
numDays=150
SOMini=150  # total SOM, only used for priming
resp=0
GMAX=1.24

DEATH=0.05
rRESP=0.03  
KS=0.05  # concentration for half speed growth, for growin on DOM
MCN=0.8
CN_bact=4 
pH=3.5
Nmin=0.00000001 # was 0.0005
SOMCN=20
PV=15 # volume of micropores 
primingIntensity=5
CN_DOM_RS=CN_DOM_RSinput # set inital DOM CN equal to input
SOM=SOMini
MAOM=MAOMini

for d in range(numDays):
      time_d.append(d)  #store days in an array for plotting
      DOM_N=DOM_RS/CN_DOM_RS
      if (d%14)==0:
          DOM_RS+=DOMinput/CN_DOM_RS
      DOM_N+=DOMinput/CN_DOM_RSinput
      CN_DOM_RS=DOM_RS/ DOM_N
      DOM_RS,CN_DOM_RS, bact_RS, SOM, resp, respPriming= mf.calcRhizosphere(MAOMsaturation,maxMAOM, bact_RS, DOM_RS, GMAX, DEATH,CN_bact, CN_DOM_RS, MCN, pH, rRESP, KS, SOMCN, Nmin, SOM, PV, primingIntensity)
                                                         # ((MAOMsaturation,maxMAOM,bact_RS, DOM_RS, gmax, DEATH,CN_bact, CN_DOM_RS, pCN, pH, res, Ks, fCN, CN_SOM, Nmin, SOM,PVstruct,  primingIntensity)        # MAOM formation
      MAOMsaturation=mf.calcMAOMsaturation (maxMAOM,MM_DOMtoMAOM, MAOMsaturation, MAOMmaxrate, DEATH, PV,bact_RS, surface_RS, DOM_RS)
      outMAOM.append(MAOMsaturation*maxMAOM)
      outMAOMsaturation.append(MAOMsaturation)
      outDOM.append(DOM_RS)
      outBact_RS.append(bact_RS)
      outPRIMING.append((SOM-MAOMsaturation*maxMAOM)-(SOMini-MAOMini))
      outResp.append(resp)
      outRespPriming.append(respPriming)
# plt.plot(outMAOM)      
# plt.plot (outDOM)
# plt.plot (outMAOMsaturation) 
# plt.plot(outBact_RS)   
# plt.plot(outPRIMING) 
# plt.plot(outResp)
# plt.plot(outRespPriming)   


def Dailyplot(outBact_RS,outResp, outRespPriming):
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
    ps[0].set_title("Respiration, gC-CO2 g-1 h-1")
    ps[1].set_title("Microbial biomass, gC m-3")
   
    p1.plot(time_d, outResp, label="substrate-derived")
    p1.plot(time_d, outRespPriming, label="soil-derived")
    ps[0].legend(loc=(1.01, 0), shadow=True) #loc='upper right',

    p2.plot(time_d, outBact_RS, label="bacterial biomass")
    # plt.legend(loc=(1.01, 0), shadow=True) #loc='upper right',

Dailyplot(outBact_RS,outResp, outRespPriming)