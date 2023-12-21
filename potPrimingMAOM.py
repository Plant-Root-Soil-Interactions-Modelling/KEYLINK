# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:08:20 2023

@author: gdeckmyn
"""
import keylink_functions as mf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#initialization
#output objects initialized as empty
time_d=[]
outMAOM=[]
outMAOMsaturation=[]
outPOM=[]
outDOMadded=[]
outDOM=[]
outDOM_CN=[]
outBact_RS=[]
outBact_Baseline=[]
outFungi_Baseline=[]
outRespSubstrate=[]
outRespSoil=[]
outRespSoilBaseline=[]

#variables 
availability=np.zeros(3)
bact_total = 5 #total biomass of bacteria [gC/m3]
bact_RS_rel = 0.2 #proportion of rhizosphere bacteria, i.e. those that are able to feed on DOM
bact_RS = bact_total*bact_RS_rel #biomass of bacteria growing on DOM [gC/m3], final noadd average from Jílková2022, was 61, *0.2 because most was li-ving on SOM, only 20% on new DOM
bact_RS_sub = 0 # proportion of bacterial carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
bact = bact_total * (1-bact_RS_rel) #80% of my bacteria are decomposing only POM and MAOM
bact_sub = 0 # proportion of this bacterial carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
BD=800   # bulk density kg m³
claySA=800000 #m²/kg was 8000000 cm²/g
CN_bact=4 #CN of bacteria, Jílková2022 initial is 10
CN_DOMinput=80  #CN of the daily input [unitless], Jílková2022: leachates 80, exudates 6, was 50
CN_DOM=CN_DOMinput # set inital DOM CN equal to input
CN_fungi=8 #KEYLINK
CN_POM=24 #CN of SOM, Jílková2022, was 20
CN_MAOM=15 
DEATH=0.05 #death rate [gC/(gC day)], KEYLINK
DEATHfungi=0.02 #death rate [gC/(gC day)], KEYLINK
DOM=0  # rhizosphere DOM [gC/m3]
DOMinput=10 #DOM added in each addition [gC/m3], Jílková2022, was 0.5
DOM_EC = 5 # DOM energetic quality = energy stored per one gram of DOM [J/g]
DOM_N=DOM/CN_DOM
DOM_sub = 0 #propate [gC/(gC day)], KEYLINK

fungi=1 #based on final noadd in Jílková et al. 2022
fungi_sub = 0 # proportion of this bacterial carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
fClay=0.17 #weight fraction
fSilt=0.24 #weight fraction
GMAX=1.24 #growth rate (DOM) for half speed growth, for growing on DOM [gC/m3], related to substrate quality was 0.001
GMAXfungi=0.6  #growth rate [gC/(gC day)], KEYLINK
k=0.3 #decay rate of negative exponential decay curve of decay price
KS=5  # concentration of ?ortion of DOM carbon that is substrate derived in contrast to soil-derived / values 0 to 1/ is a ratio between substrate-derived C and total C in DOM
KSfungi=20000  # for decaying SOM
KSbact=38000 #for decaying SOM
kPOM_MAOM = 8 #ratio of POM to MAOM decayed / overall SOM decay is partitioned using this fixed ratios really unavailable 
mRecBact=0.5  # how sensitive bact are to recalcitrance
mRecFungi=0.5

MAOMsaturation=0.9 #proportion of MAOM capacity filled [unitless], based on MAOM content from Jílková2022, was 0.5
maxMAOM=  28208 #MAOM capacity [gC/m3] based on silt+clay from Jílková2022 and Georgiou et al. 2022, was 100
MAOMini=MAOMsaturation*maxMAOM #initial MAOM concentration [gC/m3]
MM_DOM_MAOM=0.025  # DOM concentration for speed being half max speed (Michaelis Menten) [gC/m3]
MAOMmaxrate=0.2 # max proportion of DOM stabilized in MAOM per day [unitless]
MAOM=MAOMini
MAOM_sub = 0 # proportion of MAOM that is substrate derived in contrast to soil-derived / values 0 to 1/
maxMAOM = 0.86 * (fClay + fSilt)*100 * BD #[gC/m3] maximum MAOM, 28208 for Jílková et al. 2022 Georgiou et al. 2022: 86 ± 9 and 48 ± 6 mg C/g silt+clay mineral for HM and LM,
MAOMratioSP = 2 #ratio of secondary to primary
MAOMp = MAOM/(MAOMratioSP+1) #primary MAOM [gC/m3]
maxMAOMp = maxMAOM/(MAOMratioSP+1)  # maximum primary MAOM
maxMAOMs = maxMAOM-maxMAOMp  # maximum primary MAOM
MAOMs = MAOM-MAOMp
MAOMpmaxrate = 0.1 #maximum rate of primary MAOM formation
MAOMsmaxrate = 0.1 #maximum rate of secondary MAOM formation

maxEffectBactMAOM=0.9 #(half as slow when no bacteria, rate becomes 1-value)
MM_BactMAOM = 1.10  #bact amount for half speed
maxEffectN_MAOM = 0.9
MM_N_MAOM =  1
maxEffectSA_MAOM = 0.9
MM_SA_MAOM = 0.001  # ratio of SA of hyphae/roots to claysiltSA where half max speed of MAOM formation is reached
Nmin=0.00000001 # was 0.0005 [gN/m3] was 0.00000001
numDays=150 #number of days of incubation experiment/how long to run the model, Jílková2022
pCN=0.8 #???, KEYLINK
PV=15 # volume of micropores [l/m3]
pH=4.1 #Jílková2022, was 3.5
Pmax=10 #maximum decay price [J/gC]
POMini=38400  # total SOM, only used for priming [gC/m3], Jílková2022, was 150
POM=POMini
resp=0 #respiration [gC/m3]
recMAOM= 0.9
rRESP=0.05  #respiration rate resp, [gC/(gC day)], KEYLINK
rRESPfungi=0.03 #respiration rate resp, [gC/(gC day)], KEYLINK
PV=np.array([45,37,37,200,6])  #pore volume distribution
PRadius=np.array([0.05,0.525,8,382.5,875])  #in µm
PSA=np.zeros(5)
PSA=mf.calcPoreSurfaceArea(PV, PRadius, PSA)
MAOMunavail = (PSA[0]/sum(PSA))*MAOMini #the portion of MAOM stored in the smallest pores is really unavailable 
PW=np.array([45/1000,37/1000,37/1000,200/1000,6/1000]) #pore water volume, assume all pores filled , but water is in m³ while volume was in l
siltSA=45.4 #m²/kg
SAclaySilt=claySA*BD*fClay+siltSA*BD*fSilt #total surface area of clay and silt in m²/m³
surface_RS= 10000   # surface area of all roots/hyphae [m2/m3] ??unit correct / look up roots surface area equivalent to that amount of DOM input
fractionSA = surface_RS/SAclaySilt #used in calcMAOMformation? fraction of mineral surface area occupied by roots/hyphae




#run the daily calculations
for d in range(numDays):
    time_d.append(d)  #store days in an array for plotting
    DOM_added = 0
    #on day 0 and then every 14 days, add DOM
    if d==0 or (d%14)==0: #where does this if end?
        DOM_added=DOMinput #to keep track of the additions
        DOM_sub_abs= DOM*DOM_sub #absolute substrate derived C in DOM [gC/m3]
        DOM+=DOMinput #add input to the DOM carbon pool
        DOM_sub_abs+=DOMinput #add all input as substrate derived C
        DOM_sub= DOM_sub_abs/DOM #update relative substrate derived C in DOM
        DOM_N+=DOMinput/CN_DOMinput #add equivalent amount of N to DON pool
        CN_DOM=DOM/DOM_N #calculate new CN of DOM pool
          
    # microbial growth on DOM and priming
    DOM,DOM_N, bact_RS, POM, MAOM, resp, respPriming= mf.calcRhizosphere(POM, CN_POM, MAOM, CN_MAOM, bact_RS, CN_bact, DOM, CN_DOM, GMAX, DEATH, pCN, pH, rRESP, KSbact, DOM_EC, Pmax, k, kPOM_MAOM)
                                                                   # ((MAOMsaturation,maxMAOM,bact_RS, DOM, gmax, DEATH,CN_bact, CN_DOM, pCN, pH, res, Ks, fCN, CN_SOM, Nmin, SOM,PVstruct,  primingIntensity)        # MAOM formation

    #MAOM formation
    DOM =mf.calcMAOMformation (bact_RS, DOM_N, fractionSA, MAOMp, maxMAOMp, DOM, MAOMs, maxMAOMs, MAOMsmaxrate, MAOMpmaxrate, MM_DOM_MAOM,maxEffectBactMAOM,MM_BactMAOM, maxEffectN_MAOM,MM_N_MAOM, maxEffectSA_MAOM,MM_SA_MAOM)
    
    # baseline microbial growth on SOM (without substrate DOM additions)
    MAOMunavail = max((PSA[0]/sum(PSA))*MAOM,MAOMunavail)   #unavail can only go up (unless you disrupt which is not implemented yet), so if more MAOM is stored, equivalent amount will be be always unavailable
    availability=mf.calcAvailPot(PV, PW) #calculates availability of SOM decomposition by bacteria and fungi, separately, from pore size distribution and soil water
    #calculate maximal growth (gmax) for bacteria/fungi on POM/MAOM separately
    gmaxbPOM = mf.calcgmaxmod(CN_bact, CN_POM, pCN, 0.0, 0, pH, 1)*GMAX #gmax for bact on POM
    gmaxfPOM = mf.calcgmaxmod(CN_fungi, CN_POM, pCN, 0.0, 0, pH, 2)*GMAXfungi #gmax for fungi on POM
    gmaxbMAOM = mf.calcgmaxmod(CN_bact, CN_MAOM, pCN, recMAOM, mRecBact, pH, 1)*GMAX #gmax for bact on MAOM
    gmaxfMAOM = mf.calcgmaxmod(CN_fungi, CN_MAOM, pCN, recMAOM, mRecFungi, pH, 2)*GMAXfungi #gmax for fungi on MAOM
    #calculate substrate derived C in bact and fungi
    bact_sub_abs = bact * bact_sub  #absolute substrate derived C in bacteria [gC/m3]
    fungi_sub_abs = fungi * fungi_sub  #absolute substrate derived C in fungi [gC/m3]
    
    #growth equations (dB/dt) for each functional group and for variations in C and N pools
    bactPOMgrowth = mf.calcgrowth(bact, POM, availability[0], gmaxbPOM, KSbact)
    bactMAOMgrowth = mf.calcgrowth(bact, MAOM-MAOMunavail, availability[0], gmaxbMAOM, KSbact)
    dbact = bactPOMgrowth + bactMAOMgrowth - DEATH*bact - rRESP*bact
    
    bact_sub_abs += bactMAOMgrowth*MAOM_sub - DEATH*bact*bact_sub - rRESP*bact*bact_sub #add the corresponding part of growth on MAOM as substrate derived C, subtract correspodning part of death and respiration
    
    
    fungiPOMgrowth = mf.calcgrowth(fungi, POM, availability[1], gmaxfPOM, KSfungi)
    fungiMAOMgrowth = mf.calcgrowth(fungi,MAOM-MAOMunavail, availability[1], gmaxfMAOM, KSfungi)
    dfungi =  fungiPOMgrowth + fungiMAOMgrowth - DEATHfungi*fungi - rRESPfungi*fungi 
    fungi_sub_abs+=fungiMAOMgrowth*MAOM_sub - DEATHfungi*fungi*fungi_sub - rRESPfungi*fungi*fungi_sub #add the corresponding part of growth on MAOM as substrate derived C, subtract death and respiration
    
       
    POM+=-mf.calcgrowth(bact, POM, availability[0], gmaxbPOM, KSbact)-mf.calcgrowth(fungi, POM, availability[1], gmaxfPOM, KSfungi)  
    DOM+=DEATH*bact+DEATHfungi*fungi #add dead bacteria and fungi to DOM
    DOM_sub_abs+=DEATH*bact*bact_sub+DEATHfungi*fungi*fungi_sub #add corresponding part of substrate derived C to DOM 
    DOM_sub= DOM_sub_abs/DOM #relative substrate derived C in DOM
    #print(DOM, bact)
    DOM_N+=DEATH*bact/CN_bact+DEATHfungi*fungi/CN_fungi
    MAOM+=-mf.calcgrowth(bact, MAOM-MAOMunavail, availability[0], gmaxbMAOM, KSbact)-   \
        mf.calcgrowth(fungi,MAOM-MAOMunavail, availability[1], gmaxfMAOM, KSfungi)
    baselineRespBact=rRESP*bact
    bact+=dbact
    bact_sub= bact_sub_abs/bact #update relative substrate derived C in bacteria
    baselineRespFungi=rRESP*fungi
    fungi+=dfungi
    fungi_sub= fungi_sub_abs/fungi #update relative substrate derived C in fungi       

 
    #add up soil-derived respiration

    baselineResp = baselineRespBact + baselineRespFungi
    respSoil = baselineResp + respPriming
    outDOMadded.append(DOM_added)
    outMAOM.append(MAOM)
    outMAOMsaturation.append(MAOMsaturation)
    outPOM.append(POM)
    outDOM.append(DOM)
    outBact_RS.append(bact_RS)
    outBact_Baseline.append(bact)
    outFungi_Baseline.append(fungi)
    outRespSubstrate.append(resp)
    outRespSoilBaseline.append(baselineResp)
    outRespSoil.append(respSoil)
       
      
# plt.plot(outMAOM)      
# plt.plot(outDOM)
# plt.plot(outPOM)
# plt.plot(outMAOMsaturation) 
# plt.plot(outBact_RS)   
# plt.plot(outResp)
# plt.plot(outRespPriming)   
# plt.plot(outBact_Baseline) 
# plt.plot(outBact_RS) 
# plt.plot(outFungi_Baseline) 

#change units to easily understandable for the plot
outDOMadded2 = np.divide(outDOMadded,0.8) #change units from gC/m3 µgC/g soil
outBact_RS2=np.divide(outBact_RS,0.8) #change units from gC/m3 µgC/g soil
outBact_Baseline=np.divide(outBact_Baseline,0.8)
outFungi_Baseline=np.divide(outBact_Baseline,0.8)
outRespSubstrate2 =np.divide(outRespSubstrate, 0.8*24) #change units from gC/m3/day to µg CO2-C/g soil/h
outRespSoilBaseline2=np.divide(outRespSoilBaseline, 0.8*24)
outRespSoil2=np.divide(outRespSoil,0.8*24)
outDOM2 = np.divide(outDOM,0.8) #change units from gC/m3 µgC/g soil
outPOM2 = np.divide(outPOM,0.8*1000) #change units from gC/m3 mgC/g soil
outMAOM2 = np.divide(outMAOM,0.8*1000) #change units from gC/m3 mgC/g soil

#combine output arrays into a dataframe and save it to csv
df = pd.DataFrame({"DOMaddition" : outDOMadded2,
                   "DOM" : outDOM2,
                   "bact_RS" : outBact_RS2,
                   "bact_Baseline" : outBact_Baseline,
                   "resp_substrate" : outRespSubstrate2,
                   "resp_soil_baseline" : outRespSoilBaseline2,
                   "resp_soil" : outRespSoil2,
                   "POM" : outPOM2,
                   "MAOM" : outMAOM2
                   })
df.to_csv(".\output\data\Output.csv", index=False)

# def Dailyplot(outDOMadded, outDOM, outBact_RS, outRespSubstrate, outRespSoil, outRespSoilBaseline, outPOM, outMAOM): #plot in original KEYLINK units
#     # df2 = pd.DataFrame(df)
#     # x = []
#     # y = []
#     fig, ((p1, p2), (p3,p4),(p5,p6)) = plt.subplots(nrows=3,
#                                                                                       ncols=2,
#                                                                                       figsize=(10, 12))
#     fig.tight_layout(pad=2.0)
#     ps = (p1, p2, p3, p4, p5, p6)
#     # counter = count(0, 1)
#     # columns = list(df)
#     ps[0].set_title("DOM additions, gC m-3 day -1")
#     ps[1].set_title("DOM, gC m-3") 
#     ps[2].set_title("Microbial biomass, gC m-3")
#     ps[3].set_title("Respiration, gC m-3 day-1")
#     ps[4].set_title("POM, gC m-3")
#     ps[5].set_title("MAOM, gC m-3")





#plot in adjusted units matching the data
def Dailyplot2(outDOMadded2, outDOM2, outBact_RS2, outRespSubstrate2, outRespSoil2, outRespSoilBaseline2, outPOM2, outMAOM2): #plot in original KEYLINK units
    # df2 = pd.DataFrame(df)
    # x = []
    # y = []
    # fig, ((p1, p2), (p3,p4),(p5,p6)) = plt.subplots(nrows=3,
    #                                                                                   ncols=2,
    #    figsize=(5, 6))
    fig, ((p1, p2, p3), (p4, p5, p6)) = plt.subplots(nrows=2, ncols=3,figsize=(12, 7))#was 10,12
    fig.tight_layout(pad=2.0)
    ps = (p1, p2, p3, p4, p5, p6)
    # counter = count(0, 1)
    # columns = list(df)
    ps[0].set_title("DOM additions, µgC g-1 soil")
    ps[1].set_title("DOM, µgC g-1 soil") 
    ps[2].set_title("Microbial biomass, µgC g-1 soil")
    ps[3].set_title("Respiration, µg C-CO2 g-1 soil h-1")
    # ps[4].set_title("POM, mgC g-1 soil")
    ps[4].set_title("SOM, mgC g-1 soil")
    
    p1.plot(time_d, outDOMadded2)
    p2.plot(time_d, outDOM2)
    p3.plot(time_d, outBact_RS2, label="bacterial biomass")
    p4.plot(time_d, outRespSubstrate2, label="substrate-derived")
    p4.plot(time_d, outRespSoil2, label="soil-derived incl. priming")
    p4.plot(time_d, outRespSoilBaseline2, label="soil-derived baseline")
    ps[3].legend(loc=(0.03, 0.7), shadow=True) #loc='upper left',
    p5.plot(time_d, outPOM2, label="POM")
    p5.plot(time_d, outMAOM2, label="MAOM")
    ps[4].legend(loc=(0.03, 0.03), shadow=True) #loc='upper left',
    # p6.plot(time_d, outMAOM2)
    

# Dailyplot(outDOMadded, outDOM, outBact_RS, outRespSubstrate, outRespSoil, outRespSoilBaseline, outPOM, outMAOM)
Dailyplot2(outDOMadded2, outDOM2, outBact_RS2, outRespSubstrate2, outRespSoil2, outRespSoilBaseline2, outPOM2, outMAOM2)
plt.savefig("output/figures/Dailyplot2.png")