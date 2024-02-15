# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:08:20 2023

@author: gdeckmyn
"""
import keylink_functions as mf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#output dataframe list
outDataframes=[]

#********************************************************************************
#calibrated parameters (= that do NOT change during run but need calibration)

# the ones we want to calibrate
bact_DOM_rel = 0.2 #proportion of bacteria that have access to feeding on DOM (e.g. that are present in rhizophere)
DOM_EC = 5 # DOM energetic quality = energy stored per one gram of DOM [J/g]
kpriming=0.3 #decay rate of negative exponential decay curve of decay price
KS=5  # C content required to get half the maximal growth of bacteria when decaying DOM [gC/m3]
KSfungi=20000  # C content required to get half the maximal growth of fungi when decaying SOM [gC/m3]
KSbact=38000 # C content required to get half the maximal growth of bacteria when decaying SOM [gC/m3]
kPOM_MAOM = 8 #ratio of POM to MAOM decayed / overall SOM decay is partitioned using this fixed ratios really unavailable, is k-POM/k_MAOM in israel code 
MAOMpmaxrate = 0.1 #maximum rate of primary MAOM formation
MAOMsmaxrate = 0.1 #maximum rate of secondary MAOM formation
MAOMmaxrate=0.2 # max proportion of DOM stabilized in MAOM per day [unitless]
MAOMratioSP = 2 #ratio of secondary to primary MAOM
maxEffectBactMAOM=0.9 #(half as slow when no bacteria, rate becomes 1-value)
maxEffectSA_MAOM = 0.9
maxEffectN_MAOM = 0.9
MM_N_MAOM =  1
MM_Bact_MAOM = 1.10  #bact amount for half speed
MM_SA_MAOM = 0.001  # ratio of SA of hyphae/roots to claysiltSA where half max speed of MAOM formation is reached
MM_DOM_MAOM=0.025  # DOM concentration for speed being half max speed (Michaelis Menten) [gC/m3]
Priming_max=10 #maximum decay price [J/gC]


# the ones we use from other calibration

GMAX=1.24 #maximal growth rate for bacteria [gC/(gC day)], KEYLINK
GMAXfungi=0.6  #maximal growth rate for fungi [gC/(gC day)], KEYLINK
mRecBact=0.5  # how sensitive bact are to recalcitrance
mRecFungi=0.5 #
# resp=0.01 #respiration rate for bacteria growing on DOM / ??do we really need a different one? it was set to 0 decided to ditch it and just the next one
rRESP=0.05  #respiration rate resp, [gC/(gC day)], KEYLINK
rRESPfungi=0.03 #respiration rate resp, [gC/(gC day)], KEYLINK
DEATH=0.05 #death rate for bacteria [gC/(gC day)], KEYLINK
DEATHfungi=0.02 #death rate for fungi [gC/(gC day)], KEYLINK
pCN=0.8 #sensitivity to CN ratio of consumed substrate, values 0-1, taken from KEYLINK (value for bacteria)
recMAOM= 0.9 #recalcitrance of MAOM, (recalcitrance of POM assumed 0)

# input parameters that do not change but vary for different treatments of experiment (and runs)
CN_DOMinput=6  #CN of the daily input [unitless], Jílková2022: leachates 80, exudates 6

# input parameters that do not change (=measurable) and are not calibrated, just 'start situation"
#Nmin=0.00000001 # was 0.0005 [gN/m3] was 0.00000001 not used in the model currently
BD=800   # bulk density [kg/m³]
claySA=800000 # surface area of clay [m²/kg] was 8000000 cm²/g
CN_bact=4 #CN of bacteria, from KEYLINK, in Jílková2022 initial CN of microbial biomass is 10
CN_fungi=8 #KEYLINK
CN_POM=24 #CN of SOM, Jílková2022
CN_MAOM=15 #estimated but we don't know the true value
DOMinput=10 #DOM added in each addition [gC/m3], Jílková2022
fClay=0.17 #weight fraction [g/g], Jílková2022
fSilt=0.24 #weight fraction [g/g], Jílková2022
maxMAOM = 0.86 * (fClay + fSilt)*100 * BD #[gC/m3] maximum MAOM, 28208 for Jílková et al. 2022 Georgiou et al. 2022: 86 ± 9 and 48 ± 6 mg C/g silt+clay mineral for HM and LM,
maxMAOMp = maxMAOM/(MAOMratioSP+1)  # maximum primary MAOM

siltSA=45.4 #m²/kg
maxSurfaceArea=claySA*BD*fClay+siltSA*BD*fSilt #total surface area of clay and silt in m²/m³
PV=15 #!TODO volume of micropores [l/m3] but overwritten by next line
PV=np.array([45,37,37,200,6])  #pore volume for each pore size class [l/m3] / todo: estimate for our soils
PRadius=np.array([0.05,0.525,8,382.5,875])  #average radius of each pore size class [µm], defined by KEYLINK
PSA=np.zeros(5)
PSA=mf.calcPoreSurfaceArea(PV, PRadius, PSA) #pore surface area for each pore size class, calculated from  KEYLINK function
PW=np.array([45/2,37/2,37/2,200/2,6/2]) #pore water volume, assume all pores half filled , but water is in m³ while volume was in l
RootHyphaeSurface= 10000   # surface area of all roots/hyphae [m2/m3] ??unit correct / look up roots surface area equivalent to that amount of DOM input
fractionSA = RootHyphaeSurface/maxSurfaceArea #used in calcMAOM/ fraction of mineral surface area occupied by roots/hyphae
numDays=150 #number of days of incubation experiment/how long to run the model, Jílková2022
pH=4.1 #Jílková2022


Priming=1 #flag to enable Priming effect
Plotting=1 # flag 1 to enable making of plots, so that this can be turned off during sensitivity analysis etc.

# sensitivityAnalyses==False
# if sensitivityAnalyses==False
# else # sensitivity Analysis
# #run the daily calculations

#perform three experimental runs, one for each of the three treatments
#these differ in DOM input amount and CN of DOM input

DOMinput_treatments=np.array([10,10,0]) #exudates, leachates, control
CN_DOMinput_treatments = np.array([6, 80, 80])  #exudates, leachates, control, CN od control DOMinput can't be zero because of dividing by it in DOM_N calculation
treatments = np.array(["exudates", "leachates", "control"])

for i in range(len(DOMinput_treatments)):
    DOMinput = DOMinput_treatments[i]
    CN_DOMinput = CN_DOMinput_treatments[i]
    
    #initialization
    #output objects initialized as empty
    time_d=[]
    treatment=[]
    outMAOM=[]
    outMAOMp=[]
    outMAOMs=[]
    outPOM=[]
    outDOMadded=[]
    outDOM=[]
    outDOM_CN=[]
    outbact_DOM=[]
    outBact=[]
    outFungi=[]
    outRespSubstrate=[]
    outRespSoil=[]
    outRespSoilBaseline=[]

#*************************************************************************
# variables (what changes during run)

    availability=np.zeros(3)
    bact_total = 5 #total biomass of bacteria [gC/m3], final noadd average from Jílková2022
    bact = bact_total * (1-bact_DOM_rel) #biomass of bacteria growing on POM and MAOM but not on DOM [gC/m3]
    bact_sub = 0 # proportion of this bacterial carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
    bact_DOM = bact_total*bact_DOM_rel #biomass of bacteria growing on DOM [gC/m3]
    bact_DOM_sub = 0 # proportion of bacterial carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
    CN_DOM=CN_DOMinput # set inital DOM CN equal to input
    DOM=0  # DOM [gC/m3]
    DOM_sub = 0 # relative substrate derived C in DOM /values 0 to 1/, portion of DOM carbon that is substrate derived in contrast to soil-derived / values 0 to 1/ is a ratio between substrate-derived C and total C in DOM
    DOM_N=DOM/CN_DOM
    fungi=1 #biomass of fungi [gC/m3] based on final noadd in Jílková et al. 2022
    fungi_sub = 0 # proportion of fungal carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
    MAOM=25368 #C in MAOM [gC/m3] average noAdd Jílková2022 
    MAOM_sub = 0 # proportion of MAOM that is substrate derived in contrast to soil-derived / values 0 to 1/
    MAOMunavail = (PSA[0]/sum(PSA))*MAOM #the portion of MAOM stored in the smallest pores is really unavailable
    MAOMp = MAOM/(MAOMratioSP+1) #primary MAOM [gC/m3] initialised at the ratio of saturation
    MAOMs = MAOM-MAOMp #secondary MAOM[gC/m3]
    POM=13032  # C in POM [gC/m3], calculated as initialSOM-MAOM using initialSOM from Jílková2022

    # function coreMAOM
    for d in range(numDays):
        treatment.append(treatments[i])
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
        # saturation of MAOMs depends on amount of MAOMp so recalculated every day
        maxMAOMs = MAOMp*MAOMratioSP  # maximum primary MAOM
        # microbial growth on DOM and priming
        DOM, DOM_sub, DOM_N, CN_DOM, bact_DOM, bact_DOM_sub, POM, MAOM, respDOM, respDOM_sub, respPriming = mf.calcRhizosphere(Priming, POM, CN_POM, MAOM, CN_MAOM, bact_DOM, bact_DOM_sub, CN_bact, DOM, DOM_sub, CN_DOM, GMAX, DEATH, pCN, pH, rRESP, KS, DOM_EC, Priming_max, kpriming, kPOM_MAOM)
        #MAOM formation
        DOM, DOM_N, CN_DOM, MAOMp, MAOMs =mf.calcMAOM(bact_DOM, DOM_N, CN_DOM, fractionSA, MAOMp, maxMAOMp, DOM, MAOMs, maxMAOMs, MAOMsmaxrate, MAOMpmaxrate, MM_DOM_MAOM,maxEffectBactMAOM,MM_Bact_MAOM, maxEffectN_MAOM,MM_N_MAOM, maxEffectSA_MAOM,MM_SA_MAOM)
        
        # baseline microbial growth on SOM (without substrate DOM additions)
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
        #only feed on secondary MAOM
        bactPOMgrowth = mf.calcgrowth(bact, POM, availability[0], gmaxbPOM, KSbact*bact)
        bactMAOMgrowth = mf.calcgrowth(bact, MAOMs, availability[0], gmaxbMAOM, KSbact*bact)
        dbact = bactPOMgrowth + bactMAOMgrowth - DEATH*bact - rRESP*bact
        bact_sub_abs += bactMAOMgrowth*MAOM_sub - DEATH*bact*bact_sub - rRESP*bact*bact_sub #add the corresponding part of growth on MAOM as substrate derived C, subtract correspodning part of death and respiration
        
        fungiPOMgrowth = mf.calcgrowth(fungi, POM, availability[1], gmaxfPOM, KSfungi*fungi)
        fungiMAOMgrowth = mf.calcgrowth(fungi,MAOMs, availability[1], gmaxfMAOM, KSfungi*fungi)
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
        baselineRespFungi=rRESPfungi*fungi
        fungi+=dfungi
        fungi_sub= fungi_sub_abs/fungi #update relative substrate derived C in fungi       
    
        #add up substrate derived respiration
        respSubstrate = respDOM*respDOM_sub
        #add up soil-derived respiration
        baselineResp = baselineRespBact + baselineRespFungi
        respSoil = baselineResp + respPriming
        outDOMadded.append(DOM_added)
        outMAOM.append(MAOM)
        outMAOMp.append(MAOMp)
        outMAOMs.append(MAOMs)
        outPOM.append(POM)
        outDOM.append(DOM)
        outbact_DOM.append(bact_DOM)
        outBact.append(bact)
        outFungi.append(fungi)
        outRespSubstrate.append(respSubstrate)
        outRespSoilBaseline.append(baselineResp)
        outRespSoil.append(respSoil)
       
      
# plt.plot(outMAOM)      
# plt.plot(outDOM)
# plt.plot(outPOM)
# plt.plot(outbact_DOM)   
# plt.plot(outResp)
# plt.plot(outRespPriming)   
# plt.plot(outBact_Baseline) 
# plt.plot(outbact_DOM) 
# plt.plot(outFungi_Baseline) 

    #change units to easily understandable for the plot
    outDOMadded2 = np.divide(outDOMadded,0.8) #change units from gC/m3 µgC/g soil
    outbact_DOM2=np.divide(outbact_DOM,0.8) #change units from gC/m3 µgC/g soil
    outBact2=np.divide(outBact,0.8)
    outFungi2=np.divide(outFungi,0.8)
    outRespSubstrate2 =np.divide(outRespSubstrate, 0.8*24) #change units from gC/m3/day to µg CO2-C/g soil/h
    outRespSoilBaseline2=np.divide(outRespSoilBaseline, 0.8*24)
    outRespSoil2=np.divide(outRespSoil,0.8*24)
    outDOM2 = np.divide(outDOM,0.8) #change units from gC/m3 µgC/g soil
    outPOM2 = np.divide(outPOM,0.8*1000) #change units from gC/m3 mgC/g soil
    outMAOM2 = np.divide(outMAOM,0.8*1000) #change units from gC/m3 mgC/g soil
    outMAOMp2 = np.divide(outMAOMp,0.8*1000) #change units from gC/m3 mgC/g soil
    outMAOMs2 = np.divide(outMAOMs,0.8*1000) #change units from gC/m3 mgC/g soil
    
    #combine output arrays into a dataframe and save it to csv
    df = pd.DataFrame({"treatment" : treatment,
                       "DOMaddition" : outDOMadded2,
                       "DOM" : outDOM2,
                       "bact_DOM" : outbact_DOM2,
                       "bact" : outBact2,
                       "fungi" : outFungi2,
                       "resp_substrate" : outRespSubstrate2,
                       "resp_soil_baseline" : outRespSoilBaseline2,
                       "resp_soil" : outRespSoil2,
                       "POM" : outPOM2,
                       "MAOM" : outMAOM2,
                       "MAOMp" : outMAOMp2,
                       "MAOMs" : outMAOMs2
                       })
    outDataframes.append(df)
    
     
    # def Dailyplot(outDOMadded, outDOM, outbact_DOM, outRespSubstrate, outRespSoil, outRespSoilBaseline, outPOM, outMAOM): #plot in original KEYLINK units
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
    def Dailyplot2(outDOMadded2, outDOM2, outbact_DOM2, outBact2, outFungi2, outRespSubstrate2, outRespSoil2, outRespSoilBaseline2, outPOM2, outMAOMp2, outMAOMs): #plot in original KEYLINK units
        # df2 = pd.DataFrame(df)
        # x = []
        # y = []
        # fig, ((p1, p2), (p3,p4),(p5,p6)) = plt.subplots(nrows=3,
        #                                                                                   ncols=2,
        #    figsize=(5, 6))
        fig, ((p1, p2, p3), (p4, p5, p6)) = plt.subplots(nrows=2, ncols=3,figsize=(12, 7))#was 10,12
        fig.suptitle(treatments[i], size=16)
        fig.tight_layout(pad=2.0)
        ps = (p1, p2, p3, p4, p5, p6)
        # counter = count(0, 1)
        # columns = list(df)
        ps[0].set_title("DOM additions, µgC g-1 soil")
        ps[1].set_title("DOM, µgC g-1 soil") 
        ps[2].set_title("Microbial biomass, µgC g-1 soil")
        ps[3].set_title("Respiration, µg C-CO2 g-1 soil h-1")
        ps[4].set_title("POM, mgC g-1 soil")
        # ps[4].set_title("SOM, mgC g-1 soil")
        ps[5].set_title("MAOM, mgC g-1 soil")
        
        p1.plot(time_d, outDOMadded2)
        p2.plot(time_d, outDOM2)
        p3.plot(time_d, outbact_DOM2, label="bacteria DOM feeding")
        p3.plot(time_d, outBact2, label="bacteria")
        p3.plot(time_d, outFungi2, label="fungi")
        ps[2].legend(loc=(0.4, 0.03), shadow=True) #loc='bottom right',
        p4.plot(time_d, outRespSubstrate2, label="substrate-derived")
        p4.plot(time_d, outRespSoil2, label="soil-derived incl. priming")
        p4.plot(time_d, outRespSoilBaseline2, label="soil-derived baseline")
        ps[3].legend(loc=(0.25, 0.03), shadow=True) #loc='bottom right',
        p5.plot(time_d, outPOM2, label="POM")
        ps[4].legend(loc=(0.03, 0.03), shadow=True) #loc='upper left',
        p6.plot(time_d, outMAOMp2, label="primary MAOM")
        p6.plot(time_d, outMAOMs2, label="secondary MAOM")    
        ps[5].legend(loc=(0.4, 0.03), shadow=True) #loc='bottom right',
    
    # Dailyplot(outDOMadded, outDOM, outbact_DOM, outRespSubstrate, outRespSoil, outRespSoilBaseline, outPOM, outMAOM)
    if Plotting == 1: #if you want plotting to be active
        Dailyplot2(outDOMadded2, outDOM2, outbact_DOM2, outBact2, outFungi2, outRespSubstrate2, outRespSoil2, outRespSoilBaseline2, outPOM2, outMAOMp2, outMAOMs2)
        plt.savefig("output/figures/Dailyplot_" + treatments[i] +".png")

#after running the outermost loop (for three different treatments)
#merge the three dataframes to create a data output containing all three treatments
dfAll=pd.concat(outDataframes)
dfAll.to_csv(".\output\data\Output.csv", index=False)

