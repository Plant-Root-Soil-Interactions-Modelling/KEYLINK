# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:08:20 2023

@author: gdeckmyn
"""
import keylink_functions as mf
import matplotlib.pyplot as plt
import numpy as np

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
DOM_RS=0  # rhizosphere DOM [gC/m3]
bact_RS=61*0.2 #bacterial biomass [gC/m3], noadd average from Jílková2022, was 2.4, *0.2 because most was li-ving on SOM, only 20% on new C
CN_DOM_RSinput=80  #CN of the daily input [unitless], Jílková2022: leachates 80, exudates 6, was 50
MAOMsaturation=0.9 #proportion of MAOM capacity filled [unitless], based on MAOM content from Jílková2022, was 0.5
maxMAOM=  28208 #MAOM capacity [gC/m3] based on silt+clay from Jílková2022 and Georgiou et al. 2022, was 100
MAOMini=MAOMsaturation*maxMAOM #initial MAOM concentration [gC/m3]
MM_DOMtoMAOM=0.025  # DOM concentration for speed being half max speed (Michaelis Menten) [gC/m3]
MAOMmaxrate=0.2 # max proportion of DOM stabilized in MAOM per day [unitless]
surface_RS= 10000   # surface area of all roots/hyphae [m2/m3] ??unit correct
DOMinput=10 #DOM added in each addition [gC/m3], Jílková2022, was 0.5
numDays=150 #number of days of incubation experiment/how long to run the model, Jílková2022
POMini=38400  # total SOM, only used for priming [gC/m3], Jílková2022, was 150
resp=0 #respiration [gC/m3]
GMAX=1.24 #growth rate [gC/(gC day)], KEYLINK
DEATH=0.05 #death rate [gC/(gC day)], KEYLINK
rRESP=0.05  #respiration rate resp, [gC/(gC day)], KEYLINK
KS=5  # concentration of ?substrate (DOM) for half speed growth, for growing on DOM [gC/m3], related to substrate quality was 0.001
MCN=0.8 #???, KEYLINK
CN_bact=4 #CN of bacteria, Jílková2022 initial is 10
pH=4.1 #Jílková2022, was 3.5
Nmin=0.00000001 # was 0.0005 [gN/m3] was 0.00000001
POMCN=24 #CN of SOM, Jílková2022, was 20
PV=15 # volume of micropores [l/m3]
CN_fungi=8 #KEYLINK
GMAXfungi=0.6  #growth rate [gC/(gC day)], KEYLINK
rRESPfungi=0.03 #respiration rate resp, [gC/(gC day)], KEYLINK
DEATHfungi=0.02 #death rate [gC/(gC day)], KEYLINK
primingIntensity=5 #ratio of POM decayed for DOM decayed [gC/gC, unitless], depends on DOM quality, we know DOM CN which is something else
mRecBact=0.5  # how sensitive bact are to recalcitrance
mRecFungi=0.5
recMAOM= 0.9
MAOM_CN=15
claySA=800000 #m²/kg was 8000000 cm²/g
siltSA=45.4 #m²/kg
fSilt=0.24 #weight fraction
fCLay=0.17
BD=800   # bulk density kg m³
CN_DOM_RS=CN_DOM_RSinput # set inital DOM CN equal to input
POM=POMini
MAOM=MAOMini
DOM_N=DOM_RS/CN_DOM_RS
bact=61 * 0.8 #80% of my bacteria are doing the baseline
fungi=4
KSfungi=20000  # for decaying SOM
KSbact=38000 #for decaying SOM
PV=np.array([45,37,37,200,6]) 
PRadius=np.array([0.05,0.525,8,382.5,875])  #in µm
PSA=np.zeros(5)
PW=np.array([45/1000,37/1000,37/1000,200/1000,6/1000]) #☺assume all pores filled , but water is in m³ while volume was in l
SAclaySilt=claySA*BD*fCLay+siltSA*BD*fSilt #total surface area of clay and silt in m²/m³
PSA=mf.calcPoreSurfaceArea(PV, PRadius, PSA)
availability=np.zeros(3)
MAOMunavail = (PSA[0]/sum(PSA))*MAOMini #MAOM stored in the smallest pores is really unavailable, 

#on day 0 and then every 14 days, add DOM
for d in range(numDays):
      time_d.append(d)  #store days in an array for plotting
      DOM_added = 0
      if d==0 or (d%14)==0: #where does this if end?
          DOM_added=DOMinput #to keep track of the additions
          DOM_RS+=DOMinput 
          DOM_N+=DOMinput/CN_DOM_RSinput
      CN_DOM_RS=DOM_RS/ DOM_N
  # growth in soil matrix 
      MAOMunavail = max((PSA[0]/sum(PSA))*MAOM,MAOMunavail)   #unavail can only go up unless you disrupt
      availability=mf.calcAvailPot(PV, PW)
  
      gmaxbPOM = mf.calcgmaxmod(CN_bact, POMCN, MCN, 0.0, 0, pH, 1)*GMAX #gmax for bact on POM
  #    gmaxflit = mf.calcgmaxmod(CN[1], litterCN, MCN[1], recLit, MREC[1], pH, 2)* GMAX[1] #gmax for fung on litter
      gmaxfPOM = mf.calcgmaxmod(CN_fungi, POMCN, MCN, 0.0, 0, pH, 2)*GMAXfungi #gmax for fung on SOM
      gmaxbMAOM = mf.calcgmaxmod(CN_bact, MAOM_CN, MCN, recMAOM, mRecBact, pH, 1)*GMAX #gmax for bact on POM
 #    gmaxflit = mf.calcgmaxmod(CN[1], litterCN, MCN[1], recLit, MREC[1], pH, 2)* GMAX[1] #gmax for fung on litter
      gmaxfMAOM = mf.calcgmaxmod(CN_fungi, MAOM_CN, MCN, recMAOM, mRecFungi, pH, 2)*GMAXfungi #gmax for fung on SOM
    

      #growth equations (dB/dt) for each functional group and for variations in C pools
      dbact = mf.calcgrowth(bact, POM, availability[0], gmaxbPOM, KSbact)+ \
              mf.calcgrowth(bact, MAOM-MAOMunavail, availability[0], gmaxbMAOM, KSbact)+ \
              - DEATH*bact - rRESP*bact
      dfungi = mf.calcgrowth(fungi, POM, availability[1], gmaxfPOM, KSfungi) - DEATHfungi*fungi - rRESPfungi*fungi \
               + mf.calcgrowth(fungi,MAOM-MAOMunavail, availability[1], gmaxfMAOM, KSfungi) - DEATHfungi*fungi - rRESPfungi*fungi
      POM+=-mf.calcgrowth(bact, POM, availability[0], gmaxbPOM, KSbact)-mf.calcgrowth(fungi, POM, 1, gmaxfPOM, KSfungi)  
      DOM_RS+=DEATH*bact+DEATHfungi*fungi
      print(DOM_RS, bact)
      DOM_N+=DEATH*bact/CN_bact+DEATHfungi*fungi/CN_fungi
      MAOM+=-mf.calcgrowth(bact, MAOM-MAOMunavail, availability[0], gmaxbMAOM, KSbact)-   \
          mf.calcgrowth(fungi,MAOM-MAOMunavail, availability[1], gmaxfMAOM, KSfungi)
      baselineRespBact=rRESP*bact
      bact+=dbact
      baselineRespFungi=rRESP*fungi
      fungi+=dfungi
              

  # growth in rhizosphere    
      DOM_RS,DOM_N, bact_RS, POM, resp, respPriming= mf.calcRhizosphere(MAOMsaturation,maxMAOM, bact_RS, DOM_RS, GMAX, DEATH,CN_bact, CN_DOM_RS, MCN, pH, rRESP, KS, POMCN, Nmin, POM, PV, primingIntensity, MAOM_CN)
                                                         # ((MAOMsaturation,maxMAOM,bact_RS, DOM_RS, gmax, DEATH,CN_bact, CN_DOM_RS, pCN, pH, res, Ks, fCN, CN_SOM, Nmin, SOM,PVstruct,  primingIntensity)        # MAOM formation
      MAOMsaturation, DOM_RS,DOM_N=mf.calcMAOMsaturation (maxMAOM,MM_DOMtoMAOM, MAOMsaturation, MAOMmaxrate, DEATH, SAclaySilt,bact_RS, surface_RS, DOM_RS, DOM_N, CN_DOM_RS)
      MAOM=MAOMsaturation*maxMAOM   
      #add up soil-derived respiration
      baselineResp = baselineRespBact + baselineRespFungi
      respSoil = baselineResp + respPriming
      outDOMadded.append(DOM_added)
      outMAOM.append(MAOM)
      outMAOMsaturation.append(MAOMsaturation)
      outPOM.append(POM)
      outDOM.append(DOM_RS)
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
outDOM2 = np.divide(outDOM,0.8) #change units from gC/m3 µgC/g soil
outBact_RS2=np.divide(outBact_RS,0.8) #change units from gC/m3 µgC/g soil
outBact_Baseline=np.divide(outBact_Baseline,0.8)
outFungi_Baseline=np.divide(outBact_Baseline,0.8)
outRespSubstrate2 =np.divide(outRespSubstrate, 0.8*24) #change units from gC/m3/day to µg CO2-C/g soil/h
outRespSoilBaseline2=np.divide(outRespSoilBaseline, 0.8*24)
outRespSoil2=np.divide(outRespSoil,0.8*24)
outDOM2 = np.divide(outDOM,0.8) #change units from gC/m3 µgC/g soil
outPOM2 = np.divide(outPOM,0.8*1000) #change units from gC/m3 mgC/g soil
outMAOM2 = np.divide(outMAOM,0.8*1000) #change units from gC/m3 mgC/g soil

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
    
#     p1.plot(time_d, outDOMadded)
#     p2.plot(time_d, outDOM)
#     p3.plot(time_d, outBact_RS, label="bacterial biomass")
#     p4.plot(time_d, outRespSubstrate, label="substrate-derived")
#     p4.plot(time_d, outRespSoil, label="soil-derived incl. priming")
#     p4.plot(time_d, outRespSoilBaseline, label="soil-derived baseline")
#     ps[3].legend(loc=(0.1, 0.7), shadow=True) #loc='upper left',
#     p5.plot(time_d, outPOM)
#     p6.plot(time_d, outMAOM)
    # plt.legend(loc=(1.01, 0), shadow=True) #loc='upper right',

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