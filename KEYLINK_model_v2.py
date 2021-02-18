'''
KEYLINK MODEL

Created on 27.07.2016 last write 18/02/2021

@author: A Schnepf - G Deckmyn - G Cauwenberg - O Flores
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import functions_KEYLINK_v2 as mf
import random
import matplotlib
import math
import string

def import_pools(filename):
    """Load array from text file"""
    return np.loadtxt(filename + '.txt')

def export_pools(filename, array):
    """Save an array in a readable format, compatible with R"""
    np.savetxt(filename + '.txt', array)

pfaec = np.zeros(7)
PVstruct=np.zeros(5)
drainage = 0
runoff = 0

B = import_pools('KL_initC_Pools') #initial biomass in each C pool
(GMAX, KS, DEATH, RESP, FAEC, CN, REC, MCN,
 MREC, T_MIN, T_OPT, T_MAX, Q10) = import_pools('KL_FaunalParams') #parameters for each functional group

iniSOM=B[10] #initial SOM in g C/m3

B.resize((22,))

(TEMP, SUNH) = import_pools('KL_climateParams') #monthly climate data
d, BD, alpha, n, m, Ksat, pH, litterCN, SOMCN, drainmax, PVstruct[0], PVstruct[1], PVstruct[2], PVstruct[3], PVstruct[4] = import_pools('KL_initSoil') #soil parameters
ratioPVBeng, fPVB, tPVB,  PVBmax, frag, pfaec[5], pfaec[6], bioturbRate, moveRate= import_pools('KL_engineerParams') #parameters for engineer activity

# Parameters CNlit=of daily litter, litterCN=total litter pool
tStop, initWater, Nmin, rrg, rootTO, inputLit, CNlit, recLit, CtoMyc, NmyctoPlant, ee = import_pools('KL_runparams')
PW = initWater/100*PVstruct #fraction of pore volume filled with water

Nfauna=sum(B[:9]/CN) #N in food web functional groups

Ntot=Nfauna+Nmin+B[10]/SOMCN+B[9]/litterCN #initial total N of the system

gr = np.array([4000,2000,500,50,0.1]) #graph ranges, for five biomass graphs with ranges from 0 to these values
#graph labels for the populations and C pools
popname = np.array(['Bacteria','Fungi','Mycorrhiza','bacterivores','fungivores','saprotrophs','engineers','herbivores','predators','litter','SOM','roots','CO2'])
#graph colours for the populations and C pools
popcolour = np.array(['blue','red','darkcyan','cyan','orange','purple','darkgreen','magenta','black','brown','grey','chartreuse','yellow'])


'''
    Plot population size as a function of time
'''
def show_plot(soln, pwt, pvt):
    
    plt.clf()
    
    plt.subplot(2, 1, 1)
    plt.plot(t, pwt, label="SW")
    plt.plot(t, pvt, label="porosity")

    plt.ylabel('Water, l/m2')
    plt.ylim(0,400)
    plt.legend(loc=(1.01, 0), shadow=True) #loc='upper right',
    
    plt.subplot(2, 1, 2)
    for p in range(13): # it only shows populations with max biomass value higher than gr[1]
        if max(soln[:,p])>gr[1]:
            plt.plot(t, soln[:, p], label=popname[p], c=popcolour[p])
    
    plt.xlabel('Time,days')
    plt.ylabel('Biomass, gC/m3')
    plt.ylim(0, gr[0])
    plt.legend(loc=(1.01, 0), shadow=True)
    
    plt.show()
    
    plt.clf()
    
    plt.subplot(2, 1, 1)
    for p in range(13): # shows populations with maximum or mean biomass values between gr[1] and gr[2]
        maxvalue=(max(soln[:,p]))
        meanvalue=((sum(soln[:,p]))/tStop)
        if ((maxvalue>gr[2]) and (maxvalue<gr[1])) or (meanvalue>gr[2] and meanvalue<gr[1]):
            plt.plot(t, soln[:, p], label=popname[p], c=popcolour[p])    
    
    plt.ylabel('Biomass, gC/m3')
    plt.ylim(0, gr[1])
    plt.legend(loc=(1.01, 0), shadow=True)
    
    plt.subplot(2, 1, 2)
    for p in range(13): # shows populations with maximum or mean biomass values between gr[2] and gr[3]
        maxvalue=(max(soln[:,p]))
        meanvalue=((sum(soln[:,p]))/tStop)
        if ((maxvalue>gr[3]) and (maxvalue<gr[2])) or (meanvalue>gr[3] and meanvalue<gr[2]):
            plt.plot(t, soln[:, p], label=popname[p], c=popcolour[p])    
    
    plt.xlabel('Time,days')
    plt.ylabel('Biomass, gC/m3')
    plt.ylim(0, gr[2])
    plt.legend(loc=(1.01, 0), shadow=True)
    
    plt.show()    
    
    plt.clf()
    
    plt.subplot(2, 1, 1)
    for p in range(13): # shows populations with maximum or mean biomass values between gr[3] and gr[4]
        maxvalue=(max(soln[:,p]))
        meanvalue=((sum(soln[:,p]))/tStop)
        if ((maxvalue>gr[4]) and (maxvalue<gr[3])) or (meanvalue>gr[4] and meanvalue<gr[3]):
            plt.plot(t, soln[:, p], label=popname[p], c=popcolour[p])
    
    plt.ylabel('Biomass, gC/m3')
    plt.ylim(0, gr[3])
    plt.legend(loc=(1.01, 0), shadow=True)
    
    plt.subplot(2, 1, 2)
    te=(tStop-int(round((tStop/5))))
    # shows populations with last 20% days mean biomass values lower than gr[4] (assumed local extinction)
    for p in range(13):
        extvalue=((sum(soln[te:tStop,p]))/(tStop-te))
        if (extvalue<gr[4]):
            plt.plot(t, soln[:, p], label=popname[p], c=popcolour[p])    
    
    plt.xlabel('Time,days')
    plt.ylabel('Biomass, gC/m3')
    plt.ylim(0, gr[4])
    plt.legend(loc=(1.01, 0), shadow=True)

    plt.show()


HI=0 #heat index (is calculated with the daily mean temperature of all months)
for m in range(12):
    HI = HI + (TEMP[m]/5)**1.514
alfa = 0.000000675*(HI**3) - 0.0000771*(HI**2) + 0.01792*HI + 0.49239 #alpha for Thornthwaite equation (for evapotranspiration)


# function to be integrated daily solving the carbon pools 'B' ifo time
def f(B, t, avail, modt, GMAX, litterCN,SOMCN):
    (availSOMbact, availSOMfungi, availSOMeng, availSOMsap, availbbvores,
     availffvores, availfvorespred, availbvorespred, availhvorespred,
     availsappred, availengpred ,SOMunavail) = avail

    # 0=bact, 1=fungi, 2=myc, 3=bvores, 4=fvores, 5=sap,
    # 6=eng, 7=hvores, 8=pred, 9=litter, 10=SOM, 11=roots, 12=CO2

    # update GMAX for bacteria, fung and myc GMAX is modified for SOM
    # and litter seperately depending on CN (and possibly recalcitrance)
    #for bact if CN source too high they can't grow   
    gmaxblit = mf.calcgmaxmod(CN[0], litterCN, MCN[0], recLit, MREC[0], pH, 1)*GMAX[0] #gmax for bact on litter
    gmaxbSOM = mf.calcgmaxmod(CN[0], SOMCN, MCN[0], 0.0, MREC[0], pH, 1)*GMAX[0] #gmax for bact on SOM
    gmaxflit = mf.calcgmaxmod(CN[1], litterCN, MCN[1], recLit, MREC[1], pH, 2)* GMAX[1] #gmax for fung on litter
    gmaxfSOM = mf.calcgmaxmod(CN[1], SOMCN, MCN[1], 0.0, MREC[1], pH, 2)* GMAX[1] #gmax for fung on SOM
    gmaxEng = min(mf.calcgmaxEng(GMAX[6],pH),GMAX[6]) #gmax for engineers
    
    #update faeces for  SAP and engineers
    faeclitEng = min(1,mf.calcFaec(gmaxEng, FAEC[6], pfaec[6], litterCN, CN[6], rRESP[6]))
    faeclitSAP = min(1,mf.calcFaec(GMAX[5], FAEC[5], pfaec[5], litterCN, CN[5], rRESP[5]))

    #growth equations for each functional group and for variations in C pools
    bact = (modt[0]*(mf.calcgrowth(B[0], B[10]-SOMunavail, availSOMbact, gmaxbSOM, KS[0])
            + mf.calcgrowth(B[0], B[9], availSOMbact, gmaxblit, KS[0]))
            - DEATH[0]*B[0] - rRESP[0]*B[0]
            - modt[3]*mf.calcgrowth(B[3], B[0], availbbvores, GMAX[3], KS[3]))

    fungi = (modt[1]*(mf.calcgrowth(B[1], B[10]-SOMunavail, availSOMfungi, gmaxfSOM, KS[1])
             + mf.calcgrowth(B[1], B[9], availSOMbact, gmaxflit, KS[1]))
             - DEATH[1]*B[1] - rRESP[1]*B[1]
             - modt[4]*mf.calcgrowth(B[4], B[1], availffvores, GMAX[4], KS[4]))

    myc = (mf.inputCtoMyc(CtoMyc)
           + modt[2]*(mf.calcgrowth(B[2], B[9], availSOMbact, gmaxflit, KS[2])
           + mf.calcgrowth(B[2], B[10]-SOMunavail, availSOMfungi, gmaxfSOM, KS[2]))
           - DEATH[2]*B[2] - rRESP[2]*B[2]
           - modt[4]*mf.calcgrowth(B[4], B[2], availffvores, GMAX[4], KS[4]))
           # myc (being a fungi) has the same availability and gmax than fungi
    bvores = (modt[3]*mf.calcgrowth(B[3], B[0], availbbvores, GMAX[3], KS[3])
              - modt[8]*(1+FAEC[8])*mf.calcgrowth(B[8], B[3], availbvorespred, GMAX[8], KS[8])
              - DEATH[3]*B[3] - rRESP[3]*B[3])

    fvores = (modt[4]*(mf.calcgrowth(B[4], B[1], availffvores, GMAX[4], KS[4])
              + mf.calcgrowth(B[4], B[2], availffvores, GMAX[4], KS[4]))
              - modt[8]*(1+FAEC[8])*mf.calcgrowth(B[8], B[4], availfvorespred, GMAX[8], KS[8])
              - DEATH[4]*B[4] - rRESP[4]*B[4])

    sap = (modt[5]*(mf.calcgrowth(B[5], B[9],availSOMbact, GMAX[5], KS[5])
           + mf.calcgrowth(B[5], B[10]-SOMunavail, availSOMsap, GMAX[5], KS[5]))
           -modt[8]*(1+FAEC[8])*mf.calcgrowth(B[8], B[5], availsappred, GMAX[8], KS[8])
           - DEATH[5]*B[5] - rRESP[5]*B[5])

    eng = (modt[6]*(mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6])
           + mf.calcgrowth(B[6], B[10]-SOMunavail, availSOMeng, gmaxEng, KS[6]))
           - modt[8]*(1+FAEC[8])*mf.calcgrowth(B[8], B[6], availengpred, GMAX[8], KS[8])
           - DEATH[6]*B[6] - rRESP[6]*B[6])

    #roots are avaialble because larger than herbivores
    hvores = (modt[7]*mf.calcgrowth(B[7], B[11], 1, GMAX[7], KS[7])
              - modt[8]*(1+FAEC[8])*mf.calcgrowth(B[8], B[7], availhvorespred, GMAX[8], KS[8])
              - DEATH[7]*B[7] - rRESP[7]*B[7])

    pred = (modt[8]*(mf.calcgrowth(B[8], B[3], availbvorespred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[4], availfvorespred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[5], availsappred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[6], availengpred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[7], availhvorespred, GMAX[8], KS[8]))
            - DEATH[8]*B[8] - rRESP[8]*B[8])

    litter = (           
              -modt[0]*mf.calcgrowth(B[0], B[9], availSOMbact, gmaxblit, KS[0]) #eaten by bact
              -modt[1]*mf.calcgrowth(B[1], B[9], availSOMbact, gmaxflit, KS[1]) #eaten by fungi
              -modt[2]*mf.calcgrowth(B[2], B[9], availSOMbact, gmaxflit, KS[2]) #eaten by myc
              -modt[5]*(1+faeclitSAP)*mf.calcgrowth(B[5], B[9], availSOMbact, GMAX[5], KS[5]) #eaten by SAP
              -modt[6]*(1+faeclitEng)*mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6]) # eaten by engineers
              + DEATH[5]*B[5] + DEATH[6]*B[6]+ DEATH[7]*B[7] + DEATH[8]*B[8])

    som = (mf.exudation() 
           - modt[0]*mf.calcgrowth(B[0], B[10]-SOMunavail, availSOMbact, gmaxbSOM, KS[0])
           - modt[1]*mf.calcgrowth(B[1], B[10]-SOMunavail, availSOMfungi, gmaxfSOM, KS[1])
           - modt[2]*mf.calcgrowth(B[2], B[10]-SOMunavail, availSOMfungi, gmaxfSOM, KS[2])
           - modt[5]*mf.calcgrowth(B[5], B[10]-SOMunavail, availSOMsap, GMAX[5], KS[5]) #eaten by SAP
           - modt[6]*mf.calcgrowth(B[6], B[10]-SOMunavail, availSOMeng, gmaxEng, KS[6]) # eaten by engineers
           + modt[8]*FAEC[8] * (mf.calcgrowth(B[8], B[3], availbvorespred, GMAX[8], KS[8])
                        + mf.calcgrowth(B[8], B[4], availfvorespred, GMAX[8], KS[8])
                        + mf.calcgrowth(B[8], B[5], availsappred, GMAX[8], KS[8])
                        + mf.calcgrowth(B[8], B[6], availengpred, GMAX[8], KS[8])
                        + mf.calcgrowth(B[8], B[7], availhvorespred, GMAX[8], KS[8]))
           + modt[5]*faeclitSAP*mf.calcgrowth(B[5], B[9], availSOMbact, GMAX[5], KS[5])
           + modt[6]*faeclitEng*mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6])
           + modt[7]*FAEC[7]*mf.calcgrowth(B[7], B[11], 1, GMAX[7], KS[7])
           + DEATH[0]*B[0]+DEATH[1]*B[1]+DEATH[2]*B[2]+DEATH[3]*B[3]+DEATH[4]*B[4])

    roots = (- modt[7]*(1+FAEC[7])*mf.calcgrowth(B[7], B[11], 1, GMAX[7], KS[7]))
   
    co2 = (rRESP[0]*B[0]+rRESP[1]*B[1]+rRESP[2]*B[2]+rRESP[3]*B[3] #CO2 emissions from respiration
           +rRESP[4]*B[4]+rRESP[5]*B[5]+rRESP[6]*B[6]+rRESP[7]*B[7]+rRESP[8]*B[8])
           
    bactResp=rRESP[0]*B[0] #respiration of bacteria
    funResp=rRESP[1]*B[1] #respiration of fungi
    EMresp=rRESP[2]*B[2] #respiration of mycorrhizal fungi
    bactGrowthSOM=modt[0]*mf.calcgrowth(B[0], B[10]-SOMunavail, availSOMbact, gmaxbSOM, KS[0]) #growth of bact from eaten SOM
    bactGrowthLit=modt[0]*mf.calcgrowth(B[0], B[9], availSOMbact, gmaxblit, KS[0]) #growth of bact from eaten litter
    SOMeaten=modt[0]*mf.calcgrowth(B[0], B[10]-SOMunavail, availSOMbact, gmaxbSOM, KS[0]) #SOM eaten by bact
    + modt[1]*mf.calcgrowth(B[1], B[10]-SOMunavail, availSOMfungi, gmaxfSOM, KS[1]) #eaten by fungi
    + modt[2]*mf.calcgrowth(B[2], B[10]-SOMunavail, availSOMfungi, gmaxfSOM, KS[2]) #eaten by myc
    + modt[5]*mf.calcgrowth(B[5], B[10]-SOMunavail, availSOMsap, GMAX[5], KS[5]) #eaten by SAP
    + modt[6]*mf.calcgrowth(B[6], B[10]-SOMunavail, availSOMeng, gmaxEng, KS[6]) # eaten by engineers
    LITeaten=modt[0]*mf.calcgrowth(B[0], B[9], availSOMbact, gmaxblit, KS[0]) #Litter eaten by bact
    +modt[1]*mf.calcgrowth(B[1], B[9], availSOMbact, gmaxflit, KS[1]) #eaten by fungi
    +modt[2]*mf.calcgrowth(B[2], B[9], availSOMbact, gmaxflit, KS[2]) #eaten by myc
    +modt[5]*(1+faeclitSAP)*mf.calcgrowth(B[5], B[9], availSOMbact, GMAX[5], KS[5]) #eaten by SAP
    +modt[6]*(1+faeclitEng)*mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6]) # eaten by engineers      
    LITeatenEng=modt[6]*(1+faeclitEng)*mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6]) #only litter eaten by enginners
 
    return [bact, fungi, myc, bvores, fvores, sap,
            eng, hvores, pred, litter, som, roots, co2,
            bactResp,funResp,EMresp,bactGrowthSOM,bactGrowthLit, SOMeaten, LITeaten, LITeatenEng]   

    
PVt = np.zeros(tStop)
PWt = np.zeros(tStop)
psoln = np.zeros((tStop, 22))
avail = np.zeros(11)
modt = np.zeros(9)
rRESP = np.zeros(9)

std=0 #starting day
pv= np.zeros(5)
pv=PVstruct #initialise to structural

climatefile=open('PrecipKMIBrass.txt') #meteorology data
titls=climatefile.readline() # first line are titles

# this is the actual core model routine over time steps i
for i in range(int(std), int(std+tStop)):
    # calculate PSD (array of % of five size classes of pores) and aggregation
    ag = mf.calcAg(B[1], B[2], B[10]) #calculates aggregation fraction
    pv = mf.calcPVD(PVstruct, pv, ag, ratioPVBeng, fPVB, tPVB, PVBmax, d, B)
    pvd = pv*100/sum(pv)
    
    nd = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) # number of days in each month
    
    zip=climatefile.readline().split() #daily meteorological data
    year=int(zip[1])
    month=(int(zip[2])-1)
    ly=0 #it will became 1 in leap-years
    if (year%4)==0:
        ly=1
    if (year%100)==0:
        ly=0
    if (year%400)==0:
        ly=1
        
    if ly==1: #in leap-years, this changes the number of days of february to 29
        nd[1]=nd[1]+1

    precip=float(zip[4])
    temp=float(zip[5])
    for j in range (9):
        modt[j] = mf.calcmodt(temp, T_OPT[j], T_MIN[j], T_MAX[j])
        rRESP[j]=mf.calcresp(temp, T_OPT[j], RESP[j], Q10[j])
    
    #calculate the Potential Evapotranspiration (pet)
    pet = mf.PET(temp, nd[month], SUNH[month], HI, alfa)

    # water calculations                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    # (pvd (pore volume density = array of % filled of each porse size class)
    SatW = sum(pv[0:4])/(d-(pv[4])/sum(pv[0:5])) #saturated water content (excluding macropores)
    PW, drainage, runoff = mf.calcPW(pv, precip, PW, drainmax,d, (pv[0]/d), SatW, alpha, n, m, Ksat)
   
    PWt[i-std] = sum(PW) #soil water content in the current day
    PVt[i-std] = sum(pv) #total soil porosity in the current day

    #calculate availability of each pool to each relevant biota
    avail = mf.calcAvail(pv, PW, iniSOM)

    #update litter CN
    litterCN=(B[9]+mf.inputLitter(inputLit, CNlit)[0]+mf.rootTurnover(B[11],rootTO))/(B[9]/litterCN+(mf.inputLitter(inputLit, CNlit)[0]+mf.rootTurnover(B[11],rootTO))/CNlit)
    
    #add litter (plant and root) and update total soil N
    B[9]=B[9]+mf.inputLitter(inputLit, CNlit)[0]+mf.rootTurnover(B[11],rootTO)
    Ntot=Ntot+(mf.inputLitter(inputLit, CNlit)[0]+mf.rootTurnover(B[11],rootTO))/CNlit
    
    #root growth
    B[11]=B[11]+mf.rootgrowth(rrg)-mf.rootTurnover(B[11],rootTO)
    
    #interaction between mycorrhizal fungi and plants
    B[2]=B[2]+mf.inputCtoMyc(CtoMyc)
    
    #update soil C and N in soil
    CNsoil=(B[9]+B[10])/(B[9]/litterCN+B[10]/SOMCN)
    Ntot=Ntot-mf.mycNtoPlant(NmyctoPlant, CtoMyc, litterCN, CtoMyc, CNsoil)    
    
    #plant N uptake    
    Ntot=Ntot-min(mf.plantNuptake(litterCN, inputLit, NmyctoPlant), Nmin)
    Nmin=Nmin-min(mf.plantNuptake(litterCN, inputLit, NmyctoPlant), Nmin)
    
    # Call the ODE solver for day i
#    day = odeint(f, B, [i, i+1], args=(avail,modt,))
    day = odeint(f, B, [i, i+1], args=(avail, modt, GMAX, litterCN, SOMCN))

    # Second column is end value for day i, start value for day i + 1
    psoln[i-std] = day[1, :]
    B = day[1, :]
    et = ee*pet #evapotranspiration (rate of effective evapotranspiration * potential evapotranspiration)
    PW = mf.wl(PW,et) # water lost by evapotranspiration: we do this at the end of the day (otherwise soil is always dry)

    # close N budget by adding up all N and putting 'restvalue' in SOM but not part by bact (goes into mineralised)
    Nmin=Nmin-(B[16]+B[17]-B[13])/CN[0]+B[16]/litterCN+B[17]/SOMCN    # adding bact resp - growth /CN for lit and SOM
    if Nmin<0:   #Bact use more N then they 'eat' this needs to come from somewhere so from SOM (but needs to be corrected)
        Nneg=Nmin
        Nmin=0
    else:
        Nneg=0
    Nfauna=sum(B[0:9]/CN)
    NSOM=Ntot-Nfauna-Nmin+Nneg
    SOMCN=B[10]/NSOM

    
    #move SOM by engineers
    SOMdown=mf.calcBioturb(B[6], bioturbRate, B[10])
    B[10]=B[10]-SOMdown
    Ntot=Ntot-SOMdown/SOMCN

    #move litter by engineers
    Litdown=mf.calcLittermove(B[6], moveRate, B[9])
    B[9]=B[9]-Litdown
    Ntot=Ntot-Litdown/litterCN
    
    #fragmentation    
    B[9]=B[9]-frag*psoln[i,21]
    B[10]=B[10]+frag*psoln[i,21]  #LITeatenEng
     
    for s in range (0, 11):
        if B[s]<=0: #security codes to avoid errors by negative biomasses
            B[s]=0.001
    
climatefile.close

# save output
export_pools('keylinkoutput', psoln)
t = np.arange(0., tStop)

show_plot(psoln, PWt, PVt)