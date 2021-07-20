"""
FUNCTIONS FOR KEYLINK MODEL

Created on 01.06.2017 last write 24/05/2021

@author: A Schnepf - G Deckmyn - G Cauwenberg - O Flores
"""

import numpy as np
import random
import math
from scipy.integrate import odeint
import time

def calcKD(pH,fClay): #Cempirically from orchideeSOM model Cammino-serrano et al 2018
    Kd=0.001226-0.000212*pH+0.00374*fClay
    return Kd
    
def calcAvail(PV, PW, SOMini):
    #units in water and porosity volumes: l/m3
    mwater = np.zeros(5)
    if sum(PW)/sum(PV) < 0.5:
        mwatertot = 4*sum(PW)/sum(PV)*(1-sum(PW)/sum(PV))
    else:
        mwatertot = 1
    for i in range(4):
        if PW[i]/PV[i] == 1:   # pore size that is filled
            mwater[i] = 0
        elif PW[i]/PV[i] > 0:
            mwater[i] = mwatertot*PV[i]/(PV[i]+PW[i])
            if(i+1) < 5:
                mwater[i+1] = mwatertot*PW[i]/(PV[i]+PW[i])
                
    availSOMbact    = np.sum([1, 1, 1, 1]*PV[1:]*mwater[1:]/sum(PV[1:])) #SOM availability to bact
    availSOMfungi   = np.sum([0, 1, 1, 1]*PV[1:]*mwater[1:]/sum(PV[1:])) #SOM availability to fung
    availSOMeng     = np.sum([1, 1, 1, 1]*PV[1:]*mwater[1:]/sum(PV[1:])) #SOM availability to engineers
    availSOMsap     = np.sum([0, 0, 1, 1]*PV[1:]*mwater[1:]/sum(PV[1:])) #SOM availability to detritivores
    availbbvores    = np.sum([0, 1, 1, 1]*PV[1:]*mwater[1:]/sum(PV[1:])) #bact availability to bacterivores
    availffvores    = np.sum([1, 1, 1]*PV[2:]*mwater[2:]/sum(PV[2:])) #fung availability to fungivores
    availfvorespred = np.sum([0, 1, 1]*PV[2:]*mwater[2:]/sum(PV[2:])) #fungivores availability to predators
    availbvorespred = np.sum([0, 1, 1]*PV[2:]*mwater[2:]/sum(PV[2:])) #bacterivores  availability to predators
    availhvorespred = np.sum([0, 1, 1]*PV[2:]*mwater[2:]/sum(PV[2:])) #herbivores  availability to predators
    availsappred    = np.sum([1, 1]*PV[3:]*mwater[3:]/sum(PV[3:])) #sap  availability to predators (sap only in larger pores)
    availengpred    = 1 #engineers availability to predators (earthworms can't hide)
    SOMunavail      = (PV[0]/sum(PV))*SOMini #SOM in inaccesible pores
    availability    = np.array([availSOMbact, availSOMfungi, availSOMeng,
                                availSOMsap, availbbvores, availffvores,
                                availfvorespred, availbvorespred,
                                availhvorespred, availsappred, availengpred, SOMunavail])
    return availability


def inputLitter(inLit, CNlit):  # input from plant / litter to soil surface
    value = [inLit,CNlit]
    return value

def inputCtoMyc(inC):     # input from plant / litter to mycorrhizal fungi
    value = inC
    return np.real(value)
    
def mycNtoPlant(NmyctoPl, inlit, CNlit, inC, CNsoil): #N through the myc to the plant
    value=min(NmyctoPl*inlit/CNlit, 5*inC/CNsoil) #% of plant requirement is traded by myc
    return np.real(value)

def rootgrowth(rrg):
    value = rrg
    return np.real(value)
    
def rootTurnover(Broot,TurnoverRate):
    deadRoot=Broot*TurnoverRate
    return deadRoot

def plantNuptake(CNlit, inlit, NmyctoPl):
    plantNup=(1-NmyctoPl)*inlit/CNlit
    return plantNup
    
def calcAg(fungi, myc, SOM):
    '''
    Fraction aggregation of SOM in function of the fungal and EM biomass
    '''
    Ag = min(1, 10*(fungi+myc)/SOM)
    return np.real(Ag)

def calcPVD(PVstruct, pv, Ag, ratioPVBeng, fPVB, tPVB, PVBmax, d, b):
    #all volumes in l/m3
    Beng = b[6] #engineer biomass
    PVB = min(PVBmax*d, ratioPVBeng*Beng*d)
    PVB=max(PVB, (pv[4]-PVstruct[4])*(1-tPVB))
    # absolute volume of burrows ifo engineer biomass (to max)
    # burrows can increase total porosity or push material
    # so moving pore space, the ration between the 2 = parameter fPVB
    pv = np.zeros(5)
    pv[4] = PVstruct[4] + PVB
    pv[1] = PVstruct[1] + Ag*PVstruct[3]/4 + Ag*b[10]*d/1000
    pv[2] = PVstruct[2] + Ag*PVstruct[3]/4 + Ag*b[10]*d/1000
    pv[3] = PVstruct[3] - Ag*PVstruct[3]/2 - (1-fPVB)*PVB
    pv[0] = PVstruct[0]
    
    return pv

    #check this
def calcmodt(temp, topt, tmin, tmax): #calculate temperature modifier
    if temp < tmin:
        value = 0
    elif temp < topt:
        value = 2**((temp-topt)/10)
    elif temp < tmax:
        value = 1
    else:
        value = 0
    return np.real(value)

# returns pressure head at a given volumetric water content according to the van genuchten model
def pressure_head(theta, R, S, alpha, n, m): 
    theta = min(theta,S) # saturated water conent is the maximum 
    return - pow( pow( (S - R)/(theta - R), (1./m))-1., 1./n) / alpha

# returns the volumetric water content at a given pressure head  according to the van genuchten model (Eqn 21)
def water_content(h, R, S, alpha, n, m):
    return R + (S-R)*pow(1. + pow(alpha*abs(h),n),-m)

# returns the effective saturation according to the van genuchten model (dimensionless water content, Eqn 2)
def effective_saturation(h, R, S, alpha, n, m):
    h = min(h,0) # pressure head is negative, zero the maximum
    theta = water_content(h, R, S, alpha, n, m)
    se = (theta-R)/(S-R)
    return se

# returns the hydraulic conductivity according to the van genuchten model (Eqn 8)
def hydraulic_conductivity(h, R, S, alpha, n, m, Ksat):
    se = effective_saturation(h, R, S, alpha, n, m) 
    K = Ksat*math.sqrt(se)*( (1. - pow(1. - pow(se, 1/m),m)) ** 2 )
    return K

def calcPW(PV, precip, PW, drainmax, d, R, S, alpha, n, m, Ksat): #PV=pore volume in l/m3 (absolute value)
    SWtot = sum(PW)    # total soil water, assume in absolute values so m3 water per m2 soil in that layer
    h=pressure_head(SWtot, R, S, alpha, n, m)
    ImaxMat = hydraulic_conductivity(h, R, S, alpha, n, m, Ksat) #maximal infriltration in matrix (mm/day), excluding macropores
        
    SAmacro = PV[4]/(1000*d) #surface area of the macropores
    ImaxTot = ImaxMat+PV[4]+drainmax #assume burrow depth = soil layer depth d (maximal infriltrarion/day)
    if precip > min(ImaxTot, sum(PV)-SWtot+drainmax): #more rain than space below including drainage
        runoff = precip-min(ImaxTot, sum(PV)-SWtot+drainmax)
        pnet = precip-runoff
    else:
        pnet = precip # all precipitation infiltrates in soil
        runoff = 0

    if pnet*(1-SAmacro) < min(ImaxMat, sum(PV[0:4])-sum(PW[0:4])): #all rain can go to matrix
        PW[4] = PW[4]+pnet*SAmacro  #macropores take part of rain always
        pnet = pnet*(1-SAmacro)
        #then fill from smallest pores upwards
        for i in range(0, 4):
            pin = (PV[i]-PW[i])
            if pnet >= pin:
                PW[i] = PW[i]+pin
                pnet = pnet-pin
            else:
                PW[i] = PW[i]+pnet
                pnet = 0
        #but rain that fell on macropores still does drain, and mespores also drain
        if PW[3]+PW[4] < drainmax:
            drain = PW[3]+PW[4]
            PW[3] = 0
            PW[4] = 0
        elif PW[4] < drainmax:
            PW[3] = PW[3]-(drainmax-PW[4])
            PW[4] = 0
            drain = drainmax
        else:
            PW[4] = PW[4]-drainmax            
            drain=drainmax
            
    else: # more rain than space in matrix
        #  do the same for the rain that can get into the matrix
        pmat = ImaxMat
        pnet = pnet-pmat
        #then fill from smallest pores upwards
        for i in range(0, 4):
            pin = (PV[i]-PW[i])
            if pmat >= pin:
                PW[i] = PW[i]+pin
                pmat = pmat-pin
            else:
                PW[i] = PW[i]+pmat
                pmat = 0
        #all other goes to macropores
        PW[4] = PW[4]+pnet

        # macropores and mespores also drain
        if PW[3]+PW[4] < drainmax:
            drain = PW[3]+PW[4]
            PW[3] = 0
            PW[4] = 0
        elif PW[4] < drainmax: #macropores drain, mesopores partially
            PW[3] = PW[3]-(drainmax-PW[4])
            PW[4] = 0
            drain = drainmax
        else:  #saturated soil, macropores and mesopores fill
            if (PV[3]-PW[3]) > PW[4]-drainmax:  #macropores drain to below and mesopores
                PW[3] = PW[3]+PW[4]-drainmax #assume in saturate soil water goes within 1 day to mesopores
                PW[4] = 0
                drain = drainmax
            else:
                PW[4] = PW[4]-(PV[3]-PW[3]) #assume in saturate soil water goes within 1 day to mesopores
                PW[3] = PV[3] #full
                PW[4] = PW[4]-drainmax
                drain = drainmax
    return PW, drain, runoff


    #check here.
def calcgmaxmod(CNbiomass, CNsource, pCN, rec, prec, pH, id):
# effect of CN when N not limiting
#effect of CN when limiting (only for bacteria!) is in the main code
# id is te identity of the organisms (1=bacteria, 2=fungi)
      #effect of recalcitrance
   mRec = 1-prec*rec/100
      #effect of pH
   if id < 2: #for bacteria
       if pH < 3:
           mpH = min(1, 1/((3-pH)*10))
       else:
           mpH = 1
   elif pH > 8: #for fungi
       mpH = min(1, 1/((pH-8)*10))
   else:
       mpH = 1
   
   mCN = min(1, (CNbiomass/CNsource)**pCN) #effect of CN
   value = min(1,mCN*mpH*mRec) #assuming complete additivity
   
   return np.real(value) #changing value/mCN we can enable or not effects of pH and recalcitrance

def calcgmaxEng(GM,pH): #calculate engineer gmax ifo pH
    if pH<3:
        gmaxEng=0
    elif 3<=pH<5:
        gmaxEng=GM/2*(pH-3)
    else:    
        gmaxEng=GM
    return np.real(gmaxEng)
    
def calcgrowth(biomass, source, avail, gmaxmod, Ks): #Monod kinetic equation of growth
    value = min(gmaxmod*avail*source*biomass/(Ks+source), (avail*source/2))
    return np.real(value)
    
def calcresp(temp, T1, R1, Q10): #calculate respiration ifo temperature
    value = R1*(Q10**((temp-T1)/10))
    return value

def calcFaec(gmax, faec, pfaec, CNsource, CNbiomass, R): #calculate faeces
    faecm = faec+pfaec*(CNsource-CNbiomass)/CNsource*faec # ffaecEff in paper
    faecshort = CNsource/CNbiomass+R/gmax-1 # ffaecCN in paper
    value = max(faecm, faecshort)
    return value

def exudation(): #from roots
    value = 0.1
    return np.real(value)

def calcBioturb(worms, biotRate, OM): #calculate bioturbation
    soilDown=worms*biotRate/100*OM
    return soilDown

def calcLittermove(worms, moveRate, Lit): #calculate litter move by engineers
    litterDown=worms*moveRate/100*Lit
    return litterDown


def PET(temp, nd, sunh, hi, alfa):  # Potential Evapotranspiration, Thornthwaite equation (1948)
    sun = sunh/nd # average day length (hours)
    pet = (16*(sun/12)*(nd/30)*(((10*max(0,temp))/hi)**alfa))/nd # in l/m2 at daily scale
    return pet
    
def wl(PW, et): # water lost by evapotranspiration
    if et<PW[4]:
        PW[4]=PW[4]-et
    elif et<(PW[3]+PW[4]):
        PW[3]=PW[3]-(et-PW[4])
        PW[4]=0
    elif et<(PW[2]+PW[3]+PW[4]):
        PW[2]=PW[2]-(et-PW[3]-PW[4])
        PW[4]=0
        PW[3]=0
    elif et<(PW[1]+PW[2]+PW[3]+PW[4]):
        PW[1]=PW[1]-(et-PW[2]-PW[3]-PW[4])
        PW[4]=0
        PW[3]=0
        PW[2]=0
    else:
        PW[4]=0
        PW[3]=0
        PW[2]=0
        PW[1]=0

    return PW


    #check recalcitrance, pH and temperature adaptation
    #hypothesis: how well can fungi & microbes adapt to plants?
    #hypothesis 2: does a disparity in adaptation speed cause long/midterm issues.

def OptimiseParamater(argNum, orgGroup, i, f, B, avail, modt, litterCN, SOMCN, GMAX): #Takes an initial value and tries to optimise for that value, selecting on B to run on animal group orgGroup
    argTable = [avail, modt, GMAX] #Sets the argument arrays into one nested array for code compaction
    currentArgVal = argTable[argNum] #Selects the correct argument array based on the value of argNum

    organismInterval = [[0,10]] #Calculates the stepSize based off of the organism interval.
    lowerBound = organismInterval[orgGroup][0]
    upperBound = organismInterval[orgGroup][1]
    stepSize = upperBound*0.005
    
    #print(currentArgVal[orgGroup])
    ArgValUp = currentArgVal.copy() #initialise tables for increase and decrease
    ArgValDown = currentArgVal.copy()
    ArgValUp[orgGroup] = min(upperBound,currentArgVal[orgGroup]*(1+stepSize)) #Change the paramatervalue for that organism group only
    ArgValDown[orgGroup] = max(lowerBound, currentArgVal[orgGroup]*(1-stepSize))
    
    #Here code to select the argument that will be modified.
    #if selected arg == true then
    #selected Arg = currentArgVal

    
    #Run odeint 3 times, one with the current value, one with the increased and one with the decreased value.
    if argNum == 2: #changes maxgrowth
        Bday = odeint(f, B, [i, i+1], args=(avail, modt, currentArgVal, litterCN, SOMCN))
        Bcurrent = Bday[1, :]
                          
        Bdayup = odeint(f, B, [i, i+1], args=(avail, modt, ArgValUp, litterCN, SOMCN))
        Bup = Bdayup[1, :]
                     
        Bdaydown = odeint(f, B, [i, i+1], args=(avail, modt, ArgValDown, litterCN, SOMCN))
        Bdown = Bdaydown[1, :]
    elif argNum ==1: #changes temperature. Change this to change the three temperature parameters 
        Bday = odeint(f, B, [i, i+1], args=(avail, currentArgVal, GMAX, litterCN, SOMCN))
        Bcurrent = Bday[1, :]
        print(Bcurrent)                  
        
        Bdayup = odeint(f, B, [i, i+1], args=(avail, ArgValUp, GMAX, litterCN, SOMCN))
        Bup = Bdayup[1, :]
        print(Bup)             
        
        Bdaydown = odeint(f, B, [i, i+1], args=(avail, ArgValDown, GMAX, litterCN, SOMCN))
        Bdown = Bdaydown[1, :]
        
    
    else: #changes avail
        Bday = odeint(f, B, [i, i+1], args=(currentArgVal, modt, GMAX, litterCN, SOMCN))
        Bcurrent = Bday[1, :]
                          
        Bdayup = odeint(f, B, [i, i+1], args=(ArgValUp, modt, GMAX, litterCN, SOMCN))
        Bup = Bdayup[1, :]
                     
        Bdaydown = odeint(f, B, [i, i+1], args=(ArgValDown, modt, GMAX, litterCN, SOMCN))
        Bdown = Bdaydown[1, :]

    #debug

    #print(Bday)
    #print(Bdayup)
    #print(Bdaydown)
    
    bestResult = max(Bcurrent[orgGroup], Bup[orgGroup], Bdown[orgGroup])
    #Select the highest value of B and return the argument value associated with it
    if bestResult == Bcurrent[orgGroup]: return currentArgVal[orgGroup]
    elif bestResult == Bup[orgGroup]: return ArgValUp[orgGroup]
    else: return ArgValDown[orgGroup]
    
#breeder maken van temperatuur plateau.
#Als temperatuur ver boven optimum valt, traag beginnen adapteren?
#als groep dood is: geen adapttie
    