"""
FUNCTIONS FOR KEYLINK MODEL

Created on 01.06.2017 last write 19/7/2021

@author: A Schnepf - G Deckmyn - G Cauwenberg - O Flores
"""

import numpy as np
import random
import math
from scipy.integrate import odeint
import time

    #mathematical function to get y value of a normal distribution with given mean and standard deviation.
def calcNormalDistribution(height,mean,deviation,x):  
    answer = height*np.exp(-(((x-mean)**2)/(2*deviation**2)))
    return answer

#fs = function(x,epsilon,delta) dnorm(sinh(delta*asinh(x)-epsilon))*delta*cosh(delta*asinh(x)-epsilon)/sqrt(1+x^2)
#    function(x,epsilon,delta) dnorm(sinh(a*asinh(x)-b))*a*cosh(a*asinh(x)-b)/sqrt(1+x^2)

#def sinh_archsinh_transformation(x,epsilon,delta):
#    return norm.pdf(np.sinh(delta*np.arcsinh(x)-epsilon))*delta*np.cosh(delta*np.arcsinh(x)-epsilon)/np.sqrt(1+np.power(x,2))


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
def calcmodtOld(temp, topt, tmin, tmax): #calculate temperature modifier. old version for comparison purposes.
    if temp < tmin:
        value = 0
    elif temp < topt:
        value = 2**((temp-topt)/10)
    elif temp < tmax:
        value = 1
    else:
        value = 0
    return np.real(value)

def calctempdata(longlist,dailyTemp,maxLength,intervalPercent): #calculate specified metadata using recent temperature and puts them in a list
    reactionTime = 10 #np.log(generationTime) at a later point?
    if len(longlist) > maxLength:               #If the list exceeds the needed data range the first element is popped.
        del longlist[0]        
    longlist.append(dailyTemp)                  #Append today's temperature to the end of the list

    #'global' data of the timeset, currently unused.
    tempRange = max(longlist) - min(longlist)       #calculate the delta temperature
    tempMean = sum(longlist)/len(longlist)          #calculate the average temperature of that 
    lowBound = np.quantile(longlist, ((100-intervalPercent)/200))
    highBound = np.quantile(longlist,((100+intervalPercent)/200))
    
    
    #Checks for differing trends in temperatures compared to last year, only works if at least one year of data + 10 day buffer exists
    if len(longlist) > 375:
        compareCheck = 1 #set to 1 so that medium term and longterm changes can be made
        #short term reaction system, uses the long list but only takes the last few days [-N:] gets the last N elements from a list
        tempRangeNew = max(longlist[-reactionTime:]) - min(longlist[-reactionTime:])       #calculate the delta temperature
        tempMeanNew = sum(longlist[-reactionTime:])/len(longlist[-reactionTime:])          #calculate the average temperature of that 
        lowBoundNew = np.quantile(longlist[-reactionTime:], ((100-intervalPercent)/200))
        highBoundNew = np.quantile(longlist[-reactionTime:],((100+intervalPercent)/200))
        
        #data for same days, last year (ignores leapyears for simplification) 
        tempRangeOld = max(longlist[-(reactionTime+365):-365]) - min(longlist[-(reactionTime+365):-365])
        tempMeanOld = sum(longlist[-(reactionTime+365):-365])/len(longlist[-(reactionTime+365):-365])          #calculate the average temperature of that 
        lowBoundOld = np.quantile(longlist[-(reactionTime+365):-365], ((100-intervalPercent)/200))
        highBoundOld = np.quantile(longlist[-(reactionTime+365):-365],((100+intervalPercent)/200))
    #code comparing short vs long data  
        if tempRangeNew == 0:   #exception for 0 to catch div by 0 errors. 
            tempRangeRatio = 0   
        elif tempRangeOld == 0:
            tempRangeRatio = 10  #if the old range is 0 then that means the ratio opens up a lot, 10 arbitrarily chosen as a big number
        else:
            tempRangeRatio = tempRangeOld/tempRangeNew
        
        tempMeanDiff = tempMeanOld - tempMeanNew
        lowBoundDiff = lowBoundOld - lowBoundNew
        highBoundDiff = highBoundOld - highBoundNew
    else:
        tempRangeRatio, tempMeanDiff, lowBoundDiff, highBoundDiff, compareCheck = 0,0,0,0,0
    
    return longlist, tempRangeRatio, tempMeanDiff, lowBoundDiff, highBoundDiff, compareCheck



def calcmodt(tempList,temp, topt, tmin, tmax, toptl, tminl, tmaxl, absmin, absmax):
    #calculate a bunch of stats to be used for medium and long-termed behaviour
    tempList, rangeRatio, meanDiff, lowBoundDiff, highBoundDiff, compareCheck = calctempdata(tempList, temp, 1000, 80) #temporary values for logging time and interval check.

    #Only run medium and longterm calculations when meaningful data was collected
    #these are currently static numbers but will be based on thermal tolerance ranges
    if compareCheck == 1:
        #Medium term calculations
        if meanDiff >1:
            tmaxl += 0.02
        elif meanDiff <-1:
            tminl -= 0.02
        
        #longterm calculations
        #slowly shift STDEV if there is sufficienent temperature change. Set to 1 currently but can be turned into a variable based on sensitivity?
        #gain is slower than loss, but loss only happens if temperature anomalies are common. 
        #initial temperatures are NOT recorded because that would assume starting situation is "ideal"
        if highBoundDiff >2:
            absmax += 0.01
            #Also should run a check for consistency of temperature anomalies here
            absmin += 0.05
        elif lowBoundDiff <-2:
            absmin -= 0.01
            #same here
            absmax -= 0.05
        

    #part that calculates current Gmax
    if temp<absmin:
        value = 0
    elif temp>absmax:
        value = 0
    else:
        value = calcShortTemp(temp, topt, tmin, tmax)
    
    return value, tempList, topt, tmin, tmax, tmaxl, tminl, toptl, absmin, absmax #also needs to return the other values, not done yet

def calcShortTemp(temp, topt, tmin, tmax): #generates a Gmax value for the current thermal curve
    deviationleft = (topt - tmin)/5     #set first deviation for the lower bound, assuming tmin as 5 sigma.
    deviationright = (tmax - topt)/5    #set the second deviation for the higher bound
    height = 10
    if temp <= topt: #if loop to make the function work as a split normal distribution function.
        value = calcNormalDistribution(height,topt,deviationleft,temp)
    else:
        value = calcNormalDistribution(height,topt,deviationright,temp)
    return value

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
    
def optimiseTemperature():
    return
    
#breeder maken van temperatuur plateau.
#Als temperatuur ver boven optimum valt, traag beginnen adapteren?
#als groep dood is: geen adapttie
    
def fCompSpecies(B, t, avail, modt, GMAX, litterCN,SOMCN, mf, CN, MCN, MREC, pH, recLit, FAEC, pfaec, rRESP,
         KS, DEATH, CtoMyc):
    (availSOMbact, availSOMfungi, availSOMeng, availSOMsap, availbbvores,
     availffvores, availfvorespred, availbvorespred, availhvorespred,
     availsappred, availengpred ,SOMunavail) = avail

    # alternative version for '10 groups': 0=bact, 1=fungi, 2=myc, 3=bvores, 4=fvores, 5=sap,
    # 6=eng, 7=hvores, 8=pred, 9=litter, 10=SOM, 11=roots, 12=CO2, 23 = extra myc (less mutualistic), 24 = pathogenic fungi

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
            bactResp,funResp,EMresp,bactGrowthSOM,bactGrowthLit, SOMeaten, LITeaten, LITeatenEng,0]
