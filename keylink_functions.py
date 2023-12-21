"""
FUNCTIONS FOR KEYLINK MODEL

Created on 01.06.2017 last write 18/11/2021

@author: A Schnepf - G Deckmyn - G Cauwenberg - O Flores
"""

import numpy as np
import math

def calcKD(pH,fClay): #Cempirically from orchideeSOM model Cammino-serrano et al 2018
    Kd=0.001226-0.000212*pH+0.00374*fClay
    return Kd

def calcPoreSurfaceArea (PV, PRadius, PSA):
    for i in range(4):
        PSA[i] = 2* PV[i]/PRadius[i]
    
    return PSA


def calcAvailPot(PV, PW):
     
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
            if i == 3:
                mwater[i+1] = mwatertot*PW[i]/(PV[i]+PW[i])
    availSOMbact    = np.sum([1, 1, 1, 1]*PV[1:]*mwater[1:]/sum(PV[1:])) #SOM availability to bact
    availSOMfungi   = np.sum([0, 1, 1, 1]*PV[1:]*mwater[1:]/sum(PV[1:])) #SOM availability to fung
         
    
#    MAOMavail = MAOM - MAOMunavail
    availability    = np.array([availSOMbact, availSOMfungi])
    return availability
def calcAvail(PV, PW, SOMini, PSA):
 #TODO add MAOM saturation to avaialbility
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
            if i == 3:
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
#    SOMunavail      = (PV[0]/sum(PV))*SOMini * MAOMsaturation #SOM in inaccesible pores
    MAOMunavail = (PSA[0]/sum(PSA))*SOMini
#    MAOMavail = MAOM-MAOMunavail
    availability    = np.array([availSOMbact, availSOMfungi, availSOMeng,
                                availSOMsap, availbbvores, availffvores,
                                availfvorespred, availbvorespred,
                                availhvorespred, availsappred, availengpred, MAOMunavail])
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
    pv[0] = PVstruct[0]   # this is MAOM, goed to a max 
    
    return pv

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

def calcPriming(POM, CN_POM, MAOMs, CN_MAOM, bact_RS, CN_bact, DOM,CN_DOM, ExtraGrowth, DOM_EC, Pmax, k, kPOM_MAOM):
        #not sure how to use ExtraGrowth below but I have a feeling I should:D
        #how much nitrogen can be released from SOM with the energy in remaining DOM:
        DOM_E = ExtraGrowth/DOM_EC # total energy stored in the remaining DOM pool [J]
        #SOMdecayed = [gC] how much gC in POM or MAOM can be decayed with energy in DOM (DOM_E)
        #DecayCost = how much energy will be spent on SOM decay, definite integral of a decay price function [J]
              
        #numerically solve this definite integral by increasing x (decayed SOM) until the result,y = DecayCost become slightly bigger than energy available
        #this gives me how much SOM I can decay with the energy available in DOM:
        SOMdecayed = 0.001 # starting value
        DecayCost=0 #starting value
        P0 = Pmax*(0 + (1/k)*math.exp(-k*0)) 
        while DecayCost < DOM_E:
            DecayCost = Pmax*(SOMdecayed + (1/k)*math.exp(-k*SOMdecayed)) - P0
            #print(SOMdecayed, DecayCost)
            SOMdecayed += 0.001
        
        #how much of SOMdecayed will be from POM and how much from MAOM
        
        POMdecayed=SOMdecayed * (kPOM_MAOM/(kPOM_MAOM+1))
        MAOMdecayed=SOMdecayed-POMdecayed
        #how much will this SOM decay provide N
        NavailPOM=POMdecayed/CN_POM
        NavailMAOM=MAOMdecayed/CN_MAOM
        Navail=NavailPOM+NavailMAOM
        #how much bacterial biomass can be grown from this N
        PrimingGrowth = Navail*CN_bact
        respPrim=0
        #if there is enough DOM C around to build new biomass thanks to priming
        if PrimingGrowth > 0: #this should always be true, but let's check
            bact_RS += PrimingGrowth #grow new microbes thanks to priming, but where does this C come from? from POM/MAOM?
            respPrim=SOMdecayed-PrimingGrowth #C for new growth is taken from SOM, so only the rest is respired, 
            #??at which point should we let the new RS bacteria biomass respire?
            DOM-=ExtraGrowth  # quick fix not to divide by zero I burnt off all C in DOM to get N?? is that okay?? this respiration unaccounted for yet
            POM-=POMdecayed 
            MAOMs-=MAOMdecayed
            print("how much was grown from potential growth", PrimingGrowth, ExtraGrowth) #let's see if we always realize all 
        else:

            print("priming should be active but is not", PrimingGrowth, ExtraGrowth)
       
         
                   
        # if Nshortage>0:
        #     if NavailPOM>=Nshortage:   #enough Energy to decay all required POM
        #         Cbact_RS=Cbact_RS+ExtraGrowth
        #         DOM=DOM-ExtraGrowth
        #         DOM_N=DOM_N-ExtraGrowth/CN_DOM
        #         if primingIntensity*ExtraGrowth<POM: #should never go below 0 
        #             POM=POM-wantedPOMMAOMdecay
        #             SOM=SOM-wantedPOMMAOMdecay
        #             respPrim=wantedPOMMAOMdecay
        #     elif NavailPOM+NavailMAOM>=Nshortage: # limited by N, decay all DOM, and as much possible POM, and required MAOM
        #         Cbact_RS=Cbact_RS+ExtraGrowth
        #         DOM=DOM-ExtraGrowth # use all avaialble DOM, is an assumption
        #         DOM_N=DOM_N-ExtraGrowth/CN_DOM
        #         POM=POM-potentialPOMdecay
        #         SOM=SOM-potentialPOMdecay
        #         MAOM=MAOM-(wantedPOMMAOMdecay-potentialPOMdecay)
        #         respPrim=wantedPOMMAOMdecay  # assumption: all C from primed pool is respired
        #     else: # not enough N, decay all potential MAOM and POM
        #         Cbact_RS=Cbact_RS+(Navail+NavailPOM+NavailMAOM)*CNbact
        #         DOM=DOM-ExtraGrowth # use all avaialble DOM, is an assumption
        #         DOM_N=DOM_N-ExtraGrowth/CN_DOM
        #         POM=POM-potentialPOMdecay
        #         SOM=SOM-potentialPOMdecay
        #         MAOM=MAOM-potentialMAOMdecay
        #         respPrim=potentialPOMdecay+potentialMAOMdecay  # assumption: all C from primed pool is respired
        # else:
        #     respPrim=0
        #     Cbact_RS=Cbact_RS+ExtraGrowth
        return DOM, POM, MAOMs, bact_RS, respPrim

def calcRhizosphere (POM, CN_POM, MAOM, CN_MAOM, bact_RS, CN_bact, DOM, CN_DOM, GMAX, DEATH, pCN, pH, res, KSbact, DOM_EC, Pmax, k, kPOM_MAOM):  

    # rhizosphere bacterial gorwth on DOM
    DOM_Nini=DOM/CN_DOM
    #gmaxbPOM = mf.calcgmaxmod(CN_bact, CN_POM, pCN, 0.0, 0, pH, 1)*GMAX #gmax for bact on POM
    gmaxmod= calcgmaxmod(CN_bact, CN_DOM, pCN, 0, 0, pH, 1)*GMAX  #growth per unit of bacterial biomass g/(g day)
    growth= calcgrowth(bact_RS, DOM, 1, gmaxmod, KSbact*bact_RS) #Monod kinetic equation of growth  # g day net
    BactTurnover=DEATH*bact_RS
    bact_RS+=growth-BactTurnover-res*bact_RS
    DOM+=-growth+BactTurnover
    # print(bact_RS, growth, BactTurnover)
    DOM_N=DOM_Nini-growth/CN_DOM+BactTurnover/CN_bact
    CN_DOM=DOM/DOM_N
    resp=res*bact_RS  # from eating DOM
    mCN = min(1, (CN_bact/CN_DOM)**pCN) #effect of CN
    ExtraGrowth=(1-mCN)*growth  # what didn't yet grow in g/day because of N shortage
    if mCN<1:  # if there was a shortage
        print ('priming active')
        DOM, POM, MAOM, bact_RS, respPriming = calcPriming(POM, CN_POM, MAOM, CN_MAOM, bact_RS, CN_bact, DOM,CN_DOM, ExtraGrowth, DOM_EC, Pmax, k, kPOM_MAOM)
    #calcPriming(MAOM,CNbact,fCN, DOM,CN_DOM, SOM, CN_SOM, gmaxmodCN, Nmin, Cbact_RS, resp, primingIntensity)
    else:
       respPriming=0 
    return  DOM,DOM_N, bact_RS, POM, MAOM, resp, respPriming    
    
# def calcMAOMsaturation (maxMAOM,MM_DOMtoMAOM, MAOMsaturation, MAOMmaxrate, bactTurnover, SAclaySilt,RSbact, RSsurface, DOM, DOM_N, CN_DOM):
# #    Microporessaturated= PV[0]*MAOMsaturation/sum(PV[:])
#     MAOM=MAOMsaturation*maxMAOM
# #    EmptyMicropores=PV*(1-MAOMsaturation)
#        #random, needs to be surface area clay & silt
#     FractionRS = RSsurface/SAclaySilt  # to find?
#     dMAOM=FractionRS*DOM*MAOMmaxrate*(1-MAOM/maxMAOM)/ (MM_DOMtoMAOM + DOM)
#     MAOM+=dMAOM
#     MAOMsaturation=MAOM/maxMAOM
#     DOM=DOM-dMAOM
#     DOM_N=DOM_N-dMAOM/CN_DOM
#     # gradual decrease in CN dom because bact respire and turnover
#     return MAOMsaturation,DOM, DOM_N 

def calcMAOMformation (bact_RS, DOM_N, fractionSA, MAOMp, maxMAOMp, DOM, MAOMs, maxMAOMs, MAOMsmaxrate, MAOMpmaxrate, MM_DOM_MAOM,maxEffectBactMAOM,MM_BactMAOM, maxEffectN_MAOM,MM_N_MAOM, maxEffectSA_MAOM,MM_SA_MAOM):
    # MAOM formation towards saturation
    # depending on available decaying FOM and DOM, mineral N, bacteria and the size of the rhizosphere/surface area
         # Flow from disolved (DOM) to MAOM (mineral associated) organic matter
         # TODO DOM to MAOMo,chek function! is towards max, not decay but Michaelis-Menten
     # MAOMp = primary, needs to be formed first from DOM, MAOMs = secondary, depends on MAOMp
         fBact=max(0,1-maxEffectBactMAOM*MM_BactMAOM/(bact_RS+MM_BactMAOM))  # between 0 and 1
         fN=max(0, 1-maxEffectN_MAOM*MM_N_MAOM/(MM_N_MAOM+DOM_N)) # between 0 and 1
                             
     #fraction of the soil layer rooted/hyphenated, effect of surface area included using that of hyphae as max
         #maxSurfaceArea= soil_input.get('layerThickness') * plant_input.get('maxRootDensity') * soilbiota_input.get('HyphalExploration') * 2 * variables_df.get('PlantWaterFraction') / ((1-variables_df.get('PlantWaterFraction')) * 1000 * soilbiota_input.get('HyphalRadius'))
 #        if (variables_df.get('AMvolume')[i])>0:
         fRhizosphere=max(0,1-maxEffectSA_MAOM*MM_SA_MAOM/(fractionSA+MM_SA_MAOM)) # between 0 and 1
             
    #     else:
     #        fRhizosphere=min(1,max(0.0001,0.01*variables_df.get('RootSurfaceLayer')[i]/maxSurfaceArea))
     
     # both MAOMp and MAOMs saturate, calculate fSatMAOMp that is 0 when reaches saturation
     #first calculate the potential current rate of formation as michaelis menten equation so depending on DOM concentration and going to a maximum   
     
         fSatMAOMp= 1 - (MAOMp / maxMAOMp)
         Pot_dMAOMp_dt = DOM * fBact * fRhizosphere * fN * fSatMAOMp * MAOMpmaxrate * DOM / (DOM + MM_DOM_MAOM)
         #if (variables_df.get('DOM')[i])>0.1

     #    if i==0:
    #         print(variables_df.get('DOM')[0])

    # then the same for MAOMs, not influenced by fN
         
         fSatMAOMs = 1-(MAOMs/maxMAOMs)
         Pot_dMAOMs_dt = DOM * fSatMAOMs * fBact * fRhizosphere * MAOMsmaxrate * DOM / (DOM + MM_DOM_MAOM)
             
                 
    # if not yet saturated 
         if Pot_dMAOMs_dt >0:
             MAOMs = MAOMs + Pot_dMAOMs_dt
             DOM = DOM-Pot_dMAOMs_dt
         if Pot_dMAOMp_dt >0:
             MAOMp = MAOMp + Pot_dMAOMp_dt
             DOM = DOM - Pot_dMAOMp_dt
         return DOM, DOM_N, MAOMp, MAOMs   