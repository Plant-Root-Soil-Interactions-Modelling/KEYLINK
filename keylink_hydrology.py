"""
KEYLINK hydrology: daily water balance of the soil layer, evapotranspiration and the van Genuchten
soil hydraulic functions (Flores et al. 2021, 'Water flow', Eqs. 3-8), plus the DOM adsorption
coefficient (Eq. 37) of the leaching concept.

Units: pore volumes PV and water contents PW in l m-3 of soil (per class: inaccessible, bacterial,
micro, meso, macro pores, see keylink_names), precipitation / infiltration / drainage in mm day-1
(= l m-2 day-1, identical to l m-3 for the 1 m layer of the example set-up).
"""
import math

import numpy as np


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
    '''Daily water balance of the layer (Flores et al. 2021, Eqs. 4-8): infiltration into the matrix and
    the macropores, runoff, drainage; PW (water per pore class, l/m3) is updated in place and returned.
    ImaxMat is the maximal infiltration rate of the soil matrix ('measured or derived', Eq. 5 ff.); here
    the saturated conductivity Ksat (mm/day) is used.'''
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


def calcKD(pH,fClay): #Cempirically from orchideeSOM model Cammino-serrano et al 2018
    Kd=0.001226-0.000212*pH+0.00374*fClay
    return Kd
    
