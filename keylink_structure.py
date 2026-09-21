"""
KEYLINK soil structure: pore size distribution (textural porosity, aggregation by fungi and
mycorrhiza, biopores of the ecosystem engineers; Flores et al. 2021, Eqs. 27-34), the availability
of every C pool to every consumer given the pore size distribution and the water content of the pore
classes (Eqs. 14-16), and the bioturbation fluxes of the engineers.
"""
import numpy as np

from keylink_names import ENG, SOM, INACC, BACTPORES, MICRO, MESO, MACRO


def calcAg(fungi, myc, SOM):
    '''
    Fraction aggregation of SOM in function of the fungal and EM biomass
    '''
    Ag = min(1, 10*(fungi+myc)/SOM)
    return np.real(Ag)


def calcPVD(PVstruct, pv, Ag, ratioPVBeng, fPVB, tPVB, PVBmax, d, b):
    #all volumes in l/m3
    Beng = b[ENG] #engineer biomass
    PVB = min(PVBmax*d, ratioPVBeng*Beng*d)
    PVB=max(PVB, (pv[MACRO]-PVstruct[MACRO])*(1-tPVB))   # Flores et al. 2021, Eq. 27
    # absolute volume of burrows ifo engineer biomass (to max)
    # burrows can increase total porosity or push material
    # so moving pore space, the ration between the 2 = parameter fPVB
    pv = np.zeros(5)                                                                    # Eqs. 31-34:
    pv[MACRO] = PVstruct[MACRO] + PVB
    pv[BACTPORES] = PVstruct[BACTPORES] + Ag*PVstruct[MESO]/4 + Ag*b[SOM]*d/1000
    pv[MICRO] = PVstruct[MICRO] + Ag*PVstruct[MESO]/4 + Ag*b[SOM]*d/1000
    pv[MESO] = PVstruct[MESO] - Ag*PVstruct[MESO]/2 - (1-fPVB)*PVB
    pv[INACC] = PVstruct[INACC]
    
    return pv


def calcAvail(PV, PW, SOMini):
    #units in water and porosity volumes: l/m3
    mwater = np.zeros(5)
    if sum(PW)/sum(PV) < 0.5:
        mwatertot = 4*sum(PW)/sum(PV)*(1-sum(PW)/sum(PV))
    else:
        mwatertot = 1
    # Flores et al. (2021) Eqs. 15-16: activity is in the pores that are neither water-logged nor dry,
    # i.e. in the highest pore class that holds water (the 'partially filled' class, Eq. 15) and in the
    # class one above it, which is assumed not yet completely dry (Eq. 16). Pores are filled from the
    # smallest class upwards, so this is the largest class k with PW[k] > 0. A class that is exactly
    # full counts as partially filled with water fraction 1: all activity is then in the class above
    # (e.g. at field capacity: micropores full, mesopores active). This makes the modifier continuous
    # when the water level crosses a class boundary.
    # (water in the macropores only defines the level when the mesopores below them are full)
    filled = [i for i in range(5) if PW[i] > 0 and (i < 4 or PW[3] >= PV[3])]
    if filled:
        k = max(filled)
        PvolW = PW[k]                  # water filled pore volume of that class
        PvolA = PV[k]-PW[k]            # aerated pore volume of that class
        mwater[k] = mwatertot*PvolA/(PvolA+PvolW)          # Eq. 15
        if k < 4:
            mwater[k+1] = mwatertot*PvolW/(PvolA+PvolW)    # Eq. 16
                
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


def calcBioturb(worms, biotRate, OM): #calculate bioturbation
    soilDown=worms*biotRate/100*OM
    return soilDown


def calcLittermove(worms, moveRate, Lit): #calculate litter move by engineers
    litterDown=worms*moveRate/100*Lit
    return litterDown
