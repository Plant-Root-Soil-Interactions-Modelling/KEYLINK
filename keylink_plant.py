"""
KEYLINK plant interface: the external inputs from the (not modelled) vegetation - litter fall,
root growth and turnover, root exudation, plant C to mycorrhiza - and the N exchange between soil
and plant (via mycorrhiza and direct uptake). These are the hooks for coupling to a vegetation
model; in the stand-alone version they are constants from KL_runparams.txt.
"""
import numpy as np


def inputLitter(inLit, CNlit):  # input from plant / litter to soil surface
    value = [inLit,CNlit]
    return value


def inputCtoMyc(inC):     # input from plant / litter to mycorrhizal fungi
    value = inC
    return np.real(value)
    


def rootgrowth(rrg):
    value = rrg
    return np.real(value)
    


def rootTurnover(Broot,TurnoverRate):
    deadRoot=Broot*TurnoverRate
    return deadRoot


def exudation(): #from roots
    value = 0.1
    return np.real(value)


def mycNtoPlant(NmyctoPl, inlit, CNlit, inC, CNsoil): #N through the myc to the plant
    value=min(NmyctoPl*inlit/CNlit, 5*inC/CNsoil) #% of plant requirement is traded by myc
    return np.real(value)


def plantNuptake(CNlit, inlit, NmyctoPl):
    plantNup=(1-NmyctoPl)*inlit/CNlit
    return plantNup
    
