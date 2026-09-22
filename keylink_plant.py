"""
KEYLINK plant interface: the external inputs from the (not modelled) vegetation - litter fall,
root growth and turnover, root exudation, plant C to mycorrhiza - and the N exchange between soil
and plant (via mycorrhiza and direct uptake). These are the hooks for coupling to a vegetation
model; in the stand-alone version they are constants from KL_runparams.txt.
"""
import numpy as np


def inputCtoMyc(inC):     # input from plant / litter to mycorrhizal fungi
    """Plant C allocated to the mycorrhizal fungi, g C m-3 day-1 (input: data or plant model).

    inC: the daily allocation, g C m-3 day-1.
    """
    value = inC
    return np.real(value)
    


def rootgrowth(rrg):
    """Fine root growth, g C m-3 day-1 (input, from data or a plant model).

    rrg: the daily root growth, g C m-3 day-1.
    """
    value = rrg
    return np.real(value)
    


def rootTurnover(Broot,TurnoverRate):
    """Fine root litter production, g C m-3 day-1.

    Broot: fine root biomass, g C m-3; TurnoverRate: root turnover, day-1.
    """
    deadRoot=Broot*TurnoverRate
    return deadRoot


def exudation(): #from roots
    """Root exudation into the SOM pool, g C m-3 day-1 (constant placeholder for a plant model)."""
    value = 0.1
    return np.real(value)


def mycNtoPlant(NmyctoPl, inlit, CNlit, inC, CNsoil): #N through the myc to the plant
    """N delivered to the plant by the mycorrhizal fungi, g N m-3 day-1.

    NmyctoPl: fraction of the plant N uptake traded by the mycorrhiza, 0-1; inlit: litter
    input used as proxy of the plant N demand, g C m-3 day-1; CNlit: C:N of that litter, -;
    inC: plant C allocated to the mycorrhiza, g C m-3 day-1; CNsoil: C:N of litter+SOM, -.
    Limited by what the mycorrhiza can trade for the C it receives (5*inC/CNsoil).
    """
    value=min(NmyctoPl*inlit/CNlit, 5*inC/CNsoil) #% of plant requirement is traded by myc
    return np.real(value)


def plantNuptake(CNlit, inlit, NmyctoPl):
    """Direct plant N uptake from the mineral N pool, g N m-3 day-1.

    CNlit: C:N of the litter input, -; inlit: litter input as proxy of the plant N demand,
    g C m-3 day-1; NmyctoPl: fraction of the demand delivered by the mycorrhiza, 0-1 (the
    rest is taken up directly).
    """
    plantNup=(1-NmyctoPl)*inlit/CNlit
    return plantNup
    
