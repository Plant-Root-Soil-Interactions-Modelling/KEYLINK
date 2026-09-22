"""
Compatibility facade: the functions of the former single module keylink_functions are now in
    keylink_hydrology   water balance, evapotranspiration, van Genuchten functions, DOM adsorption
    keylink_structure   pore size distribution, availability, bioturbation
    keylink_foodweb     growth modifiers, growth/respiration/turnover, fCompSpecies (ODE right-hand side)
    keylink_plant       external plant inputs and plant N uptake
This module re-exports all of them so that `import keylink_functions as mf` keeps working.
New code should import from the specific modules.
"""
from keylink_hydrology import calcPW, PET, wl
from keylink_structure import calcAg, calcPVD, calcAvail, calcBioturb, calcLittermove
from keylink_foodweb import (calcmodt, calcresp, calcgmaxmod, calcgmaxEng, calcgrowth, calcDeath, calcFaec,
                             fCompSpecies)
from keylink_plant import inputCtoMyc, rootgrowth, rootTurnover, exudation, mycNtoPlant, plantNuptake

__all__ = ['calcPW', 'PET', 'wl', 'calcAg', 'calcPVD', 'calcAvail', 'calcBioturb', 'calcLittermove', 
           'calcmodt', 'calcresp', 'calcgmaxmod', 'calcgmaxEng', 'calcgrowth', 'calcDeath', 'calcFaec',
           'fCompSpecies', 'inputCtoMyc', 'rootgrowth', 'rootTurnover', 'exudation', 'mycNtoPlant', 
           'plantNuptake']
