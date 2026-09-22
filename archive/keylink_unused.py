"""
Functions that are no longer called by the KEYLINK model, kept here for reference only.
NOT imported by the model; do not add imports of this module to the model code.

* calcKD: adsorption coefficient of dissolved organic matter (Flores et al. 2021, Eq. 37,
  after Camino-Serrano et al. 2018). Part of the leaching concept (Eqs. 35-37), which is not
  implemented: SOM is a single pool without a DOM fraction in this model version.
* inputLitter: returned [litter input, C:N of the litter input] as a list; the model uses the
  configuration fields cfg.inputLit and cfg.CNlit directly.
"""

def calcKD(pH,fClay): #Cempirically from orchideeSOM model Cammino-serrano et al 2018
    Kd=0.001226-0.000212*pH+0.00374*fClay
    return Kd
    


def inputLitter(inLit, CNlit):  # input from plant / litter to soil surface
    value = [inLit,CNlit]
    return value
