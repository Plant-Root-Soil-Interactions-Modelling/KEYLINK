"""
KEYLINK food web: growth, respiration, turnover and consumption of the nine functional groups and
the resulting changes of litter, SOM, roots and CO2 (Flores et al. 2021, 'C flow', Eqs. 9-13), the
modifiers of the maximal growth rate (temperature, pH, C:N, recalcitrance; Eqs. 17-26), the faeces
fraction and the density dependent turnover. fCompSpecies is the right-hand side integrated daily
by keylink_core.KeylinkModel.
"""
import numpy as np

from keylink_names import (BACT, FUNGI, MYC, BVORES, FVORES, SAP, ENG, HVORES, PRED, LITTER, SOM, ROOTS, CO2,
                           NGROUPS)
from keylink_plant import inputCtoMyc, exudation


def calcmodt(temp, topt, tmin, tmax): #calculate temperature modifier
    """Temperature modifier of the growth rate, 0-1 (Flores et al. 2021, Eq. 17).

    temp, topt, tmin, tmax: actual, optimal, minimum and maximum temperature, degC.
    Q10-shaped below the optimum, a plateau of 1 between optimum and maximum, 0 outside
    [tmin, tmax). The Q10 of growth is 2 here, while respiration uses the Q10 of the
    parameter file (see calcresp).
    """
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


def calcresp(temp, T1, R1, Q10): #calculate respiration ifo temperature
    """Respiration rate at the temperature of the day, g C g C-1 day-1 (Eq. 11).

    temp: temperature of the day, degC; T1: reference temperature of R1 (the optimum
    temperature), degC; R1: respiration rate at T1, g C g C-1 day-1; Q10: -.
    Unlike the growth modifier this curve has no plateau, so above the optimum respiration
    keeps increasing while growth does not.
    """
    value = R1*(Q10**((temp-T1)/10))
    return value


def calcgmaxmod(CNbiomass, CNsource, pCN, rec, prec, pH, id):
# effect of CN when N not limiting
#effect of CN when limiting (only for bacteria!) is in the main code
# id is te identity of the organisms (1=bacteria, 2=fungi)
      #effect of recalcitrance
   """Combined modifier of gmax for a microbial group on one substrate, 0-1 (Eqs. 18-25).

   CNbiomass: C:N of the consumer, -; CNsource: C:N of the substrate, -; pCN: sensitivity
   to the C:N ratio pmC:N, - (Eqs. 21-22); rec: recalcitrance of the substrate, % (0 for
   SOM); prec: sensitivity to recalcitrance pmRec, - (Eqs. 23-24); pH: -; id: 1 for
   bacteria, >= 2 for fungi and mycorrhiza (selects the pH response of Eqs. 18-19).
   The modifiers are multiplied (Eq. 25) and capped at 1.
   """
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
    """Maximal growth rate of the engineers as a function of pH, g C g C-1 day-1 (Eq. 20).

    GM: unmodified gmax, g C g C-1 day-1; pH: -. Zero below pH 3, linear between 3 and 5,
    GM above.
    """
    if pH<3:
        gmaxEng=0
    elif 3<=pH<5:
        gmaxEng=GM/2*(pH-3)
    else:    
        gmaxEng=GM
    return np.real(gmaxEng)
    


def calcgrowth(biomass, source, avail, gmaxmod, Ks): #Monod kinetic equation of growth
    """Growth of a consumer on one substrate, g C m-3 day-1 (Monod kinetics, Eq. 9).

    biomass: consumer biomass B, g C m-3; source: substrate pool S, g C m-3; avail:
    available fraction fa of that substrate, 0-1 (keylink_structure.calcAvail); gmaxmod:
    maximal growth rate including its modifiers, g C g C-1 day-1; Ks: half saturation
    constant, g C m-3. Returns gmaxmod*fa*S/(Ks+S)*B, limited to half of the available
    substrate per day (a safeguard that is not part of Eq. 9).
    """
    value = min(gmaxmod*avail*source*biomass/(Ks+source), (avail*source/2))
    return np.real(value)
    


def calcDeath(biomass, d, beta=1.0, Bref=1.0):
    '''Turnover (death) flux, g C m-3 day-1. With beta == 1 this is the first-order term d*B of
    Flores et al. (2021) Eq. 12. With beta > 1 the turnover is density dependent following
    Georgiou et al. (2017), Dt = d*B*(B/Bref)**(beta-1): at the reference biomass Bref (e.g. the
    observed biomass) the specific turnover equals d, above it the turnover increases and below
    it decreases. This self-limitation damps the boom-bust behaviour of the exponential growth
    equation (suggested as an add-on in Flores et al. 2021, 'Simulating the variability in gmax').'''
    if beta == 1.0:
        return d*biomass
    return d*biomass*(max(biomass, 0.0)/Bref)**(beta-1.0)


def calcFaec(gmax, faec, pfaec, CNsource, CNbiomass, R): #calculate faeces
    """Fraction of the consumed C that becomes faeces, 0-1 (Eq. 26 plus a C:N constraint).

    gmax: maximal growth rate of the consumer, g C g C-1 day-1; faec: basal faeces fraction
    ffaec, 0-1; pfaec: sensitivity of ffaec to the C:N ratio mfaec, -; CNsource: C:N of the
    food, -; CNbiomass: C:N of the consumer, -; R: respiration rate, g C g C-1 day-1.
    The first term is Eq. 26 (more faeces on poorer food); the second (faecshort) is the
    fraction needed to cover the consumer N demand from food of that C:N and is not in the
    paper. The larger of the two is returned; the caller caps the result at 1.
    """
    faecm = faec+pfaec*(CNsource-CNbiomass)/CNsource*faec # ffaecEff in paper
    faecshort = CNsource/CNbiomass+R/gmax-1 # ffaecCN in paper
    value = max(faecm, faecshort)
    return value


def fCompSpecies(B, t, avail, modt, GMAX, litterCN, SOMCN, CN, MCN, MREC, pH, recLit, FAEC, pfaec, rRESP,
                 KS, DEATH, CtoMyc, BETA=None, BREF=None):
    '''Right-hand side of the daily ODE system (dB/dt of the 13 C pools and the 8 flux outputs, g C m-3
    day-1; Flores et al. 2021, Eqs. 9-13 with the modifiers of Eqs. 17-26). B is the state vector (see
    keylink_names), avail the availability vector of keylink_structure.calcAvail, modt the temperature
    modifiers and rRESP the temperature-corrected respiration rates per group; the remaining arguments
    are the parameters of KeylinkConfig.'''
    # BETA, BREF: optional density-dependent turnover parameters per group (see calcDeath); None = first order
    (availSOMbact, availSOMfungi, availSOMeng, availSOMsap, availbbvores,
     availffvores, availfvorespred, availbvorespred, availhvorespred,
     availsappred, availengpred ,SOMunavail) = avail

    # alternative version for '10 groups': 0=bact, 1=fungi, 2=myc, 3=bvores, 4=fvores, 5=sap,
    # 6=eng, 7=hvores, 8=pred, 9=litter, 10=SOM, 11=roots, 12=CO2, 23 = extra myc (less mutualistic), 24 = pathogenic fungi

    # update GMAX for bacteria, fung and myc GMAX is modified for SOM
    # and litter seperately depending on CN (and possibly recalcitrance)
    #for bact if CN source too high they can't grow   
    gmaxblit = calcgmaxmod(CN[BACT], litterCN, MCN[BACT], recLit, MREC[BACT], pH, 1)*GMAX[BACT] #gmax for bact on litter
    gmaxbSOM = calcgmaxmod(CN[BACT], SOMCN, MCN[BACT], 0.0, MREC[BACT], pH, 1)*GMAX[BACT] #gmax for bact on SOM
    gmaxflit = calcgmaxmod(CN[FUNGI], litterCN, MCN[FUNGI], recLit, MREC[FUNGI], pH, 2)* GMAX[FUNGI] #gmax for fung on litter
    gmaxfSOM = calcgmaxmod(CN[FUNGI], SOMCN, MCN[FUNGI], 0.0, MREC[FUNGI], pH, 2)* GMAX[FUNGI] #gmax for fung on SOM
    # mycorrhiza: same pH response as fungi (id=2) but its own gmax, CN and recalcitrance parameters
    gmaxmlit = calcgmaxmod(CN[MYC], litterCN, MCN[MYC], recLit, MREC[MYC], pH, 2)* GMAX[MYC] #gmax for myc on litter
    gmaxmSOM = calcgmaxmod(CN[MYC], SOMCN, MCN[MYC], 0.0, MREC[MYC], pH, 2)* GMAX[MYC] #gmax for myc on SOM
    gmaxEng = min(calcgmaxEng(GMAX[ENG],pH),GMAX[ENG]) #gmax for engineers
    
    #update faeces for  SAP and engineers
    faeclitEng = min(1,calcFaec(gmaxEng, FAEC[ENG], pfaec[ENG], litterCN, CN[ENG], rRESP[ENG]))
    faeclitSAP = min(1,calcFaec(GMAX[SAP], FAEC[SAP], pfaec[SAP], litterCN, CN[SAP], rRESP[SAP]))

    #turnover (death) flux of each group; density dependent if BETA/BREF are given
    if BETA is None:
        dead = [DEATH[k]*B[k] for k in range(NGROUPS)]
    else:
        dead = [calcDeath(B[k], DEATH[k], BETA[k], BREF[k]) for k in range(NGROUPS)]

    #growth equations for each functional group and for variations in C pools
    bact = (modt[BACT]*(calcgrowth(B[BACT], B[SOM]-SOMunavail, availSOMbact, gmaxbSOM, KS[BACT])
            + calcgrowth(B[BACT], B[LITTER], availSOMbact, gmaxblit, KS[BACT]))
            - dead[BACT] - rRESP[BACT]*B[BACT]
            - modt[BVORES]*calcgrowth(B[BVORES], B[BACT], availbbvores, GMAX[BVORES], KS[BVORES]))

    fungi = (modt[FUNGI]*(calcgrowth(B[FUNGI], B[SOM]-SOMunavail, availSOMfungi, gmaxfSOM, KS[FUNGI])
             + calcgrowth(B[FUNGI], B[LITTER], availSOMbact, gmaxflit, KS[FUNGI]))
             - dead[FUNGI] - rRESP[FUNGI]*B[FUNGI]
             - modt[FVORES]*calcgrowth(B[FVORES], B[FUNGI], availffvores, GMAX[FVORES], KS[FVORES]))

    myc = (inputCtoMyc(CtoMyc)
           + modt[MYC]*(calcgrowth(B[MYC], B[LITTER], availSOMbact, gmaxmlit, KS[MYC])
           + calcgrowth(B[MYC], B[SOM]-SOMunavail, availSOMfungi, gmaxmSOM, KS[MYC]))
           - dead[MYC] - rRESP[MYC]*B[MYC]
           - modt[FVORES]*calcgrowth(B[FVORES], B[MYC], availffvores, GMAX[FVORES], KS[FVORES]))
           # myc (being a fungi) has the same availability as fungi, but its own gmax
    bvores = (modt[BVORES]*calcgrowth(B[BVORES], B[BACT], availbbvores, GMAX[BVORES], KS[BVORES])
              - modt[PRED]*(1+FAEC[PRED])*calcgrowth(B[PRED], B[BVORES], availbvorespred, GMAX[PRED], KS[PRED])
              - dead[BVORES] - rRESP[BVORES]*B[BVORES])

    fvores = (modt[FVORES]*(calcgrowth(B[FVORES], B[FUNGI], availffvores, GMAX[FVORES], KS[FVORES])
              + calcgrowth(B[FVORES], B[MYC], availffvores, GMAX[FVORES], KS[FVORES]))
              - modt[PRED]*(1+FAEC[PRED])*calcgrowth(B[PRED], B[FVORES], availfvorespred, GMAX[PRED], KS[PRED])
              - dead[FVORES] - rRESP[FVORES]*B[FVORES])

    sap = (modt[SAP]*(calcgrowth(B[SAP], B[LITTER],availSOMbact, GMAX[SAP], KS[SAP])
           + calcgrowth(B[SAP], B[SOM]-SOMunavail, availSOMsap, GMAX[SAP], KS[SAP]))
           -modt[PRED]*(1+FAEC[PRED])*calcgrowth(B[PRED], B[SAP], availsappred, GMAX[PRED], KS[PRED])
           - dead[SAP] - rRESP[SAP]*B[SAP])

    eng = (modt[ENG]*(calcgrowth(B[ENG], B[LITTER], availSOMbact, gmaxEng, KS[ENG])
           + calcgrowth(B[ENG], B[SOM]-SOMunavail, availSOMeng, gmaxEng, KS[ENG]))
           - modt[PRED]*(1+FAEC[PRED])*calcgrowth(B[PRED], B[ENG], availengpred, GMAX[PRED], KS[PRED])
           - dead[ENG] - rRESP[ENG]*B[ENG])

    #roots are avaialble because larger than herbivores
    hvores = (modt[HVORES]*calcgrowth(B[HVORES], B[ROOTS], 1, GMAX[HVORES], KS[HVORES])
              - modt[PRED]*(1+FAEC[PRED])*calcgrowth(B[PRED], B[HVORES], availhvorespred, GMAX[PRED], KS[PRED])
              - dead[HVORES] - rRESP[HVORES]*B[HVORES])

    pred = (modt[PRED]*(calcgrowth(B[PRED], B[BVORES], availbvorespred, GMAX[PRED], KS[PRED])
            + calcgrowth(B[PRED], B[FVORES], availfvorespred, GMAX[PRED], KS[PRED])
            + calcgrowth(B[PRED], B[SAP], availsappred, GMAX[PRED], KS[PRED])
            + calcgrowth(B[PRED], B[ENG], availengpred, GMAX[PRED], KS[PRED])
            + calcgrowth(B[PRED], B[HVORES], availhvorespred, GMAX[PRED], KS[PRED]))
            - dead[PRED] - rRESP[PRED]*B[PRED])

    litter = (           
              -modt[BACT]*calcgrowth(B[BACT], B[LITTER], availSOMbact, gmaxblit, KS[BACT]) #eaten by bact
              -modt[FUNGI]*calcgrowth(B[FUNGI], B[LITTER], availSOMbact, gmaxflit, KS[FUNGI]) #eaten by fungi
              -modt[MYC]*calcgrowth(B[MYC], B[LITTER], availSOMbact, gmaxmlit, KS[MYC]) #eaten by myc
              -modt[SAP]*(1+faeclitSAP)*calcgrowth(B[SAP], B[LITTER], availSOMbact, GMAX[SAP], KS[SAP]) #eaten by SAP
              -modt[ENG]*(1+faeclitEng)*calcgrowth(B[ENG], B[LITTER], availSOMbact, gmaxEng, KS[ENG]) # eaten by engineers
              + dead[SAP] + dead[ENG]+ dead[HVORES] + dead[PRED])

    som = (exudation() 
           - modt[BACT]*calcgrowth(B[BACT], B[SOM]-SOMunavail, availSOMbact, gmaxbSOM, KS[BACT])
           - modt[FUNGI]*calcgrowth(B[FUNGI], B[SOM]-SOMunavail, availSOMfungi, gmaxfSOM, KS[FUNGI])
           - modt[MYC]*calcgrowth(B[MYC], B[SOM]-SOMunavail, availSOMfungi, gmaxmSOM, KS[MYC])
           - modt[SAP]*calcgrowth(B[SAP], B[SOM]-SOMunavail, availSOMsap, GMAX[SAP], KS[SAP]) #eaten by SAP
           - modt[ENG]*calcgrowth(B[ENG], B[SOM]-SOMunavail, availSOMeng, gmaxEng, KS[ENG]) # eaten by engineers
           + modt[PRED]*FAEC[PRED] * (calcgrowth(B[PRED], B[BVORES], availbvorespred, GMAX[PRED], KS[PRED])
                        + calcgrowth(B[PRED], B[FVORES], availfvorespred, GMAX[PRED], KS[PRED])
                        + calcgrowth(B[PRED], B[SAP], availsappred, GMAX[PRED], KS[PRED])
                        + calcgrowth(B[PRED], B[ENG], availengpred, GMAX[PRED], KS[PRED])
                        + calcgrowth(B[PRED], B[HVORES], availhvorespred, GMAX[PRED], KS[PRED]))
           + modt[SAP]*faeclitSAP*calcgrowth(B[SAP], B[LITTER], availSOMbact, GMAX[SAP], KS[SAP])
           + modt[ENG]*faeclitEng*calcgrowth(B[ENG], B[LITTER], availSOMbact, gmaxEng, KS[ENG])
           + modt[HVORES]*FAEC[HVORES]*calcgrowth(B[HVORES], B[ROOTS], 1, GMAX[HVORES], KS[HVORES])
           + dead[BACT]+dead[FUNGI]+dead[MYC]+dead[BVORES]+dead[FVORES])

    roots = (- modt[HVORES]*(1+FAEC[HVORES])*calcgrowth(B[HVORES], B[ROOTS], 1, GMAX[HVORES], KS[HVORES]))
   
    co2 = (rRESP[BACT]*B[BACT]+rRESP[FUNGI]*B[FUNGI]+rRESP[MYC]*B[MYC]+rRESP[BVORES]*B[BVORES] #CO2 emissions from respiration
           +rRESP[FVORES]*B[FVORES]+rRESP[SAP]*B[SAP]+rRESP[ENG]*B[ENG]+rRESP[HVORES]*B[HVORES]+rRESP[PRED]*B[PRED])
           
    bactResp=rRESP[BACT]*B[BACT] #respiration of bacteria
    funResp=rRESP[FUNGI]*B[FUNGI] #respiration of fungi
    EMresp=rRESP[MYC]*B[MYC] #respiration of mycorrhizal fungi
    bactGrowthSOM=modt[BACT]*calcgrowth(B[BACT], B[SOM]-SOMunavail, availSOMbact, gmaxbSOM, KS[BACT]) #growth of bact from eaten SOM
    bactGrowthLit=modt[BACT]*calcgrowth(B[BACT], B[LITTER], availSOMbact, gmaxblit, KS[BACT]) #growth of bact from eaten litter
    # total SOM and litter consumption: the sums MUST be enclosed in parentheses,
    # otherwise the continuation lines are separate statements and only the bacterial term is kept
    SOMeaten=(modt[BACT]*calcgrowth(B[BACT], B[SOM]-SOMunavail, availSOMbact, gmaxbSOM, KS[BACT]) #SOM eaten by bact
              + modt[FUNGI]*calcgrowth(B[FUNGI], B[SOM]-SOMunavail, availSOMfungi, gmaxfSOM, KS[FUNGI]) #eaten by fungi
              + modt[MYC]*calcgrowth(B[MYC], B[SOM]-SOMunavail, availSOMfungi, gmaxmSOM, KS[MYC]) #eaten by myc
              + modt[SAP]*calcgrowth(B[SAP], B[SOM]-SOMunavail, availSOMsap, GMAX[SAP], KS[SAP]) #eaten by SAP
              + modt[ENG]*calcgrowth(B[ENG], B[SOM]-SOMunavail, availSOMeng, gmaxEng, KS[ENG])) # eaten by engineers
    LITeaten=(modt[BACT]*calcgrowth(B[BACT], B[LITTER], availSOMbact, gmaxblit, KS[BACT]) #Litter eaten by bact
              + modt[FUNGI]*calcgrowth(B[FUNGI], B[LITTER], availSOMbact, gmaxflit, KS[FUNGI]) #eaten by fungi
              + modt[MYC]*calcgrowth(B[MYC], B[LITTER], availSOMbact, gmaxmlit, KS[MYC]) #eaten by myc
              + modt[SAP]*(1+faeclitSAP)*calcgrowth(B[SAP], B[LITTER], availSOMbact, GMAX[SAP], KS[SAP]) #eaten by SAP
              + modt[ENG]*(1+faeclitEng)*calcgrowth(B[ENG], B[LITTER], availSOMbact, gmaxEng, KS[ENG])) # eaten by engineers
    LITeatenEng=modt[ENG]*(1+faeclitEng)*calcgrowth(B[ENG], B[LITTER], availSOMbact, gmaxEng, KS[ENG]) #only litter eaten by enginners
 
    return [bact, fungi, myc, bvores, fvores, sap,
            eng, hvores, pred, litter, som, roots, co2,
            bactResp,funResp,EMresp,bactGrowthSOM,bactGrowthLit, SOMeaten, LITeaten, LITeatenEng,0]   
