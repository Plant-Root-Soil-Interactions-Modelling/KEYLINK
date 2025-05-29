"""
FUNCTIONS FOR KEYLINK MODEL

Created on 01.06.2017 last write 18/11/2021

@author: A Schnepf - G Deckmyn - G Cauwenberg - O Flores
"""

import numpy as np
import math


def calcKD(pH, fClay):  # Cempirically from orchideeSOM model Cammino-serrano et al 2018
    Kd = 0.001226 - 0.000212 * pH + 0.00374 * fClay
    return Kd


def calcPoreSurfaceArea(PV, PRadius, PSA):
    for i in range(4):
        PSA[i] = 2 * PV[i] / PRadius[i]

    return PSA


def calcAvailPot(PV, PW):

    mwater = np.zeros(5)
    if sum(PW) / sum(PV) < 0.5:
        mwatertot = 4 * sum(PW) / sum(PV) * (1 - sum(PW) / sum(PV))
    else:
        mwatertot = 1
    for i in range(4):
        if PW[i] / PV[i] == 1:  # pore size that is filled
            mwater[i] = 0
        elif PW[i] / PV[i] > 0:
            mwater[i] = mwatertot * PV[i] / (PV[i] + PW[i])
            if i == 3:
                mwater[i + 1] = mwatertot * PW[i] / (PV[i] + PW[i])
    availSOMbact = np.sum(
        [1, 1, 1, 1] * PV[1:] * mwater[1:] / sum(PV[1:])
    )  # SOM availability to bact
    availSOMfungi = np.sum(
        [0, 1, 1, 1] * PV[1:] * mwater[1:] / sum(PV[1:])
    )  # SOM availability to fung

    #    MAOMavail = MAOM - MAOMunavail
    availability = np.array([availSOMbact, availSOMfungi])
    return availability


def calcAvail(PV, PW, SOMini, PSA):
    # TODO add MAOM saturation to avaialbility
    # units in water and porosity volumes: l/m3
    mwater = np.zeros(5)
    if sum(PW) / sum(PV) < 0.5:
        mwatertot = 4 * sum(PW) / sum(PV) * (1 - sum(PW) / sum(PV))
    else:
        mwatertot = 1
    for i in range(4):
        if PW[i] / PV[i] == 1:  # pore size that is filled
            mwater[i] = 0
        elif PW[i] / PV[i] > 0:
            mwater[i] = mwatertot * PV[i] / (PV[i] + PW[i])
            if i == 3:
                mwater[i + 1] = mwatertot * PW[i] / (PV[i] + PW[i])

    availSOMbact = np.sum(
        [1, 1, 1, 1] * PV[1:] * mwater[1:] / sum(PV[1:])
    )  # SOM availability to bact
    availSOMfungi = np.sum(
        [0, 1, 1, 1] * PV[1:] * mwater[1:] / sum(PV[1:])
    )  # SOM availability to fung
    availSOMeng = np.sum(
        [1, 1, 1, 1] * PV[1:] * mwater[1:] / sum(PV[1:])
    )  # SOM availability to engineers
    availSOMsap = np.sum(
        [0, 0, 1, 1] * PV[1:] * mwater[1:] / sum(PV[1:])
    )  # SOM availability to detritivores
    availbbvores = np.sum(
        [0, 1, 1, 1] * PV[1:] * mwater[1:] / sum(PV[1:])
    )  # bact availability to bacterivores
    availffvores = np.sum(
        [1, 1, 1] * PV[2:] * mwater[2:] / sum(PV[2:])
    )  # fung availability to fungivores
    availfvorespred = np.sum(
        [0, 1, 1] * PV[2:] * mwater[2:] / sum(PV[2:])
    )  # fungivores availability to predators
    availbvorespred = np.sum(
        [0, 1, 1] * PV[2:] * mwater[2:] / sum(PV[2:])
    )  # bacterivores  availability to predators
    availhvorespred = np.sum(
        [0, 1, 1] * PV[2:] * mwater[2:] / sum(PV[2:])
    )  # herbivores  availability to predators
    availsappred = np.sum(
        [1, 1] * PV[3:] * mwater[3:] / sum(PV[3:])
    )  # sap  availability to predators (sap only in larger pores)
    availengpred = 1  # engineers availability to predators (earthworms can't hide)
    #    SOMunavail      = (PV[0]/sum(PV))*SOMini * MAOMsaturation #SOM in inaccesible pores
    MAOMunavail = (PSA[0] / sum(PSA)) * SOMini
    #    MAOMavail = MAOM-MAOMunavail
    availability = np.array(
        [
            availSOMbact,
            availSOMfungi,
            availSOMeng,
            availSOMsap,
            availbbvores,
            availffvores,
            availfvorespred,
            availbvorespred,
            availhvorespred,
            availsappred,
            availengpred,
            MAOMunavail,
        ]
    )
    return availability


def inputLitter(inLit, CNlit):  # input from plant / litter to soil surface
    value = [inLit, CNlit]
    return value


def inputCtoMyc(inC):  # input from plant / litter to mycorrhizal fungi
    value = inC
    return np.real(value)


def mycNtoPlant(NmyctoPl, inlit, CNlit, inC, CNsoil):  # N through the myc to the plant
    value = min(
        NmyctoPl * inlit / CNlit, 5 * inC / CNsoil
    )  # % of plant requirement is traded by myc
    return np.real(value)


def rootgrowth(rrg):
    value = rrg
    return np.real(value)


def rootTurnover(Broot, TurnoverRate):
    deadRoot = Broot * TurnoverRate
    return deadRoot


def plantNuptake(CNlit, inlit, NmyctoPl):
    plantNup = (1 - NmyctoPl) * inlit / CNlit
    return plantNup


def calcAg(fungi, myc, SOM):
    """
    Fraction aggregation of SOM in function of the fungal and EM biomass
    """
    Ag = min(1, 10 * (fungi + myc) / SOM)
    return np.real(Ag)


def calcPVD(PVstruct, pv, Ag, ratioPVBeng, fPVB, tPVB, PVBmax, d, b):
    # all volumes in l/m3
    Beng = b[6]  # engineer biomass
    PVB = min(PVBmax * d, ratioPVBeng * Beng * d)
    PVB = max(PVB, (pv[4] - PVstruct[4]) * (1 - tPVB))
    # absolute volume of burrows ifo engineer biomass (to max)
    # burrows can increase total porosity or push material
    # so moving pore space, the ration between the 2 = parameter fPVB
    pv = np.zeros(5)
    pv[4] = PVstruct[4] + PVB
    pv[1] = PVstruct[1] + Ag * PVstruct[3] / 4 + Ag * b[10] * d / 1000
    pv[2] = PVstruct[2] + Ag * PVstruct[3] / 4 + Ag * b[10] * d / 1000
    pv[3] = PVstruct[3] - Ag * PVstruct[3] / 2 - (1 - fPVB) * PVB
    pv[0] = PVstruct[0]  # this is MAOM, goed to a max

    return pv


def calcmodt(temp, topt, tmin, tmax):  # calculate temperature modifier
    if temp < tmin:
        value = 0
    elif temp < topt:
        value = 2 ** ((temp - topt) / 10)
    elif temp < tmax:
        value = 1
    else:
        value = 0
    return np.real(value)


# returns pressure head at a given volumetric water content according to the van genuchten model
def pressure_head(theta, R, S, alpha, n, m):
    theta = min(theta, S)  # saturated water conent is the maximum
    return -pow(pow((S - R) / (theta - R), (1.0 / m)) - 1.0, 1.0 / n) / alpha


# returns the volumetric water content at a given pressure head  according to the van genuchten model (Eqn 21)
def water_content(h, R, S, alpha, n, m):
    return R + (S - R) * pow(1.0 + pow(alpha * abs(h), n), -m)


# returns the effective saturation according to the van genuchten model (dimensionless water content, Eqn 2)
def effective_saturation(h, R, S, alpha, n, m):
    h = min(h, 0)  # pressure head is negative, zero the maximum
    theta = water_content(h, R, S, alpha, n, m)
    se = (theta - R) / (S - R)
    return se


# returns the hydraulic conductivity according to the van genuchten model (Eqn 8)
def hydraulic_conductivity(h, R, S, alpha, n, m, Ksat):
    se = effective_saturation(h, R, S, alpha, n, m)
    K = Ksat * math.sqrt(se) * ((1.0 - pow(1.0 - pow(se, 1 / m), m)) ** 2)
    return K


def calcPW(
    PV, precip, PW, drainmax, d, R, S, alpha, n, m, Ksat
):  # PV=pore volume in l/m3 (absolute value)
    SWtot = sum(
        PW
    )  # total soil water, assume in absolute values so m3 water per m2 soil in that layer
    h = pressure_head(SWtot, R, S, alpha, n, m)
    ImaxMat = hydraulic_conductivity(
        h, R, S, alpha, n, m, Ksat
    )  # maximal infriltration in matrix (mm/day), excluding macropores

    SAmacro = PV[4] / (1000 * d)  # surface area of the macropores
    ImaxTot = (
        ImaxMat + PV[4] + drainmax
    )  # assume burrow depth = soil layer depth d (maximal infriltrarion/day)
    if precip > min(
        ImaxTot, sum(PV) - SWtot + drainmax
    ):  # more rain than space below including drainage
        runoff = precip - min(ImaxTot, sum(PV) - SWtot + drainmax)
        pnet = precip - runoff
    else:
        pnet = precip  # all precipitation infiltrates in soil
        runoff = 0

    if pnet * (1 - SAmacro) < min(
        ImaxMat, sum(PV[0:4]) - sum(PW[0:4])
    ):  # all rain can go to matrix
        PW[4] = PW[4] + pnet * SAmacro  # macropores take part of rain always
        pnet = pnet * (1 - SAmacro)
        # then fill from smallest pores upwards
        for i in range(0, 4):
            pin = PV[i] - PW[i]
            if pnet >= pin:
                PW[i] = PW[i] + pin
                pnet = pnet - pin
            else:
                PW[i] = PW[i] + pnet
                pnet = 0
        # but rain that fell on macropores still does drain, and mespores also drain
        if PW[3] + PW[4] < drainmax:
            drain = PW[3] + PW[4]
            PW[3] = 0
            PW[4] = 0
        elif PW[4] < drainmax:
            PW[3] = PW[3] - (drainmax - PW[4])
            PW[4] = 0
            drain = drainmax
        else:
            PW[4] = PW[4] - drainmax
            drain = drainmax

    else:  # more rain than space in matrix
        #  do the same for the rain that can get into the matrix
        pmat = ImaxMat
        pnet = pnet - pmat
        # then fill from smallest pores upwards
        for i in range(0, 4):
            pin = PV[i] - PW[i]
            if pmat >= pin:
                PW[i] = PW[i] + pin
                pmat = pmat - pin
            else:
                PW[i] = PW[i] + pmat
                pmat = 0
        # all other goes to macropores
        PW[4] = PW[4] + pnet

        # macropores and mespores also drain
        if PW[3] + PW[4] < drainmax:
            drain = PW[3] + PW[4]
            PW[3] = 0
            PW[4] = 0
        elif PW[4] < drainmax:  # macropores drain, mesopores partially
            PW[3] = PW[3] - (drainmax - PW[4])
            PW[4] = 0
            drain = drainmax
        else:  # saturated soil, macropores and mesopores fill
            if (PV[3] - PW[3]) > PW[
                4
            ] - drainmax:  # macropores drain to below and mesopores
                PW[3] = (
                    PW[3] + PW[4] - drainmax
                )  # assume in saturate soil water goes within 1 day to mesopores
                PW[4] = 0
                drain = drainmax
            else:
                PW[4] = PW[4] - (
                    PV[3] - PW[3]
                )  # assume in saturate soil water goes within 1 day to mesopores
                PW[3] = PV[3]  # full
                PW[4] = PW[4] - drainmax
                drain = drainmax
    return PW, drain, runoff


def calcgmaxmod(CNbiomass, CNsource, pCN, rec, prec, pH, id):
    # effect of CN when N not limiting
    # effect of CN when limiting (only for bacteria!) is in the main code
    # id is te identity of the organisms (1=bacteria, 2=fungi)
    # effect of recalcitrance
    mRec = 1 - prec * rec / 100
    # effect of pH
    if id < 2:  # for bacteria
        if pH < 3:
            mpH = min(1, 1 / ((3 - pH) * 10))
        else:
            mpH = 1
    elif pH > 8:  # for fungi
        mpH = min(1, 1 / ((pH - 8) * 10))
    else:
        mpH = 1
    if(CNsource<=0):
        CNsource=CNsource
    mCN = min(1, (CNbiomass / CNsource) ** pCN)  # effect of CN
    value = min(1, mCN * mpH * mRec)  # assuming complete additivity

    return np.real(
        value
    )  # changing value/mCN we can enable or not effects of pH and recalcitrance


def calcgmaxEng(GM, pH):  # calculate engineer gmax ifo pH
    if pH < 3:
        gmaxEng = 0
    elif 3 <= pH < 5:
        gmaxEng = GM / 2 * (pH - 3)
    else:
        gmaxEng = GM
    return np.real(gmaxEng)


def calcgrowth(biomass, source, avail, gmaxmod, Ks):  # Monod kinetic equation of growth
    value = min(
        gmaxmod * avail * source * biomass / (Ks + source), (avail * source / 2)
    )
    return np.real(value)


def calcresp(temp, T1, R1, Q10):  # calculate respiration ifo temperature
    value = R1 * (Q10 ** ((temp - T1) / 10))
    return value


def calcFaec(gmax, faec, pfaec, CNsource, CNbiomass, R):  # calculate faeces
    faecm = faec + pfaec * (CNsource - CNbiomass) / CNsource * faec  # ffaecEff in paper
    faecshort = CNsource / CNbiomass + R / gmax - 1  # ffaecCN in paper
    value = max(faecm, faecshort)
    return value


def exudation():  # from roots
    value = 0.1
    return np.real(value)


def calcBioturb(worms, biotRate, OM):  # calculate bioturbation
    soilDown = worms * biotRate / 100 * OM
    return soilDown


def calcLittermove(worms, moveRate, Lit):  # calculate litter move by engineers
    litterDown = worms * moveRate / 100 * Lit
    return litterDown


def PET(
    temp, nd, sunh, hi, alfa
):  # Potential Evapotranspiration, Thornthwaite equation (1948)
    sun = sunh / nd  # average day length (hours)
    pet = (
        16 * (sun / 12) * (nd / 30) * (((10 * max(0, temp)) / hi) ** alfa)
    ) / nd  # in l/m2 at daily scale
    return pet


def wl(PW, et):  # water lost by evapotranspiration
    if et < PW[4]:
        PW[4] = PW[4] - et
    elif et < (PW[3] + PW[4]):
        PW[3] = PW[3] - (et - PW[4])
        PW[4] = 0
    elif et < (PW[2] + PW[3] + PW[4]):
        PW[2] = PW[2] - (et - PW[3] - PW[4])
        PW[4] = 0
        PW[3] = 0
    elif et < (PW[1] + PW[2] + PW[3] + PW[4]):
        PW[1] = PW[1] - (et - PW[2] - PW[3] - PW[4])
        PW[4] = 0
        PW[3] = 0
        PW[2] = 0
    else:
        PW[4] = 0
        PW[3] = 0
        PW[2] = 0
        PW[1] = 0

    return PW


def fCompSpecies(
    B,
    t,
    avail,
    modt,
    GMAX,
    litterCN,
    SOMCN,
    mf,
    CN,
    MCN,
    MREC,
    pH,
    recLit,
    FAEC,
    pfaec,
    rRESP,
    KS,
    DEATH,
    CtoMyc,
):
    (
        availSOMbact,
        availSOMfungi,
        availSOMeng,
        availSOMsap,
        availbbvores,
        availffvores,
        availfvorespred,
        availbvorespred,
        availhvorespred,
        availsappred,
        availengpred,
        SOMunavail,
    ) = avail

    # alternative version for '10 groups': 0=bact, 1=fungi, 2=myc, 3=bvores, 4=fvores, 5=sap,
    # 6=eng, 7=hvores, 8=pred, 9=litter, 10=SOM, 11=roots, 12=CO2, 23 = extra myc (less mutualistic), 24 = pathogenic fungi

    # update GMAX for bacteria, fung and myc GMAX is modified for SOM
    # and litter seperately depending on CN (and possibly recalcitrance)
    # for bact if CN source too high they can't grow
    gmaxblit = (
        mf.calcgmaxmod(CN[0], litterCN, MCN[0], recLit, MREC[0], pH, 1) * GMAX[0]
    )  # gmax for bact on litter
    gmaxbSOM = (
        mf.calcgmaxmod(CN[0], SOMCN, MCN[0], 0.0, MREC[0], pH, 1) * GMAX[0]
    )  # gmax for bact on SOM
    gmaxflit = (
        mf.calcgmaxmod(CN[1], litterCN, MCN[1], recLit, MREC[1], pH, 2) * GMAX[1]
    )  # gmax for fung on litter
    gmaxfSOM = (
        mf.calcgmaxmod(CN[1], SOMCN, MCN[1], 0.0, MREC[1], pH, 2) * GMAX[1]
    )  # gmax for fung on SOM
    gmaxEng = min(mf.calcgmaxEng(GMAX[6], pH), GMAX[6])  # gmax for engineers

    # update faeces for  SAP and engineers
    faeclitEng = min(
        1, mf.calcFaec(gmaxEng, FAEC[6], pfaec[6], litterCN, CN[6], rRESP[6])
    )
    faeclitSAP = min(
        1, mf.calcFaec(GMAX[5], FAEC[5], pfaec[5], litterCN, CN[5], rRESP[5])
    )

    # growth equations for each functional group and for variations in C pools
    bact = (
        modt[0]
        * (
            mf.calcgrowth(B[0], B[10] - SOMunavail, availSOMbact, gmaxbSOM, KS[0])
            + mf.calcgrowth(B[0], B[9], availSOMbact, gmaxblit, KS[0])
        )
        - DEATH[0] * B[0]
        - rRESP[0] * B[0]
        - modt[3] * mf.calcgrowth(B[3], B[0], availbbvores, GMAX[3], KS[3])
    )

    fungi = (
        modt[1]
        * (
            mf.calcgrowth(B[1], B[10] - SOMunavail, availSOMfungi, gmaxfSOM, KS[1])
            + mf.calcgrowth(B[1], B[9], availSOMbact, gmaxflit, KS[1])
        )
        - DEATH[1] * B[1]
        - rRESP[1] * B[1]
        - modt[4] * mf.calcgrowth(B[4], B[1], availffvores, GMAX[4], KS[4])
    )

    myc = (
        mf.inputCtoMyc(CtoMyc)
        + modt[2]
        * (
            mf.calcgrowth(B[2], B[9], availSOMbact, gmaxflit, KS[2])
            + mf.calcgrowth(B[2], B[10] - SOMunavail, availSOMfungi, gmaxfSOM, KS[2])
        )
        - DEATH[2] * B[2]
        - rRESP[2] * B[2]
        - modt[4] * mf.calcgrowth(B[4], B[2], availffvores, GMAX[4], KS[4])
    )
    # myc (being a fungi) has the same availability and gmax than fungi
    bvores = (
        modt[3] * mf.calcgrowth(B[3], B[0], availbbvores, GMAX[3], KS[3])
        - modt[8]
        * (1 + FAEC[8])
        * mf.calcgrowth(B[8], B[3], availbvorespred, GMAX[8], KS[8])
        - DEATH[3] * B[3]
        - rRESP[3] * B[3]
    )

    fvores = (
        modt[4]
        * (
            mf.calcgrowth(B[4], B[1], availffvores, GMAX[4], KS[4])
            + mf.calcgrowth(B[4], B[2], availffvores, GMAX[4], KS[4])
        )
        - modt[8]
        * (1 + FAEC[8])
        * mf.calcgrowth(B[8], B[4], availfvorespred, GMAX[8], KS[8])
        - DEATH[4] * B[4]
        - rRESP[4] * B[4]
    )

    sap = (
        modt[5]
        * (
            mf.calcgrowth(B[5], B[9], availSOMbact, GMAX[5], KS[5])
            + mf.calcgrowth(B[5], B[10] - SOMunavail, availSOMsap, GMAX[5], KS[5])
        )
        - modt[8]
        * (1 + FAEC[8])
        * mf.calcgrowth(B[8], B[5], availsappred, GMAX[8], KS[8])
        - DEATH[5] * B[5]
        - rRESP[5] * B[5]
    )

    eng = (
        modt[6]
        * (
            mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6])
            + mf.calcgrowth(B[6], B[10] - SOMunavail, availSOMeng, gmaxEng, KS[6])
        )
        - modt[8]
        * (1 + FAEC[8])
        * mf.calcgrowth(B[8], B[6], availengpred, GMAX[8], KS[8])
        - DEATH[6] * B[6]
        - rRESP[6] * B[6]
    )

    # roots are avaialble because larger than herbivores
    hvores = (
        modt[7] * mf.calcgrowth(B[7], B[11], 1, GMAX[7], KS[7])
        - modt[8]
        * (1 + FAEC[8])
        * mf.calcgrowth(B[8], B[7], availhvorespred, GMAX[8], KS[8])
        - DEATH[7] * B[7]
        - rRESP[7] * B[7]
    )

    pred = (
        modt[8]
        * (
            mf.calcgrowth(B[8], B[3], availbvorespred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[4], availfvorespred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[5], availsappred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[6], availengpred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[7], availhvorespred, GMAX[8], KS[8])
        )
        - DEATH[8] * B[8]
        - rRESP[8] * B[8]
    )

    litter = (
        -modt[0]
        * mf.calcgrowth(B[0], B[9], availSOMbact, gmaxblit, KS[0])  # eaten by bact
        - modt[1]
        * mf.calcgrowth(B[1], B[9], availSOMbact, gmaxflit, KS[1])  # eaten by fungi
        - modt[2]
        * mf.calcgrowth(B[2], B[9], availSOMbact, gmaxflit, KS[2])  # eaten by myc
        - modt[5]
        * (1 + faeclitSAP)
        * mf.calcgrowth(B[5], B[9], availSOMbact, GMAX[5], KS[5])  # eaten by SAP
        - modt[6]
        * (1 + faeclitEng)
        * mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6])  # eaten by engineers
        + DEATH[5] * B[5]
        + DEATH[6] * B[6]
        + DEATH[7] * B[7]
        + DEATH[8] * B[8]
    )

    som = (
        mf.exudation()
        - modt[0]
        * mf.calcgrowth(B[0], B[10] - SOMunavail, availSOMbact, gmaxbSOM, KS[0])
        - modt[1]
        * mf.calcgrowth(B[1], B[10] - SOMunavail, availSOMfungi, gmaxfSOM, KS[1])
        - modt[2]
        * mf.calcgrowth(B[2], B[10] - SOMunavail, availSOMfungi, gmaxfSOM, KS[2])
        - modt[5]
        * mf.calcgrowth(
            B[5], B[10] - SOMunavail, availSOMsap, GMAX[5], KS[5]
        )  # eaten by SAP
        - modt[6]
        * mf.calcgrowth(
            B[6], B[10] - SOMunavail, availSOMeng, gmaxEng, KS[6]
        )  # eaten by engineers
        + modt[8]
        * FAEC[8]
        * (
            mf.calcgrowth(B[8], B[3], availbvorespred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[4], availfvorespred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[5], availsappred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[6], availengpred, GMAX[8], KS[8])
            + mf.calcgrowth(B[8], B[7], availhvorespred, GMAX[8], KS[8])
        )
        + modt[5] * faeclitSAP * mf.calcgrowth(B[5], B[9], availSOMbact, GMAX[5], KS[5])
        + modt[6] * faeclitEng * mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6])
        + modt[7] * FAEC[7] * mf.calcgrowth(B[7], B[11], 1, GMAX[7], KS[7])
        + DEATH[0] * B[0]
        + DEATH[1] * B[1]
        + DEATH[2] * B[2]
        + DEATH[3] * B[3]
        + DEATH[4] * B[4]
    )

    roots = -modt[7] * (1 + FAEC[7]) * mf.calcgrowth(B[7], B[11], 1, GMAX[7], KS[7])

    co2 = (
        rRESP[0] * B[0]
        + rRESP[1] * B[1]
        + rRESP[2] * B[2]
        + rRESP[3] * B[3]  # CO2 emissions from respiration
        + rRESP[4] * B[4]
        + rRESP[5] * B[5]
        + rRESP[6] * B[6]
        + rRESP[7] * B[7]
        + rRESP[8] * B[8]
    )

    bactResp = rRESP[0] * B[0]  # respiration of bacteria
    funResp = rRESP[1] * B[1]  # respiration of fungi
    EMresp = rRESP[2] * B[2]  # respiration of mycorrhizal fungi
    bactGrowthSOM = modt[0] * mf.calcgrowth(
        B[0], B[10] - SOMunavail, availSOMbact, gmaxbSOM, KS[0]
    )  # growth of bact from eaten SOM
    bactGrowthLit = modt[0] * mf.calcgrowth(
        B[0], B[9], availSOMbact, gmaxblit, KS[0]
    )  # growth of bact from eaten litter
    SOMeaten = modt[0] * mf.calcgrowth(
        B[0], B[10] - SOMunavail, availSOMbact, gmaxbSOM, KS[0]
    )  # SOM eaten by bact
    +modt[1] * mf.calcgrowth(
        B[1], B[10] - SOMunavail, availSOMfungi, gmaxfSOM, KS[1]
    )  # eaten by fungi
    +modt[2] * mf.calcgrowth(
        B[2], B[10] - SOMunavail, availSOMfungi, gmaxfSOM, KS[2]
    )  # eaten by myc
    +modt[5] * mf.calcgrowth(
        B[5], B[10] - SOMunavail, availSOMsap, GMAX[5], KS[5]
    )  # eaten by SAP
    +modt[6] * mf.calcgrowth(
        B[6], B[10] - SOMunavail, availSOMeng, gmaxEng, KS[6]
    )  # eaten by engineers
    LITeaten = modt[0] * mf.calcgrowth(
        B[0], B[9], availSOMbact, gmaxblit, KS[0]
    )  # Litter eaten by bact
    +modt[1] * mf.calcgrowth(
        B[1], B[9], availSOMbact, gmaxflit, KS[1]
    )  # eaten by fungi
    +modt[2] * mf.calcgrowth(B[2], B[9], availSOMbact, gmaxflit, KS[2])  # eaten by myc
    +modt[5] * (1 + faeclitSAP) * mf.calcgrowth(
        B[5], B[9], availSOMbact, GMAX[5], KS[5]
    )  # eaten by SAP
    +modt[6] * (1 + faeclitEng) * mf.calcgrowth(
        B[6], B[9], availSOMbact, gmaxEng, KS[6]
    )  # eaten by engineers
    LITeatenEng = (
        modt[6]
        * (1 + faeclitEng)
        * mf.calcgrowth(B[6], B[9], availSOMbact, gmaxEng, KS[6])
    )  # only litter eaten by enginners

    return [
        bact,
        fungi,
        myc,
        bvores,
        fvores,
        sap,
        eng,
        hvores,
        pred,
        litter,
        som,
        roots,
        co2,
        bactResp,
        funResp,
        EMresp,
        bactGrowthSOM,
        bactGrowthLit,
        SOMeaten,
        LITeaten,
        LITeatenEng,
        0,
    ]


def calcRhizosphere(
    treatmentID,
    Priming,
    POM,
    POM_sub,
    CN_POM,
    MAOMs,
    MAOMs_sub,
    MAOMp,
    MAOMp_sub,
    CN_MAOMp,
    CN_MAOMs,
    bact_DOM,
    bact_DOM_sub,
    CN_bact,
    DOM,
    DOM_sub,
    CN_DOM,
    GMAX,
    DEATH,
    pCN,
    pH,
    rRESPbact,
    KS,
    DOM_EC,
    Priming_max,
    kpriming,
    kPOM_MAOM,
    kMAOMs_MAOMp,
    modtBact,
    fSOM
):
    # describes rhizosphere bacterial growth on DOM

    DOM_Nini = DOM / CN_DOM
    
    # print ('line 511 calcRhizosphere', 'bact_DOM=', bact_DOM, 'CN_DOM', CN_DOM)
    # calcgmaxmod(CNbiomass, CNsource, pCN, rec, prec, pH, id)
    # gmaxbPOM = mf.calcgmaxmod(CN_bact, CN_POM, pCN, 0.0, 0, pH, 1) * GMAX #gmax for bact on POM
    gmaxmod = (
        calcgmaxmod(CN_bact, CN_DOM, pCN, 0, 0, pH, 1) * GMAX
    )  # maximum growth for bacteria growing on DOM g/(g day)

    # calculate substrate derived C in bact and DOM
    bact_DOM_sub_abs = (
        bact_DOM * bact_DOM_sub
    )  # absolute substrate derived C in bacteria [gC/m3]
    DOM_sub_abs = DOM * DOM_sub  # absolute substrate derived C in DOM [gC/m3]

    # calculate growth
    # def calcgrowth(biomass, source, avail, gmaxmod, Ks):
    growth = modtBact * calcgrowth(
        bact_DOM, DOM, 1, gmaxmod, KS * bact_DOM
    )  # Monod kinetic equation of growth  # g day net

    mCN = min(1, (CN_bact / CN_DOM) ** pCN)  # effect of CN
    ExtraGrowth = (
        1 - mCN
    ) * growth  # what didn't yet grow in g/day because of N shortage
    # Priming = False

    # print(treatment, "CN_bact: ", CN_bact, "CN_DOM: ", CN_DOM, "pCN: ", pCN)
    # if treatment == 5:
    #     exit()

    if Priming is True and mCN < 1:  # if Priming is allowed and there was a shortage
        # print('priming active')
        # POM, POM_sub, CN_POM, MAOMs, MAOMs_sub, MAOMp, MAOMp_sub, CN_MAOMp, CN_MAOMs, CN_bact, ExtraGrowth, DOM_sub, DOM_EC, Priming_max, kpriming, kPOM_MAOM, kMAOMs_MAOMp
      
        POM, MAOMs, MAOMp, respPriming, respPriming_sub, PrimingGrowth, DOMusedforPriming, SOMprimed_sub = calcPriming(
            POM,
            POM_sub,
            CN_POM,
            MAOMs,
            MAOMs_sub,
            MAOMp,
            MAOMp_sub,
            CN_MAOMp,
            CN_MAOMs,
            CN_bact,
            ExtraGrowth,
            DOM_sub,
            DOM_EC,
            Priming_max,
            kpriming,
            kPOM_MAOM,
            kMAOMs_MAOMp,
            fSOM
        )
    else:
        respPriming = 0
        respPriming_sub = 0
        PrimingGrowth = 0
        ExtraGrowth = 0
        DOMusedforPriming=0
        SOMprimed_sub=0

    BactTurnover = DEATH * bact_DOM  # death of bacteria before adding today's growth
    respDOM = (
        rRESPbact * bact_DOM
    )  # respiration of DOM-feeding bacteria before adding today's growth without priming effect yet
    respDOM_sub_abs = (
        respDOM * bact_DOM_sub
    )  # what part of this respiration is substrate derived
    respDOM_sub = respDOM_sub_abs / respDOM

    bact_DOM += growth + PrimingGrowth - BactTurnover - respDOM
    
    #calculate how much of PrimingGrowth is done using carbon from SOM and how much from DOM
    PrimingGrowth_SOM = fSOM * PrimingGrowth  #fSOM is the fraction 0-1 of growth realized using carbon from primed SOM
    PrimingGrowth_DOM = (1-fSOM) * PrimingGrowth
    
    bact_DOM_sub_abs += (
        (growth + PrimingGrowth_DOM) * DOM_sub
        + PrimingGrowth_SOM * SOMprimed_sub
        - BactTurnover * bact_DOM_sub
        - respDOM * bact_DOM_sub
    )  # add the corresponding part of growth on DOM and SOMprimed as substrate derived C, subtract corresponding part of death and respiration

    # change DOM / what was eaten and what was added from dying bacteria
    DOM += -growth - DOMusedforPriming + BactTurnover
    # substrate-derived amounts
    DOM_sub_abs += (
        -growth * DOM_sub - DOMusedforPriming * DOM_sub + BactTurnover * bact_DOM_sub
    )  # subtract what has been eaten and add corresponding part of substrate derived C from dead bacteria to DOM
    DOM_sub = DOM_sub_abs / DOM  # recalculate relative substrate derived C in DOM
    # if DOM_sub > 1 :
    #     print('line550 calcgrowth DOM_sub=', DOM_sub, 'DOM_sub_abs', DOM_sub_abs, 'DOM=', DOM, 'growth', growth, 'ExtraG', ExtraGrowth, 'BactTurnover', BactTurnover)
    bact_DOM_sub = (
        bact_DOM_sub_abs / bact_DOM
    )  # recalculate relative substrate derived C in bacteria
    DOM_N = (
        DOM_Nini - (growth + DOMusedforPriming) / CN_DOM + BactTurnover / CN_bact
    )  # hopefully it's correct that when burning ExtraGrowth C, some N was lost as well
    # if DOM < 0:
    #     print ('line 520 calcRhizophere DOM=', DOM, 'bactDOM=', bact_DOM, 'modtBact=', modtBact)

    CN_DOM = DOM / DOM_N

    return (
        DOM,
        DOM_sub,
        DOM_N,
        CN_DOM,
        bact_DOM,
        bact_DOM_sub,
        POM,
        MAOMs,
        MAOMp,
        respDOM,
        respDOM_sub,
        respPriming,
        respPriming_sub,
    )


def calcPriming(
    POM,
    POM_sub,
    CN_POM,
    MAOMs,
    MAOMs_sub,
    MAOMp,
    MAOMp_sub,
    CN_MAOMp,
    CN_MAOMs,
    CN_bact,
    ExtraGrowth,
    DOM_sub,
    DOM_EC,
    Priming_max,
    kpriming,
    kPOM_MAOM,
    kMAOMs_MAOMp,
    fSOM
):
    #calculate energy content stored in the DOM that is still available for decay (assimilation) by rhizosphere microbes
    DOM_E = (
        ExtraGrowth * DOM_EC
    )  
    # what can be primed
    SOMprimable = POM + MAOMs + MAOMp  # MAOMp is primed, decision 13/8/2024

    # old version SOMprimed = [gC] how much gC in POM or MAOM can be decayed with energy in DOM (DOM_E)
    #SOMprimed = Priming_max * SOMprimable * (1 - math.exp(-kpriming * DOM_E))

    # old function, overruled by thi new one
    SOMprimed = min(Priming_max * SOMprimable, ExtraGrowth * DOM_EC)
    
    #either SOM or DOM was limiting,if SOM was limiting calculate used DOM
    DOMusedforPriming=SOMprimed/DOM_EC 
    
    # print('SOMprimed', SOMprimed, 'DOMusedforPriming', DOMusedforPriming )
    # print('Priming_max * SOMprimable', Priming_max * SOMprimable)
    # how much of SOMdecayed will be from POM and how much from MAOM? Assume according to difficulty = k and relative pool size
    factor = 1 / (1 + (kPOM_MAOM * POM / (MAOMs + MAOMp)))
    # can go below 0 by ecaying too much of the favourite
    if (1 - factor) * SOMprimed < POM:  # if POMprimed is less than POM available
        POMprimed = (1 - factor) * SOMprimed
        MAOMprimed = factor * SOMprimed
    else:
        POMprimed = (
            POM / SOMprimable
        ) * SOMprimed  # if demand for POM is bigger than available POM pool, do this so that POM does not become negative
        MAOMprimed = ((MAOMp + MAOMs) / SOMprimable) * SOMprimed
    # split the primed MAOM between MAOMp and MAOMs
    MAOMfactor = 1 / (1 + (kMAOMs_MAOMp * MAOMs / MAOMp))

    if (1 - MAOMfactor) * MAOMprimed < MAOMs:
        MAOMsprimed = (1 - MAOMfactor) * MAOMprimed
        MAOMpprimed = MAOMfactor * MAOMprimed
    else:
        MAOMsprimed = MAOMprimed * (MAOMs / (MAOMp + MAOMs))
        MAOMpprimed = MAOMprimed * (MAOMp / (MAOMp + MAOMs))

    # check
    # if SOMprimed != MAOMsprimed + MAOMpprimed + POMprimed:
    #     print('priming ', SOMprimed - (MAOMsprimed + MAOMpprimed + POMprimed))
    # how much will this SOM decay provide N
    # if (SOMprimable < SOMprimed):
    #     print('SOMprimed, SOMprimable', SOMprimed, SOMprimable)
    # print('priming467 SOMprimed',  SOMprimed)
    NavailPOM = POMprimed / CN_POM
    NavailMAOM = MAOMpprimed / CN_MAOMp + MAOMsprimed / CN_MAOMs
    Navail = NavailPOM + NavailMAOM
    # how much bacterial biomass can be grown from this N
    PotentialPrimingGrowth = Navail * CN_bact
    PrimingGrowth = min(PotentialPrimingGrowth, DOMusedforPriming)
    # print("PotentialPrimingGrowth, DOMusedforPriming",PotentialPrimingGrowth, DOMusedforPriming)
    #calculate how much of PrimingGrowth is done using carbon from SOM and how much from DOM
    PrimingGrowth_SOM = fSOM * PrimingGrowth  #fSOM is the fraction 0-1 of growth realized using carbon from primed SOM
    PrimingGrowth_DOM = (1-fSOM) * PrimingGrowth
    
    #calculate respiration
    respPrim = 0

    respPrim = (
        SOMprimed + DOMusedforPriming - PrimingGrowth
    )  # carbon from primed SOM is respired, the C used for biomass of PrimingGrowth is taken from DOM and then the rest was burnt off for mining for nitrogen
    SOMprimed_sub_abs = (
        POMprimed * POM_sub + MAOMsprimed * MAOMs_sub + MAOMpprimed * MAOMp_sub
    )  # absolute substrate derived C in primed part of SOM pools
    #old version
    # respPrim_sub_abs = (
    #     respPrim_SOMprimed_sub_abs + (DOMusedforPriming - PrimingGrowth) * DOM_sub
    # )  # total substrate derived C respired during priming (including C from burning off DOM)
    #use the calculation above for respiration also to calcualte the substrate derived fraction in SOMprimed
    SOMprimed_sub = SOMprimed_sub_abs/SOMprimed
    #to calculate actual respiration from primed SOM, subtract the carbon used for growth from total primed SOM
    #(assume that carbon is used from all SOM pools proportionally to how they were primed)
    respPrim_SOMprimed_sub_abs = SOMprimed_sub_abs - PrimingGrowth_SOM * SOMprimed_sub
    #similarly to include respiration from DOM, subtract from it what was used for priming growth
    respPrim_sub_abs = (
        respPrim_SOMprimed_sub_abs + (DOMusedforPriming - PrimingGrowth_DOM) * DOM_sub
    )  # total substrate derived C respired during priming (including C from burning off DOM)
    
    respPrim_sub = respPrim_sub_abs / respPrim

    
    POM -= POMprimed
    MAOMs -= MAOMsprimed
    # print("calcPriming line 962")
    # print('fraction primed POM, MAOMs, MAOMp', POMprimed/SOMprimed, MAOMsprimed/SOMprimed, MAOMpprimed/SOMprimed)
    # print('fraction_sub POM, MAOMs, MAOMp', POM_sub, MAOMs_sub, MAOMp_sub)
    # print('fraction SOMprimed and DOM burnt off (DOMusedforPriming - PrimingGrowth)', SOMprimed/respPrim, (DOMusedforPriming - PrimingGrowth)/respPrim)
    # print('SOMprimed_sub', respPrim_SOMprimed_sub_abs/SOMprimed)
    # print('DOM respired from priming_sub', (DOMusedforPriming - PrimingGrowth) * DOM_sub)
    # print('respPrim_sub', respPrim_sub)
    

    # if MAOMpprimed > MAOMp:
    #     print(
    #         "MAOMp v Priming před odečtením: ",
    #         MAOMp,
    #         "MAOMpprimed: ",
    #         MAOMpprimed,
    #         "kpriming: ",
    #         kpriming,
    #         "SOMprimed: ",
    #         SOMprimed,
    #         "SOMprimable: ",
    #         SOMprimable
    #     )

    MAOMp -= MAOMpprimed

    #     # print("how much was priming growth compared to priming potential growth and ExtraGrowth", PrimingGrowth, PotentialPrimingGrowth, ExtraGrowth) #let's see if we always realize all
    # # else:

    #     # print("priming should be active but is not", PrimingGrowth, ExtraGrowth)

    return POM, MAOMs, MAOMp, respPrim, respPrim_sub, PrimingGrowth, DOMusedforPriming, SOMprimed_sub


# POM, MAOMs, MAOMp, respPriming, respPriming_sub, PrimingGrowth


def calcMAOM(
    MicrobialC,
    DOM_N,
    CN_DOM,
    fractionSA,
    MAOMp,
    MAOMp_sub,
    maxMAOMp,
    DOM,
    DOM_sub,
    MAOMs,
    MAOMs_sub,
    maxMAOMs,
    MAOMsmaxrate,
    MAOMpmaxrate,
    MM_DOM_MAOM,
    maxEffectBactMAOM,
    MM_Bact_MAOM,
    maxEffectN_MAOM,
    MM_N_MAOM,
    maxEffectSA_MAOM,
    MM_SA_MAOM,
    CN_MAOMp,
    CN_MAOMs,
):
    # MAOM formation towards saturation
    # Flow from disolved (DOM) to MAOM (mineral associated) organic matter
    # depends on available DOM, bacteria, N availability (N in DOM) and the size of the rhizosphere/surface area
    fMic = max(
        0, 1 - maxEffectBactMAOM * MM_Bact_MAOM / (MicrobialC + MM_Bact_MAOM)
    )  # between 0 and 1
    fN = max(
        0, 1 - maxEffectN_MAOM * MM_N_MAOM / (DOM_N + MM_N_MAOM)
    )  # between 0 and 1

    # fraction of the soil layer rooted/hyphenated, effect of surface area included using that of hyphae as max
    # maxSurfaceArea= soil_input.get('layerThickness') * plant_input.get('maxRootDensity') * soilbiota_input.get('HyphalExploration') * 2 * variables_df.get('PlantWaterFraction') / ((1-variables_df.get('PlantWaterFraction')) * 1000 * soilbiota_input.get('HyphalRadius'))
    # if (variables_df.get('AMvolume')[i])>0:
    fRhizosphere = max(
        0, 1 - maxEffectSA_MAOM * MM_SA_MAOM / (fractionSA + MM_SA_MAOM)
    )  # between 0 and 1

    # else:
    #        fRhizosphere = min(1, max(0.0001, 0.01 * variables_df.get('RootSurfaceLayer')[i] / maxSurfaceArea))

    # MAOMp = primary, needs to be formed first from DOM, MAOMs = secondary, depends on MAOMp
    # both MAOMp and MAOMs saturate, calculate fSatMAOMp that is 0 when reaches saturation
    # first calculate the potential current rate of formation as michaelis menten equation so depending on DOM concentration and going to a maximum

    fSatMAOMp = 1 - (MAOMp / maxMAOMp)
    dMAOMp = (
        DOM
        * fMic
        * fN
        * fRhizosphere
        * fSatMAOMp
        * MAOMpmaxrate
        * DOM
        / (DOM + MM_DOM_MAOM)
    )

    # then the same for MAOMs, not influenced by fN, based on the study of Koppitke et al.2020
    fSatMAOMs = 1 - (MAOMs / maxMAOMs)

    dMAOMs = (
        DOM * fMic * fRhizosphere * fSatMAOMs * MAOMsmaxrate * DOM / (DOM + MM_DOM_MAOM)
    )

    # if dMAOMs < 0:
    #     print(
    #         "dMAOMs: ",
    #         dMAOMs,
    #         "fSatMAOMs: ",
    #         fSatMAOMs,
    #         "MAOMs: ",
    #         MAOMs,
    #         "maxMAOMs: ",
    #         maxMAOMs,
    #         "MAOMp: ",
    #         MAOMp,
    #     )
    #     exit()
    # if not yet saturated so there is still some potential rate of MAOM formation
    # MAOMs takes over CN of DOM, so CN of MAOMs changes but that of DOM does not
    # if dMAOMs <= 0 or dMAOMp <= 0: print("calcMAOM line 591", dMAOMs, dMAOMp)
    # preparation for proportions calculations
    DOM_sub_abs = DOM * DOM_sub
    MAOMs_sub_abs = MAOMs * MAOMs_sub
    MAOMp_sub_abs = MAOMp * MAOMp_sub
    # secondary MAOM formation
    # if dMAOMs >0: safety that is not needed anymore
    # print ('line 572 calcMaom DOM,MAOMs, dMAOMs, fSatMAOMs, maxMAOMs, MAOMp =', DOM,MAOMs, dMAOMs, fSatMAOMs, maxMAOMs,MAOMp)
    DOM = DOM - dMAOMs
    DOM_N -= dMAOMs / CN_DOM
    # CN_DOM = DOM/DOM_N # calculate new CN of DOM pool
    CN_MAOMs = (MAOMs + dMAOMs) / (MAOMs / CN_MAOMs + dMAOMs / CN_DOM)
    MAOMs = MAOMs + dMAOMs

    # secondary MAOM formation
    # if dMAOMp > 0: not needed
    # MAOMp is--has high N, with constant CN ratio so changes the CN ration of the DOM
    # but limited by N in DOM
    if dMAOMp / CN_MAOMp>DOM_N: # if there is not enough N   
          dMAOMp = DOM_N*0.9/CN_MAOMp
    DOM_N -= dMAOMp / CN_MAOMp
    MAOMp = MAOMp + dMAOMp
    DOM = DOM - dMAOMp
    CN_DOM = DOM / DOM_N    
    # substrate-derived proportion changes calculations
    # changes in absolute pools
    DOM_sub_abs -= (
        dMAOMs + dMAOMp
    ) * DOM_sub  # subtract was was taken away from DOM, maybe not needed?
    MAOMs_sub_abs += dMAOMs * DOM_sub  # and what was added to MAOMs
    MAOMp_sub_abs += dMAOMp * DOM_sub  # and what was added to MAOMp
    # recalculate proportions
    DOM_sub = DOM_sub_abs / DOM
    MAOMs_sub = MAOMs_sub_abs / MAOMs
    MAOMp_sub = MAOMp_sub_abs / MAOMp
    # if DOM < 0:
    #     print ('line 579 calcMaom DOM=', DOM)
    if CN_DOM<0:
        CN_DOM=CN_DOM
    return DOM, DOM_N, CN_DOM, DOM_sub, MAOMp, MAOMp_sub, MAOMs, MAOMs_sub, CN_MAOMs


# calculates total porosity from organic matter content and bulk density
def calcTotalporosity(SOM, BD):
    BD = BD / 1000  # change units from kg/m3 to g/cm3
    percSOM = SOM / BD / 10000 * 1.72  # SOM change units from [gC/m3] to % SOM
    # print('percSOM', percSOM)
    Ds = 100 / (percSOM / 1.35 + (100 - percSOM) / 2.65)  # calculate particle density
    # print('particle density', Ds)
    TP = 1 - BD / Ds  # in m3/m3
    # print('Total porosity in l/m3', TP * 1000)
    return TP


# calculates van genuchten model parameters of a water retention curve using pedotransfer function of Tian et al. (2021)
def calcVangenuchten(BD, SOM, fClay, fSand):
    BD = BD / 1000  # change units from kg/m3 to g/cm3
    OC = SOM / BD / 10000  # change units from g/m3 to %
    # print('OC', OC)
    clay = fClay * 100  # change units from g/g to %
    sand = fSand * 100  # change units from g/g to %
    # print('clay and sand', clay, sand)
    S = -0.3334 * BD + 0.0005 * clay + 0.8945
    R = 0.0115 * BD * pow(clay, 0.7489)
    alpha = (0.0012 * sand + 0.0001 * clay + 0.0089 * OC + 0.0101) * pow(BD, -2.5325)
    n = (-0.0034 * sand - 0.0186 * clay - 0.0351 * OC + 1.1477) * BD + (
        0.0068 * sand + 0.0217 * clay + 0.0047 * OC + 0.0080
    )
    m = 1 - 1 / n

    # version of pedotransfer functions without considering OC
    # S = -0.3311 * BD + 0.0005 * clay + 0.8916
    # R = 0.0112 * BD * pow(clay, 0.7550)
    # alpha = (0.0014 * sand + 0.0001 * clay + 0.0159) * pow(BD, -2.8834)
    # n = (-0.0046 * sand - 0.0212 * clay + 1.3398) * BD + (0.0079 * sand + 0.0250 * clay - 0.2617)
    # m = 1 - 1/n

    # print ('Van genuchten model parameters', R, S, alpha, n, m)
    return R, S, alpha, n, m


def calcPoresDistribution(
    R, S, alpha, n, m, TP
):  # calculates pore volume in different pore size classes from van genuchten model parameters

    PV = np.zeros(5)  # initialize empty array
    # pore volume for each pore size class [l/m3]
    PV[0] = 1000 * water_content(
        30000, R, S, alpha, n, m
    )  # vol inac pores, l/m3, pF = 4.5 h=30000 cm
    PV[1] = 1000 * (
        water_content(1500, R, S, alpha, n, m) - water_content(30000, R, S, alpha, n, m)
    )  # vol bact pores, l/m3
    PV[2] = 1000 * (
        water_content(100, R, S, alpha, n, m) - water_content(1500, R, S, alpha, n, m)
    )  # vol micro pores, l/m3
    PV[3] = 1000 * (
        water_content(2, R, S, alpha, n, m) - water_content(100, R, S, alpha, n, m)
    )  # vol meso pores, l/m3
    PV[4] = 1000 * (TP - water_content(2, R, S, alpha, n, m))  # vol macro pores, l/m3
    # saturation = water_content(1, R, S, alpha, n, m)
    # print('saturation', saturation)
    # print('PV in l/m3', PV)
    return PV
