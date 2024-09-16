# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:08:20 2023

@author: gdeckmyn
"""
import keylink_functions as mf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os


def run_model(AllParam, treatmentVar, mode_):
    # output dataframe list
    outDataframes = []
    # devide 'input' into: parametersToCalibrate, ParametersCalibrated, Inputvariables (run-specific)
    # ********************************************************************************
    # calibrated parameters (= that do NOT change during run but need calibration)

    # the ones we want to calibrate
    # bact_DOM_rel = 0.2 #proportion of bacteria that have access to feeding on DOM (e.g. that are present in rhizophere)
    # DOM_EC = 2 # DOM energetic quality = energy stored per one gram of DOM [J/g] was 5
    # kpriming=0.001 #decay rate of negative exponential decay curve of decay price, was 0., was 0.001
    # KS=5  # C content required to get half the maximal growth of bacteria when decaying DOM [gC/m3]
    # KSfungi= 200 # C content required to get half the maximal growth of fungi when decaying SOM [gC/m3] was 20000
    # KSbact=380 # C content required to get half the maximal growth of bacteria when decaying SOM [gC/m3]
    # kPOM_MAOM = 8 #ratio of POM to MAOM decayed / overall SOM decay is partitioned using this fixed ratios really unavailable, is k-POM/k_MAOM in israel code
    # kMAOMs_MAOMp = 8 #ratio of MAOMs to MAOMp decayed / overall MAOM decay is partitioned using this fixed ratios really unavailable,
    # MAOMpmaxrate = 0.1 #maximum rate of primary MAOM formation
    # MAOMmaxrate=0.2 # max proportion of DOM stabilized in MAOM per day [unitless]
    # MAOMratioSP = 2 #ratio of secondary to primary MAOM
    # maxEffectBactMAOM=0.9 #(half as slow when no bacteria, rate becomes 1-value)
    # maxEffectSA_MAOM = 0.9
    # maxEffectN_MAOM = 0.9
    # MM_N_MAOM =  1
    # MM_Bact_MAOM = 1.10  #bact amount for half speed
    # MM_SA_MAOM = 0.001  # ratio of SA of hyphae/roots to claysiltSA where half max speed of MAOM formation is reached
    # MM_DOM_MAOM=0.025  # DOM concentration for speed being half max speed (Michaelis Menten) [gC/m3]
    # Priming_max=10 #maximum decay price [J/gC] was 10
    bact_DOM_rel = AllParam["bact_DOM_rel"]
    DOM_EC = AllParam["DOM_EC"]
    kpriming = AllParam["kpriming"]
    KS = AllParam["KS"]
    KSfungi = AllParam["KSfungi"]
    KSbact = AllParam["KSbact"]
    kPOM_MAOM = AllParam["kPOM_MAOM"]
    kMAOMs_MAOMp = AllParam["kMAOMs_MAOMp"]
    MAOMpmaxrate = AllParam["MAOMpmaxrate"]
    MAOMsmaxrate = AllParam["MAOMsmaxrate"]
    MAOMratioSP = AllParam["MAOMratioSP"]
    maxEffectBactMAOM = AllParam["maxEffectBactMAOM"]
    maxEffectSA_MAOM = AllParam["maxEffectSA_MAOM"]
    maxEffectN_MAOM = AllParam["maxEffectN_MAOM"]
    MM_N_MAOM = AllParam["MM_N_MAOM"]
    MM_Bact_MAOM = AllParam["MM_Bact_MAOM"]
    MM_SA_MAOM = AllParam["MM_SA_MAOM"]
    MM_DOM_MAOM = AllParam["MM_DOM_MAOM"]
    Priming_max = AllParam["Priming_max"]

    # the ones we use from other calibration
    GMAX = AllParam[
        "GMAX"
    ]  # maximal growth rate for bacteria [gC/(gC day)], KEYLINK was 1.24
    GMAXfungi = AllParam[
        "GMAXfungi"
    ]  # maximal growth rate for fungi [gC/(gC day)], KEYLINK
    mRecBact = AllParam["mRecBact"]  # how sensitive bact are to recalcitrance
    mRecFungi = AllParam["mRecFungi"]  #
    # resp=0.01 #respiration rate for bacteria growing on DOM / ??do we really need a different one? it was set to 0 decided to ditch it and just the next one
    # =0.05  #respiration rate resp, [gC/(gC day)], KEYLINK
    DEATH = AllParam["DEATH"]  # death rate for bacteria [gC/(gC day)], KEYLINK
    DEATHfungi = AllParam["DEATHfungi"]  # death rate for fungi [gC/(gC day)], KEYLINK
    pCN = AllParam[
        "pCN"
    ]  # sensitivity to CN ratio of consumed substrate, values 0-1, taken from KEYLINK (value for bacteria)
    recMAOM = AllParam[
        "recMAOM"
    ]  # recalcitrance of MAOM, (recalcitrance of POM assumed 0)
    RESPbact = AllParam[
        "RESPbact"
    ]  # respiration rate of bacteria, [gC/(gC day)], was 0.05 KEYLINK
    RESPfungi = AllParam[
        "RESPfungi"
    ]  # respiration rate of fungi, [gC/(gC day)], was 0.03 KEYLINK
    T_MAXbact = AllParam["T_MAXbact"]
    T_MINbact = AllParam["T_MINbact"]
    T_OPTbact = AllParam["T_OPTbact"]
    T_MAXfungi = AllParam["T_MAXfungi"]
    T_MINfungi = AllParam["T_MINfungi"]
    T_OPTfungi = AllParam["T_OPTfungi"]
    Q10bact = AllParam["Q10bact"]
    Q10fungi = AllParam["Q10fungi"]

    # input parameters that do not change (=measurable) and are not calibrated, just 'start situation"

    # those which will change for different runs
    DOMinput = treatmentVar[
        "DOMinput"
    ]  # 10 #DOM added in each addition [gC/m3], Jílková2022
    CN_DOMinput = treatmentVar[
        "CN_DOMinput"
    ]  #  #CN of the daily input [unitless], Jílková2022: leachates 80, exudates 6
    CN_MAOMp = treatmentVar["CN_MAOMp"]  # 15 #☼ assumed constant
    CN_POM = treatmentVar["CN_POM"]  # 24 #CN of SOM, Jílková2022
    pH = treatmentVar["pH"]  # 4.1 #Jílková2022
    temp = treatmentVar["temp"]  # 21
    treatment = treatmentVar["treatment"]

    # those different for Jílková 2022 and experiment 2024
    d_freq = 14  # how often is substrate added, every x days, is 14 for Jílková2022, but 21 for experiment 2024
    numDays = 161  # number of days of incubation experiment/how long to run the model, 155 in Jílková2022, 161 in Jílková 2024

    # those that will be the same for all 16 runs
    BD = 800  # bulk density [kg/m³]
    claySA = 800000  # surface area of clay [m²/kg] was 8000000 cm²/g
    CN_bact = 4  # CN of bacteria, from KEYLINK, in Jílková2022 initial CN of microbial biomass is 10
    CN_fungi = 8  # KEYLINK
    fClay = 0.17  # weight fraction [g/g], Jílková2022
    fSilt = 0.24  # weight fraction [g/g], Jílková2022
    maxMAOM = (
        0.86 * (fClay + fSilt) * 100 * BD
    )  # [gC/m3] maximum MAOM, 28208 for Jílková et al. 2022 Georgiou et al. 2022: 86 ± 9 and 48 ± 6 mg C/g silt+clay mineral for HM and LM,
    maxMAOMp = maxMAOM / (MAOMratioSP + 1)  # maximum primary MAOM
    siltSA = 45.4  # m²/kg
    maxSurfaceArea = (
        claySA * BD * fClay + siltSA * BD * fSilt
    )  # total surface area of clay and silt in m²/m³
    # PV = np.array(
    #     [45, 37, 37, 200, 6]
    # )  # pore volume for each pore size class [l/m3]
    # PRadius = np.array(
    #     [0.05, 0.525, 8, 382.5, 875]
    # )  # average radius of each pore size class [µm], defined by KEYLINK
    # PSA = np.zeros(5)
    # PSA = mf.calcPoreSurfaceArea(
    #     PV, PRadius, PSA
    # )  # pore surface area for each pore size class, calculated from  KEYLINK function
    # PW = np.divide(
    #     PV, 2
    # )  # pore water volume, assume all pores half filled , but water is in m³ while volume was in l
    RootHyphaeSurface = 71  # surface area of all roots/hyphae [m2/m3] ??unit correct / look up roots surface area equivalent to that amount of DOM input
    fractionSA = (
        RootHyphaeSurface / maxSurfaceArea
    )  # used in calcMAOM/ fraction of mineral surface area occupied by roots/hyphae

    Priming = True  # flag to enable Priming effect
    Plotting = False  # flag 1 to enable making of plots, so that this can be turned off during sensitivity analysis etc.
    # Bayesian = True  # flag 1 if performing Bayesian
    # Sensitivity = False

    # create dictionary for respiration plot
    respPlot = {
        "soilControl": [],
        "soilLeachates": [],
        "soilExudates": [],
        "subControl": [],
        "subLeachates": [],
        "subExudates": [],
    }

    if mode_ == "Sensitivity" or mode_ == "Bayesian":  # safety
        Plotting = False

    # numruns = 0  # initializing the number of runs

    # run the daily calculations

    # DOMinput_treatments = np.array([10, 10, 0])  # exudates, leachates, control
    # CN_DOMinput_treatments = np.array(
    # [6, 80, 0]

    # treatments = np.array(["exudates", "leachates", "control"])

    # define different output for different modes
    if mode_ == "Sensitivity":
        column_names = [
            "treatment",
            "day",
            "DOMaddition",
            "DOM",
            "bact_DOM",
            "bact",
            "fungi",
            "resp_substrate",
            "resp_soil_baseline",
            "resp_soil",
            "POM",
            "MAOMs",
            "MAOMp",
            "MAOM",
        ]

        results_df = pd.DataFrame(columns=column_names)

    # print("initial", paramsToTestDict)  # check

    if mode_ == "Bayesian":

        # variables for which we have measured data
        # treatment, d, resp, resp_sub, POM, MAOM, bact_total, fungi, POM_sub, MAOM_sub, bact_total_sub, fungi_sub]
        column_names = [
            "treatment",
            "day",
            "resp",
            "resp_sub",
            "POM",
            "MAOM",
            "bact_total",
            "fungi",
            "POM_sub",
            "MAOM_sub",
            "bact_total_sub",
            "fungi_sub",
        ]

        results_df = pd.DataFrame(columns=column_names)

    if mode_ == "Normal":  # normal runs

        column_names = [
            "treatment",
            "day",
            "DOMaddition",
            "DOM",
            "bact_DOM",
            "bact",
            "fungi",
            "respSubstrate",
            "baselineResp",
            "resp",
            "POM",
            "MAOMs",
            "MAOMp",
            "MAOM",
        ]

        results_df = pd.DataFrame(columns=column_names)

    if Plotting:  # only need these if plotting
        time_d = []
        outtreatment = []
        outMAOM = []
        outMAOMp = []
        outMAOMs = []
        outPOM = []
        outDOMadded = []
        outDOM = []
        outDOM_CN = []
        outBact_total = []
        outbact_DOM = []
        outBact = []
        outFungi = []
        outRespSubstrate = []
        outRespSoil = []
        outRespSoilBaseline = []
        outBact_DOM_sub = []
        outBact_total_sub = []
        outBact_sub = []
        outFungi_sub = []
        outDOM_sub = []
        outPOM_sub = []
        outMAOM_sub = []
        outMAOMs_sub = []
        outMAOMp_sub = []
        outResp_sub = []

    # *************************************************************************
    # initializing variables (what changes during run)

    # variables that will be initialized differently for different runs
    bact_total = treatmentVar[
        "bact_total"
    ]  # total biomass of bacteria [gC/m3], was 6 final noadd average from PLFA from Jílková2022
    CN_MAOMs = treatmentVar[
        "CN_MAOMs"
    ]  # estimated but we don't know the true value, assumed to vary with CN_DOM
    fungi = treatmentVar[
        "fungi"
    ]  # biomass of fungi [gC/m3] based on final noadd in Jílková et al. 2022
    MAOM = treatmentVar["MAOM"]  # C in MAOM [gC/m3] average noAdd Jílková2022
    POM = treatmentVar[
        "POM"
    ]  # C in POM [gC/m3], calculated as initialSOM-MAOM using initialSOM from Jílková2022
    TP = None
    PV = None

    if TP == None:
        TP = mf.calcTotalporosity(MAOM + POM, BD)

    if PV == None:
        R, S, alpha, n, m = mf.calcVangenuchten(
            BD, MAOM + POM, fClay, 1 - fClay - fSilt
        )
        PV = mf.calcPoresDistribution(R, S, alpha, n, m, TP)
        PW = np.divide(PV, 2)  # pore water volume, assume all pores half filled  [l/m3]
        # PSA = np.zeros(5)
        # PSA = mf.calcPoreSurfaceArea(
        #     PV, PRadius, PSA
        # )  # pore surface area for each pore size class, calculated from  KEYLINK function, not used currently

    # same for all runs
    availability = np.zeros(3)
    bact = bact_total * (
        1 - bact_DOM_rel
    )  # biomass of bacteria growing on POM and MAOM but not on DOM [gC/m3]
    bact_sub = 0  # proportion of this bacterial carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
    bact_DOM = bact_total * bact_DOM_rel  # biomass of bacteria growing on DOM [gC/m3]
    bact_DOM_sub = 0  # proportion of bacterial carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
    CN_DOM = 0  #
    DOM = 0  # DOM [gC/m3]
    DOM_sub = 0  # relative substrate derived C in DOM /values 0 to 1/, portion of DOM carbon that is substrate derived in contrast to soil-derived / values 0 to 1/ is a ratio between substrate-derived C and total C in DOM
    DOM_N = 0  # set DOM N to zero

    if DOM > 0:
        DOM_N = (
            DOM / CN_DOM
        )  # but if there is some initial DOM, calculate it from CN_DOM
    else:
        DOM_N = 0

    fungi_sub = 0  # proportion of fungal carbon that is substrate derived in contrast to soil-derived / values 0 to 1/

    # MAOMunavail = (
    #     PSA[0] / sum(PSA)
    # ) * MAOM  # the portion of MAOM stored in the smallest pores is really unavailable
    MAOMp = MAOM / (
        MAOMratioSP + 1
    )  # primary MAOM [gC/m3] initialised at the ratio of saturation
    MAOMp_sub = 0  # proportion of MAOMp that is substrate derived in contrast to soil-derived / values 0 to 1/
    MAOMs = MAOM - MAOMp  # secondary MAOM[gC/m3]
    MAOMs_sub = 0  # proportion of MAOMs that is substrate derived in contrast to soil-derived / values 0 to 1/
    POM_sub = 0  # proportion of POM that is substrate derived in contrast to soil-derived / values 0 to 1/
    # check all C
    resp = 0
    DOM_added = 0
    DOM_added_all = 0
    # AllC = DOM + POM + MAOM + bact_total + fungi + resp - DOMadded
    # print(AllC)

    # function coreMAOM
    # def coreMAOM (Bayesian, Sensitivity):

    for d in range(numDays):
        # on day 0 and then every 14 days, add DOM
        if d == 0 or (d % d_freq) == 0:  # on first day and then every d_freq days
            DOM_added = DOMinput  # to keep track of the additions
            DOM_added_all += DOM_added  # keep track of sum of additions
            DOM_sub_abs = DOM * DOM_sub  # absolute substrate derived C in DOM [gC/m3]
            DOM += DOMinput  # add input to the DOM carbon pool
            DOM_sub_abs += DOMinput  # add all input as substrate derived C

            if DOMinput > 0:
                DOM_sub = (
                    DOM_sub_abs / DOM
                )  # update relative substrate derived C in DOM
                # if DOM_sub > 1:
                # print('line265 treatment=', treatment, 'd=', d, 'DOM_sub=', DOM_sub, 'DOM_sub_abs', DOM_sub_abs, 'DOM=', DOM)
                # else: probably not needed
                #     DOM_sub=0
                DOM_N += (
                    DOMinput / CN_DOMinput
                )  # add equivalent amount of N to DON pool
                # print('CN_DOMinput', CN_DOMinput)
                # if DOM_N>0:
                CN_DOM = DOM / DOM_N  # calculate new CN of DOM pool
        else:
            DOM_added = 0
        # saturation of MAOMs depends on amount of MAOMp so recalculated every day
        maxMAOMs = MAOMp * MAOMratioSP  # maximum primary MAOM

        # if maxMAOMs < 0:
        #     print("MAOMp2: ", MAOMp)

        # find t modifier
        modtBact = mf.calcmodt(temp, T_OPTbact, T_MINbact, T_MAXbact)
        modtFungi = mf.calcmodt(temp, T_OPTfungi, T_MINfungi, T_MAXfungi)
        rRESPbact = mf.calcresp(temp, T_OPTbact, RESPbact, Q10bact)
        rRESPfungi = mf.calcresp(temp, T_OPTfungi, RESPfungi, Q10fungi)
        # print(rRESPbact, rRESPfungi)

        # microbial growth on DOM and priming, only susing MAOMs
        if CN_DOM > 0:
            (
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
            ) = mf.calcRhizosphere(
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
            )
        # print('calc.Rhizo')
        #               if (MAOMs<0):
        #                   print('mainLine270 DOM, bact, fungi, MAOMs, MAOMp', DOM,bact, fungi, MAOMs, MAOMp)
        else:

            respDOM = 0
            respDOM_sub = 0
            respPriming = 0
            respPriming_sub = 0

        resp = respDOM + respPriming
        # resp_all += resp
        bact_total = bact_DOM + bact
        MAOM = MAOMs + MAOMp
        # AllC = DOM + POM + MAOM + bact_total + fungi + resp_all - DOM_added

        # MAOM formation
        MicrobialC = (
            bact + bact_DOM + fungi
        )  # all microbes contribute to MAOM formation
        if CN_DOM > 0:
            (
                DOM,
                DOM_N,
                CN_DOM,
                DOM_sub,
                MAOMp,
                MAOMp_sub,
                MAOMs,
                MAOMs_sub,
                CN_MAOMs,
            ) = mf.calcMAOM(
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
            )

        MAOM = MAOMs + MAOMp

        # baseline microbial growth on SOM (without substrate DOM additions)
        availability = mf.calcAvailPot(
            PV, PW
        )  # calculates availability of SOM decomposition by bacteria and fungi, separately, from pore size distribution and soil water
        # calculate maximal growth (gmax) for bacteria/fungi on POM/MAOM separately

        # if CN_MAOMs <= 0:
        #     print("CN_MAOMs: ", CN_MAOMs)

        gmaxbPOM = (
            mf.calcgmaxmod(CN_bact, CN_POM, pCN, 0.0, 0, pH, 1) * GMAX
        )  # gmax for bact on POM
        gmaxfPOM = (
            mf.calcgmaxmod(CN_fungi, CN_POM, pCN, 0.0, 0, pH, 2) * GMAXfungi
        )  # gmax for fungi on POM
        # we assume MAOMp can only be lost through priming, so normal growth uses MAOMs
        gmaxbMAOM = (
            mf.calcgmaxmod(CN_bact, CN_MAOMs, pCN, recMAOM, mRecBact, pH, 1) * GMAX
        )  # gmax for bact on MAOM
        gmaxfMAOM = (
            mf.calcgmaxmod(CN_fungi, CN_MAOMs, pCN, recMAOM, mRecFungi, pH, 2)
            * GMAXfungi
        )  # gmax for fungi on MAOM
        # calculate substrate derived C in bact and fungi
        DOM_sub_abs = DOM * DOM_sub  # recalculate because changesin calc.Rhizosphere
        POM_sub_abs = POM * POM_sub  # recalculate because changes in calc.Rhizosphere
        MAOMs_sub_abs = (
            MAOMs * MAOMs_sub
        )  # recalculate because changes in calc.Rhizosphere and calc.MAOM
        bact_sub_abs = (
            bact * bact_sub
        )  # absolute substrate derived C in bacteria [gC/m3]
        fungi_sub_abs = (
            fungi * fungi_sub
        )  # absolute substrate derived C in fungi [gC/m3]
        #               if (bact<0):
        #                  print('mainLine290 DOM, bact, fungi', DOM,bact, fungi)
        # growth equations (dB/dt) for each functional group and for variations in C and N pools
        # only feed on secondary MAOM
        bactPOMgrowth = modtBact * mf.calcgrowth(
            bact, POM, availability[0], gmaxbPOM, KSbact * bact
        )
        bactMAOMgrowth = modtBact * mf.calcgrowth(
            bact, MAOMs, availability[0], gmaxbMAOM, KSbact * bact
        )
        dbact = bactPOMgrowth + bactMAOMgrowth - DEATH * bact - rRESPbact * bact

        fungiPOMgrowth = modtFungi * mf.calcgrowth(
            fungi, POM, availability[1], gmaxfPOM, KSfungi * fungi
        )
        fungiMAOMgrowth = modtFungi * mf.calcgrowth(
            fungi, MAOMs, availability[1], gmaxfMAOM, KSfungi * fungi
        )
        dfungi = (
            fungiPOMgrowth + fungiMAOMgrowth - DEATHfungi * fungi - rRESPfungi * fungi
        )

        DOM += DEATH * bact + DEATHfungi * fungi  # add dead bacteria and fungi to DOM
        POM += -bactPOMgrowth - fungiPOMgrowth  # subtract what has been eaten from POM
        MAOMs += -bactMAOMgrowth - fungiMAOMgrowth  # and MAOMs

        # update CN DOM
        DOM_N += DEATH * bact / CN_bact + DEATHfungi * fungi / CN_fungi
        CN_DOM = DOM / DOM_N  # recalculate CN DOM

        #    if (-dbact>bact):
        #        print('mainLine307  bact, bactPOMgrowth, POM, bactMAOMgrowth, DEATH*bact, rRESPbact*bact', bact, bactPOMgrowth, POM, bactMAOMgrowth, DEATH*bact, rRESPbact*bact)

        DOM_sub_abs += (
            DEATH * bact * bact_sub + DEATHfungi * fungi * fungi_sub
        )  # add corresponding part of substrate derived C to DOM
        POM_sub_abs -= (bactPOMgrowth + fungiPOMgrowth) * POM_sub
        MAOMs_sub_abs -= (bactMAOMgrowth + fungiMAOMgrowth) * MAOMs_sub
        fungi_sub_abs += (
            fungiMAOMgrowth * MAOMs_sub
            + fungiPOMgrowth * POM_sub
            - DEATHfungi * fungi * fungi_sub
            - rRESPfungi * fungi * fungi_sub
        )  # add the corresponding part of growth on MAOM as substrate derived C, subtract death and respiration
        bact_sub_abs += (
            bactMAOMgrowth * MAOMs_sub
            + bactPOMgrowth * POM_sub
            - DEATH * bact * bact_sub
            - rRESPbact * bact * bact_sub
        )  # add the corresponding part of growth on MAOM as substrate derived C, subtract correspodning part of death and respiration

        baselineRespBact = rRESPbact * bact
        baselineRespBact_sub_abs = (
            baselineRespBact * bact_sub
        )  # what part of this respiration is substrate derived
        baselineRespBact_sub = baselineRespBact_sub_abs / baselineRespBact
        bact += dbact

        baselineRespFungi = rRESPfungi * fungi
        baselineRespFungi_sub_abs = (
            baselineRespFungi * fungi_sub
        )  # what part of this respiration is substrate derived
        baselineRespFungi_sub = baselineRespFungi_sub_abs / baselineRespFungi
        fungi += dfungi

        # update relative substrate derived C proportions
        DOM_sub = DOM_sub_abs / DOM  # relative substrate derived C in DOM
        POM_sub = POM_sub_abs / POM  # relative substrate derived C in DOM
        MAOMs_sub = MAOMs_sub_abs / MAOMs  # relative substrate derived C in DOM
        fungi_sub = (
            fungi_sub_abs / fungi
        )  # update relative substrate derived C in fungi
        # print(' treatment, day, fungi_sub', treatment, d, fungi_sub)
        bact_sub = (
            bact_sub_abs / bact
        )  # update relative substrate derived C in bacteria
        # add up things
        MAOM = MAOMp + MAOMs
        bact_total = bact + bact_DOM
        # baseline respiration without priming
        baselineResp = (
            baselineRespBact + baselineRespFungi + respDOM
        )  # of course this respDOM is higher if previous day DOM-feeding bacteria grew more because of priming
        # all respiration
        resp = baselineResp + respPriming

        # calculate average substrate proportions
        MAOM_sub = MAOMp_sub * (MAOMp / MAOM) + MAOMs_sub * (
            MAOMs / MAOM
        )  # average substrate proportion in MAOM
        bact_total_sub = bact_DOM_sub * (bact_DOM / bact_total) + bact_sub * (
            bact / bact_total
        )  # average substrate proportion in bacteria
        resp_sub = (
            baselineRespBact_sub * (baselineRespBact / resp)
            + baselineRespFungi_sub * (baselineRespFungi / resp)
            + respDOM_sub * (respDOM / resp)
            + respPriming_sub * (respPriming / resp)
        )

        respSubstrate = resp_sub * resp  # substrate derived respiration (absolute)
        respSoil = resp - respSubstrate  # soil-derived respiration (absolute)
        # AllC = DOM + POM + MAOM + bact_total + fungi + resp - DOMadded
        # print('line412', treatment, d, AllC)

        if Plotting:  # save data for Plotting
            outtreatment.append(treatment)
            time_d.append(d)  # store days in an array for plotting
            outDOMadded.append(DOM_added / 0.8)  # change units from gC/m3 µgC/g soil
            outMAOM.append(MAOM / (0.8 * 1000))  # change units from gC/m3 mgC/g soil
            outMAOMp.append(MAOMp / (0.8 * 1000))  # change units from gC/m3 mgC/g soil)
            outMAOMs.append(MAOMs / (0.8 * 1000))  # change units from gC/m3 mgC/g soil)
            outPOM.append(POM / (0.8 * 1000))  # change units from gC/m3 mgC/g soil
            outDOM.append(DOM / 0.8)  # change units from gC/m3 µgC/g soil)
            outBact_total.append(bact_total / 0.8)  # change units from gC/m3 µgC/g soil
            outbact_DOM.append(bact_DOM / 0.8)  # change units from gC/m3 µgC/g soil
            outBact.append(bact / 0.8)  # change units from gC/m3 µgC/g soil
            outFungi.append(fungi / 0.8)  # change units from gC/m3 µgC/g soil
            outRespSubstrate.append(
                respSubstrate / (0.8 * 24)
            )  # change units from gC/m3/day
            outRespSoilBaseline.append(
                baselineResp / (0.8 * 24)
            )  # change units from gC/m3/day
            outRespSoil.append(respSoil / (0.8 * 24))  # change units from gC/m3/day
            # substrate-derived %
            outBact_total_sub.append(bact_total_sub)
            outBact_DOM_sub.append(bact_DOM_sub)
            outBact_sub.append(bact_sub)
            outFungi_sub.append(fungi_sub)
            outDOM_sub.append(DOM_sub)
            outPOM_sub.append(POM_sub)
            outMAOM_sub.append(MAOM_sub)
            outMAOMs_sub.append(MAOMs_sub)
            outMAOMp_sub.append(MAOMp_sub)
            outResp_sub.append(resp_sub)

            if mode_ == "Normal":  # for normal runs
                if treatment == "control":
                    respPlot["soilControl"].append(respSoil / (0.8 * 24))
                    respPlot["subControl"].append(respSubstrate / (0.8 * 24))

                elif treatment == "exudates":
                    respPlot["soilExudates"].append(respSoil / (0.8 * 24))
                    respPlot["subExudates"].append(respSubstrate / (0.8 * 24))

                elif treatment == "leachates":
                    respPlot["soilLeachates"].append(respSoil / (0.8 * 24))
                    respPlot["subLeachates"].append(respSubstrate / (0.8 * 24))

            # make different output depending on the type of run
            # make different output depending on the type of run
        if mode_ == "Sensitivity":
            results_df.loc[len(results_df)] = [
                treatment,
                d,
                DOM_added / 0.8,  # change units from gC/m3 µgC/g soil
                DOM / 0.8,  # change units from gC/m3 µgC/g soil
                bact_DOM / 0.8,  # change units from gC/m3 µgC/g soil
                bact / 0.8,  # change units from gC/m3 µgC/g soil
                fungi / 0.8,  # change units from gC/m3 µgC/g soil
                respSubstrate
                / (0.8 * 24),  # change units from gC/m3/day to µg CO2-C/g soil/h
                baselineResp
                / (0.8 * 24),  # change units from gC/m3/day to µg CO2-C/g soil/h
                respSoil
                / (0.8 * 24),  # change units from gC/m3/day to µg CO2-C/g soil/h
                POM / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOMs / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOMp / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOM / (0.8 * 1000),
            ]  # change units from gC/m3 mgC/g soil

        if mode_ == "Bayesian":  # variables for which we have measured data
            results_df.loc[len(results_df)] = [
                treatment,
                d,
                resp / (0.8 * 24),  # change units from gC/m3/day to µg CO2-C/g soil/h
                resp_sub,
                POM / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOM / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                bact_total / 0.8,  # change units from gC/m3 µgC/g soil
                fungi / 0.8,  # change units from gC/m3 µgC/g soil
                POM_sub,
                MAOM_sub,
                bact_total,
                fungi_sub,
            ]

        if mode_ == "Normal":  # for normal runs
            results_df.loc[len(results_df)] = [
                treatment,
                d,
                DOM_added / 0.8,  # change units from gC/m3 µgC/g soil
                DOM / 0.8,  # change units from gC/m3 µgC/g soil
                bact_DOM / 0.8,  # change units from gC/m3 µgC/g soil
                bact / 0.8,  # change units from gC/m3 µgC/g soil
                fungi / 0.8,  # change units from gC/m3 µgC/g soil
                respSubstrate
                / (0.8 * 24),  # change units from gC/m3/day to µg CO2-C/g soil/h
                baselineResp
                / (0.8 * 24),  # change units from gC/m3/day to µg CO2-C/g soil/h
                respSoil
                / (0.8 * 24),  # change units from gC/m3/day to µg CO2-C/g soil/h
                POM / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOMs / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOMp / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOM / (0.8 * 1000),
            ]  # change units from gC/m3 mgC/g soil

        ############# end of daily run of coreMAOM   #############

        # after the total run is completed (after numDays)

        ############# Plotting   #############
        if Plotting:  # transform data for Plotting to adjusted units matching the data
            # change units to easily understandable for the plot
            # outDOMadded2 = np.divide(outDOMadded, 0.8) # change units from gC/m3 µgC/g soil
            # outBact_total2 = np.divide(outBact_total, 0.8) # change units from gC/m3 µgC/g soil
            # outbact_DOM2 = np.divide(outbact_DOM, 0.8)
            # outBact2 = np.divide(outBact, 0.8)
            # outFungi2 = np.divide(outFungi, 0.8)
            # outRespSubstrate2 = np.divide(outRespSubstrate, 0.8 * 24) # change units from gC/m3/day to µg CO2-C/g soil/h
            # outRespSoilBaseline2 = np.divide(outRespSoilBaseline, 0.8 * 24)
            # outRespSoil2 = np.divide(outRespSoil, 0.8 * 24)
            # outDOM2 = np.divide(outDOM, 0.8) # change units from gC/m3 µgC/g soil
            # outPOM2 = np.divide(outPOM, 0.8 * 1000) # change units from gC/m3 mgC/g soil
            # outMAOM2 = np.divide(outMAOM, 0.8 * 1000) # change units from gC/m3 mgC/g soil
            # outMAOMp2 = np.divide(outMAOMp, 0.8 * 1000) # change units from gC/m3 mgC/g soil
            # outMAOMs2 = np.divide(outMAOMs, 0.8 * 1000) # change units from gC/m3 mgC/g soil

            # first plot function
            def Dailyplot1(
                outDOMadded,
                outDOM,
                outbact_DOM,
                outBact,
                outFungi,
                outRespSubstrate,
                outRespSoil,
                outRespSoilBaseline,
                outPOM,
                outMAOMp,
                outMAOMs,
            ):  # plot in original KEYLINK units
                fig, ((p1, p2, p3), (p4, p5, p6)) = plt.subplots(
                    nrows=2, ncols=3, figsize=(12, 7)
                )  # was 10, 12
                fig.suptitle(treatment, size=16)
                fig.tight_layout(pad=2.0)
                ps = (p1, p2, p3, p4, p5, p6)
                # counter = count(0, 1)
                # columns = list(df)
                ps[0].set_title("DOM additions, µgC g-1 soil")
                ps[1].set_title("DOM, µgC g-1 soil")
                ps[2].set_title("Microbial biomass, µgC g-1 soil")
                ps[3].set_title("Respiration, µg C-CO2 g-1 soil h-1")
                ps[4].set_title("POM, mgC g-1 soil")
                # ps[4].set_title("SOM, mgC g-1 soil")
                ps[5].set_title("MAOM, mgC g-1 soil")

                p1.plot(time_d, outDOMadded)
                p2.plot(time_d, outDOM)
                p3.plot(time_d, outBact_total, label="bacteria")
                p3.plot(time_d, outbact_DOM, label="bacteria DOM feeding")
                p3.plot(time_d, outBact, label="bacteria only SOM feeding")
                p3.plot(time_d, outFungi, label="fungi")
                ps[2].legend(loc=(0.4, 0.03), shadow=True)  # loc='bottom right',

                p4.plot(time_d, outRespSubstrate, label="substrate-derived")
                p4.plot(time_d, outRespSoil, label="soil-derived incl. priming")
                p4.plot(
                    time_d,
                    outRespSoilBaseline,
                    label="soil-derived baseline",
                )
                ps[3].legend(loc=(0.25, 0.03), shadow=True)  # loc='bottom right',

                p5.plot(time_d, outPOM, label="POM")
                ps[4].legend(loc=(0.03, 0.03), shadow=True)  # loc='upper left',

                p6.plot(time_d, outMAOM, label="MAOM")
                p6.plot(time_d, outMAOMp, label="primary MAOM")
                p6.plot(time_d, outMAOMs, label="secondary MAOM")
                ps[5].legend(loc=(0.4, 0.03), shadow=True)  # loc='bottom right',

            # plot substrate-derived proportions
            def Dailyplot2(
                outBact_DOM_sub,
                outBact_sub,
                outFungi_sub,
                outDOM_sub,
                outPOM_sub,
                outMAOMs_sub,
                outMAOMp_sub,
            ):
                fig, (p1, p2, p3, p4) = plt.subplots(
                    nrows=4, ncols=1, figsize=(5, 14)
                )  # was 10,12
                fig.suptitle(treatment, size=16)
                fig.tight_layout(pad=3.0)
                ps = (p1, p2, p3, p4)
                # counter = count(0, 1)
                # columns = list(df)
                ps[0].set_title("substrate derived % of microbial pools")
                ps[1].set_title("substrate derived % of SOM pools")
                ps[2].set_title("substrate derived % of respiration")

                p1.plot(time_d, outBact_total_sub, label="bacteria")
                p1.plot(time_d, outBact_DOM_sub, label="bacteria DOM feeding")
                p1.plot(time_d, outBact_sub, label="bacteria SOM feeding")
                p1.plot(time_d, outFungi_sub, label="fungi")
                ps[0].legend(
                    loc="upper left", bbox_to_anchor=(1, 1), shadow=True
                )  # loc='bottom right',

                p2.plot(time_d, outDOM_sub, label="DOM")
                ps[1].legend(
                    loc="upper left", bbox_to_anchor=(1, 1), shadow=True
                )  # loc='bottom right',

                p3.plot(time_d, outPOM_sub, label="POM")
                p3.plot(time_d, outMAOM_sub, label="MAOM")
                p3.plot(time_d, outMAOMp_sub, label="primary MAOM")
                p3.plot(time_d, outMAOMs_sub, label="secondary MAOM")
                ps[2].legend(
                    loc="upper left", bbox_to_anchor=(1, 1), shadow=True
                )  # loc='bottom right',

                p4.plot(time_d, outResp_sub, label="total respiration")

            # check if there is a "figures" folder; if not, create one
            try:
                os.makedirs("./output/figures")
            except FileExistsError:
                # directory already exists
                pass

            # after each run, make a plot
            Dailyplot1(
                outDOMadded,
                outDOM,
                outbact_DOM,
                outBact,
                outFungi,
                outRespSubstrate,
                outRespSoil,
                outRespSoilBaseline,
                outPOM,
                outMAOMp,
                outMAOMs,
            )
            plt.savefig(".\output\\figures\Dailyplot1_" + treatment + ".png")

            Dailyplot2(
                outBact_DOM_sub,
                outBact_sub,
                outFungi_sub,
                outDOM_sub,
                outPOM_sub,
                outMAOMs_sub,
                outMAOMp_sub,
            )
            plt.savefig(
                ".\output\\figures\Dailyplot2_" + treatment + ".png",
                bbox_inches="tight",
            )
    ############# end of Plotting   #############

    # def drawRespPlot(respPlot):
    #     # count mean values of the modelled data
    #     soil_values_model = []
    #     soil_values_model.append(
    #         sum(respPlot["soilControl"]) / len(respPlot["soilControl"])
    #     )
    #     soil_values_model.append(
    #         sum(respPlot["soilLeachates"]) / len(respPlot["soilLeachates"])
    #     )
    #     soil_values_model.append(
    #         sum(respPlot["soilExudates"]) / len(respPlot["soilExudates"])
    #     )

    #     sub_values_model = []
    #     sub_values_model.append(
    #         sum(respPlot["subControl"]) / len(respPlot["subControl"])
    #     )
    #     sub_values_model.append(
    #         sum(respPlot["subLeachates"]) / len(respPlot["subLeachates"])
    #     )
    #     sub_values_model.append(
    #         sum(respPlot["subExudates"]) / len(respPlot["subExudates"])
    #     )

    #     labels = ["Control", "Leachates", "Exudates"]

    #     # measured values
    #     soil_values_measure = [26.23, 32.64, 27.64]
    #     sub_values_measure = [0, 6.83, 8.43]

    #     # create plot
    #     plt.figure(figsize=(8, 4))
    #     x = np.arange(len(labels))  # label locations
    #     width = 0.2  # width of the bars

    #     # create first subplot
    #     plt.subplot(1, 2, 1)
    #     plt.bar(x - width / 2, soil_values_model, width, label="Modeled", color="gray")
    #     plt.bar(
    #         x + width / 2, soil_values_measure, width, label="Measured", color="black"
    #     )
    #     plt.title("Soil derived")
    #     plt.ylabel("Respiration [µg C-CO2/g soil/h]")
    #     plt.xticks(x, labels)

    #     # create second subplot
    #     plt.subplot(1, 2, 2)
    #     plt.bar(x - width / 2, sub_values_model, width, label="Modeled", color="gray")
    #     plt.bar(
    #         x + width / 2, sub_values_measure, width, label="Measured", color="black"
    #     )
    #     plt.title("Substrate derived")
    #     plt.ylabel("")
    #     plt.xticks(x, labels)
    #     plt.legend(loc="upper left", bbox_to_anchor=(1, 1), shadow=True)

    #     plt.tight_layout()

    # drawRespPlot(respPlot)
    # plt.savefig("./output/figures/respPlot.png")

    # if (Bayesian):
    return results_df
