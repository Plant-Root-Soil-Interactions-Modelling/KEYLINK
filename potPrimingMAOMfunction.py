# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:08:20 2023

@author: gdeckmyn
"""
import keylink_functions as mf
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

import os


def run_model(AllParam, treatmentVar, mode_, Plotting, numDays, path):
    # output dataframe list
    bact_rhiz_rel = AllParam["bact_rhiz_rel"]
    fungi_rhiz_rel = AllParam["fungi_rhiz_rel"]
    DOM_EC = AllParam["DOM_EC"] #energy content of DOM J/gC
    # kpriming = AllParam["kpriming"]
    KSrhiz = AllParam["KSrhiz"]
    KSbulk = AllParam["KSbulk"]
    kPOM_MAOM = AllParam["kPOM_MAOM"]
    kMAOMs_MAOMp = AllParam["kMAOMs_MAOMp"]
    MAOMpmaxrate = AllParam["MAOMpmaxrate"]
    MAOMsmaxrate = AllParam["MAOMsmaxrate"]
    MAOMratioSP = AllParam["MAOMratioSP"]
    maxEffectMicMAOM = AllParam["maxEffectMicMAOM"]
    maxEffectSA_MAOM = AllParam["maxEffectSA_MAOM"]
    maxEffectN_MAOM = AllParam["maxEffectN_MAOM"]
    MM_N_MAOM = AllParam["MM_N_MAOM"]
    MM_Mic_MAOM = AllParam["MM_Mic_MAOM"]
    MM_SA_MAOM = AllParam["MM_SA_MAOM"]
    MM_DOM_MAOM = AllParam["MM_DOM_MAOM"]
    Priming_max = AllParam["Priming_max"]
    kpriming = 0
    # the ones we use from other calibration
    GMAXrhiz = AllParam[
        "GMAXrhiz"
    ]  # maximal growth rate for rhizeria [gC/(gC day)], KEYLINK was 1.24
    GMAXbulk = AllParam[
        "GMAXbulk"
    ]  # maximal growth rate for fungi [gC/(gC day)], KEYLINK
    mRecbulk = AllParam["mRecbulk"]  # how sensitive rhiz are to recalcitrance
    # resp=0.01 #respiration rate for rhizeria growing on DOM / ??do we really need a different one? it was set to 0 decided to ditch it and just the next one
    # =0.05  #respiration rate resp, [gC/(gC day)], KEYLINK
    DEATH = AllParam["DEATHrhiz"]  # death rate for rhizeria [gC/(gC day)], KEYLINK
    DEATHbulk = AllParam["DEATHbulk"]  # death rate for fungi [gC/(gC day)], KEYLINK
    pCN = AllParam[
        "pCN"
    ]  # sensitivity to CN ratio of consumed substrate, values 0-1, taken from KEYLINK (value for rhizeria)
    recMAOM = AllParam[
        "recMAOM"
    ]  # recalcitrance of MAOM, (recalcitrance of POM assumed 0)
    RESPrhiz = AllParam[
        "RESPrhiz"
    ]  # respiration rate of rhizeria, [gC/(gC day)], was 0.05 KEYLINK
    RESPbulk = AllParam[
        "RESPbulk"
    ]  # respiration rate of fungi, [gC/(gC day)], was 0.03 KEYLINK
    fSOM = AllParam["fSOM"] # what part of PrimingGrowth uses primed SOM as opposed to DOM, fraction 0-1
    T_MAXrhiz = AllParam["T_MAXrhiz"]
    T_MINrhiz = AllParam["T_MINrhiz"]
    T_OPTrhiz = AllParam["T_OPTrhiz"]
    T_MAXbulk = AllParam["T_MAXbulk"]
    T_MINbulk = AllParam["T_MINbulk"]
    T_OPTbulk = AllParam["T_OPTbulk"]
    Q10rhiz = AllParam["Q10rhiz"]
    Q10bulk = AllParam["Q10bulk"]

    # input parameters that do not change (=measurable) and are not calibrated, just 'start situation"

    # those which will change for different runs
    DOMinput = treatmentVar[
        "DOMinput"
    ]  # 10 #DOM added in each addition [gC/m3], Jílková2022 >> in Jilkova2024, it was 0.5 mgC, therefore it is 5 gC/m3
    CN_DOMinput = treatmentVar[
        "CN_DOMinput"
    ]  #  #CN of the daily input [unitless], Jílková2022: leachates 80, exudates 6
    CN_MAOMp = treatmentVar["CN_MAOMp"]  # 15 #☼ assumed constant
    CN_POM = treatmentVar["CN_POM"]  # 24 #CN of SOM, Jílková2022
    pH = treatmentVar["pH"]  # 4.1 #Jílková2022
    temp = treatmentVar["temp"]  # 21
    treatment = treatmentVar["treatment"] #name of the treatment
    treatmentID = treatmentVar["treatmentID"] #number of the treatment
    fClay = treatmentVar["fClay"] #  # weight fraction [g/g], 0.17 for Jílková2022
    fSilt = treatmentVar["fSilt"] # weight fraction [g/g], 0.24 for Jílková2022

    # those same for Jílková 2022 and experiment 2024
     
    d_freq = (
        14  # how often is substrate added, every x days, is 14 for Jílková2022 and 2024
    )
    # numDays = 161  # number of days of incubation experiment/how long to run the model, 155 in Jílková2022, 161 in Jílková 2024

    # those that will be the same for all 16 runs
    BD = 800  # bulk density [kg/m³]
    claySA = 800000  # surface area of clay [m²/kg] was 8000000 cm²/g
    CN_bact = 4  # CN of rhizosphere microbes, from KEYLINK, in Jílková2022 initial CN of microbial biomass is 10
    CN_fungi = 8  # KEYLINK
   
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
    # Plotting = False  # flag 1 to enable making of plots, so that this can be turned off during sensitivity analysis etc.
    # Bayesian = True  # flag 1 if performing Bayesian
    # Sensitivity = False

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
            "rhiz_DOM",
            "rhiz",
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
        # treatment, d, resp, resp_sub, POM, MAOM, rhiz_total, fungi, POM_sub, MAOM_sub, rhiz_total_sub, fungi_sub]
        column_names = [
            "treatment",
            "day",
            "respSoil",
            "respSubstrate",
            "DOMSoil",
            "DOMSubstrate",
            "POMSoil",
            "POMSubstrate",
            "MAOMSoil",
            "MAOMSubstrate",
            "CN_DOM",
            "CN_MAOM",
            "bactSoil",
            "bactSubstrate",
            "fungiSoil",
            "fungiSubstrate"
            ]
        results_df = pd.DataFrame(columns=column_names)
        
    if mode_ == "Normal":  # normal runs

        column_names = [
            "treatment",
            "day", 
            "DOM_added",
            "respSoil",
            "respSubstrate",
            "DOMSoil",
            "DOMSubstrate",
            "POMSoil",
            "POMSubstrate",
            "MAOMSoil",
            "MAOMSubstrate",
            "CN_DOM",
            "CN_MAOM",
            "bactSoil",
            "bactSubstrate",
            "fungiSoil",
            "fungiSubstrate"
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
        outrhiz = []
        outbulk = []
        outRespSubstrate = []
        outRespSoil = []
        outRespSoilBaseline = []
        outrhiz_sub = []
        outbulk_sub = []
        outDOM_sub = []
        outPOM_sub = []
        outMAOM_sub = []
        outMAOMs_sub = []
        outMAOMp_sub = []
        outResp_sub = []

    # *************************************************************************
    # initializing variables (what changes during run)

    # variables that will be initialized differently for different runs
    bact = treatmentVar[
        "bact"
    ]  # total biomass of rhizeria [gC/m3], was 6 final noadd average from PLFA from Jílková2022
    CN_MAOMs = treatmentVar[
        "CN_MAOMsini"
    ]  # estimated but we don't know the true value, assumed to vary with CN_DOM
    print(treatmentVar["treatment"])
    # if treatmentVar["treatmentID"] == 2:
    #     return
    
    print(CN_MAOMs, "CN_MAOMs")
    fungi = treatmentVar[
        "fungi"
    ]  # biomass of fungi [gC/m3] based on final noadd in Jílková et al. 2022
    DOM = treatmentVar["DOMini"]   # DOM [gC/m3]
    CN_DOM = treatmentVar["CN_DOMini"] 
    POM = treatmentVar[
        "POMini"
    ]  # C in POM [gC/m3], calculated as initialSOM-MAOM using initialSOM from Jílková2022
    MAOM = treatmentVar["MAOMini"]  # C in MAOM [gC/m3] average noAdd Jílková2022
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
    # biomass of rhizosphere microbes [gC/m3]
    bact_rhiz = bact * bact_rhiz_rel

    fungi_rhiz = fungi * fungi_rhiz_rel 
    rhiz = bact_rhiz +  fungi_rhiz
    rhiz_sub = 0  # proportion of rhizosphere microbial carbon that is substrate derived in contrast to soil-derived / values 0 to 1/
    # biomass of bulk soil microbes [gC/m3]
    bact_bulk = bact * (1 - bact_rhiz_rel)
    fungi_bulk = fungi * (1 - fungi_rhiz_rel)
    bulk = bact_bulk + fungi * fungi_bulk 
    bulk_sub = 0  # proportion of this carbon in microbes that is substrate derived in contrast to soil-derived / values 0 to 1/
    # print('bact_rhiz', 'fungi_rhiz', 'bact_bulk', 'fungi_bulk', bact_rhiz, fungi_rhiz, bact_bulk, fungi_bulk)
    
    #CN ratio for rhiz and bulk based on the proportion of bacterial and fungal biomass
    CN_rhiz = (bact_rhiz*CN_bact + fungi_rhiz*CN_fungi)/rhiz
    CN_bulk = (bact_bulk*CN_bact + fungi_bulk*CN_fungi)/bulk
    FB_rhiz = fungi_rhiz/bact_rhiz
    FB_bulk = fungi_bulk/bact_bulk
    # print('FB_rhiz', FB_rhiz,
    #       'FB_bulk', FB_bulk)
    DOM_sub = 0  # relative substrate derived C in DOM /values 0 to 1/, portion of DOM carbon that is substrate derived in contrast to soil-derived / values 0 to 1/ is a ratio between substrate-derived C and total C in DOM

    DOM_N = DOM/CN_DOM  # calculate DOM N from initial DOM and its CN

    # fungi_sub = 0  # proportion of fungal carbon that is substrate derived in contrast to soil-derived / values 0 to 1/

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
    # AllC = DOM + POM + MAOM + rhiz_total + fungi + resp - DOMadded
    # print(AllC)

    # function coreMAOM
    # def coreMAOM (Bayesian, Sensitivity):

    for d in range(numDays):
        # print("day", d)
        # on day 0 and then every 14 days, add DOM
        # if treatmentID == 5:
        #     print(treatmentID, "CN_DOM in the beginning of day: ", CN_DOM)

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

                # if treatmentID == 5:
                #     print(
                #         treatmentID,
                #         "CN_DOM recalculated in the beginning of input day: ",
                #         CN_DOM,
                #     )
        else:
            DOM_added = 0
        # saturation of MAOMs depends on amount of MAOMp so recalculated every day
        maxMAOMs = MAOMp * MAOMratioSP  # maximum primary MAOM

        # if maxMAOMs < 0:
        #     print("MAOMp2: ", MAOMp)

        # find t modifier
        modtrhiz = mf.calcmodt(temp, T_OPTrhiz, T_MINrhiz, T_MAXrhiz) #temp modifier of growth, rhizosphere microbes
        modtbulk = mf.calcmodt(temp, T_OPTbulk, T_MINbulk, T_MAXbulk) #temp modifier of growth, bulk microbes
        rRESPrhiz = mf.calcresp(temp, T_OPTrhiz, RESPrhiz, Q10rhiz) #respiration rate as modified by temperature
        rRESPbulk = mf.calcresp(temp, T_OPTbulk, RESPbulk, Q10bulk)
        # print(rRESPrhiz, rRESPfungi)

        # microbial growth on DOM and priming, only susing MAOMs
        if CN_DOM > 0:
            (
                DOM,
                DOM_sub,
                DOM_N,
                CN_DOM,
                rhiz,
                rhiz_sub,
                POM,
                MAOMs,
                MAOMp,
                respDOM,
                respDOM_sub,
                respPriming,
                respPriming_sub,
            ) = mf.calcRhizosphere(
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
                rhiz,
                rhiz_sub,
                CN_rhiz,
                DOM,
                DOM_sub,
                CN_DOM,
                GMAXrhiz,
                DEATH,
                pCN,
                pH,
                rRESPrhiz,
                KSrhiz,
                DOM_EC,
                Priming_max,
                kpriming,
                kPOM_MAOM,
                kMAOMs_MAOMp,
                modtrhiz,
                fSOM
            )
        # print('calc.Rhizo')
        #               if (MAOMs<0):
        #                   print('mainLine270 DOM, rhiz, fungi, MAOMs, MAOMp', DOM,rhiz, fungi, MAOMs, MAOMp)
        else:

            respDOM = 0
            respDOM_sub = 0
            respPriming = 0
            respPriming_sub = 0

        # if treatmentID == 5:
        #     print(treatmentID, "CN_DOM after calcRhizosphere: ", CN_DOM)
        resp = respDOM + respPriming
        # resp_all += resp
        #rhiz_total = rhiz_DOM + rhiz
        MAOM = MAOMs + MAOMp
        # AllC = DOM + POM + MAOM + rhiz_total + fungi + resp_all - DOM_added

        # MAOM formation
        MicrobialC = (
            rhiz + bulk
        )  # all microbes contribute to MAOM formation
        #print("line 448 CN_MAOMs", CN_MAOMs) 
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
                maxEffectMicMAOM,
                MM_Mic_MAOM,
                maxEffectN_MAOM,
                MM_N_MAOM,
                maxEffectSA_MAOM,
                MM_SA_MAOM,
                CN_MAOMp,
                CN_MAOMs,
            )
                
        # if treatmentID == 5:
        #     print(treatmentID, "CN_DOM after calcMAOM: ", CN_DOM)

        else: #CN_DOM negative
                CN_DOM=CN_DOM  
            
                #print("line 485 CN_MAOMs", CN_MAOMs) 
        MAOM = MAOMs + MAOMp

        # bulk soil microbial growth on SOM (without substrate DOM additions)
        availability = mf.calcAvailPot(
            PV, PW
        )  # calculates availability of SOM decomposition by rhizeria and fungi, separately, from pore size distribution and soil water
        # calculate maximal growth (gmax) for rhizeria/fungi on POM/MAOM separately

        # if CN_MAOMs <= 0:
        #     print("CN_MAOMs: ", CN_MAOMs)
        
        #first calculate maximum growth on DOM if it was unlimited, both for bulk microbes
        # do this only if there is some non-zero DOM, not to run into problems with dividing by zero
   
        
        if DOM > 0:
            
            gmaxbDOM = (
             mf.calcgmaxmod(CN_bulk, CN_DOM, pCN, 0.0, 0, pH, 1) * GMAXbulk
             ) # gmax for bulk on DOM
            
            # get resp from DOM and reduce avaialabilty
            if availability[0]* bulk > rRESPbulk * bulk:
                respDOM+=rRESPbulk * bulk
                respDOMbulk = rRESPbulk * bulk
                avail=availability[0]-rRESPbulk
                respRest=0
            else:
                respDOM=availability[0]* bulk
                avail=0
                respRest=rRESPbulk * bulk - respDOM
            #calculate realized growth on DOM (this is actually assimilation, not growth)
            bulkDOMgrowth = modtbulk * mf.calcgrowth(
                bulk, DOM, avail, gmaxbDOM, KSbulk * bulk
            )
            
        else:
            bulkDOMgrowth = 0
         
   
        #then to ensure that the sum of gmaxes from different substrates does not exceed GMAX, 
        #reduce GMAX accordingly by what growth was already realized from previous substrates
        #print("518 POM bulk CN_POM", CN_POM)
        gmaxbPOM = (
            mf.calcgmaxmod(CN_bulk, CN_POM, pCN, 0.0, 0, pH, 1) * (GMAXbulk - bulkDOMgrowth)
        )  # gmax for rhiz on POM
        
        if availability[0]* POM > respRest:
                respPOM=respRest
                avail=availability[0]-respRest/POM
                respRest=0
        else:
                respPOM=availability[0]* POM
                avail=0
                respRest=respRest-respMAOMs
        #calculate realized growth on POM
        bulkPOMgrowth = modtbulk * mf.calcgrowth(
            bulk, POM, avail, gmaxbPOM, KSbulk * bulk
        )
        
        # we assume MAOMp can only be lost through priming, so normal growth uses MAOMs
        #also reduce gmax by what was already grown on DOM and POM
        #print("530 MAOM bulk CN_MAOMs", CN_MAOMs)
        gmaxbMAOM = (
            mf.calcgmaxmod(CN_bulk, CN_MAOMs, pCN, recMAOM, mRecbulk, pH, 1) * (GMAXbulk - bulkDOMgrowth - bulkPOMgrowth)
        )  # gmax for rhiz on MAOM
        
        if availability[0]* MAOMs > respRest:
                respMAOMs=respRest
                avail=availability[0]-respRest/MAOMs
                respRest=0
        else:
                respMAOMs=availability[0]* MAOMs
                avail=0
                respRest=respRest-respMAOMs
        #calculate realized growth on MAOM
        bulkMAOMgrowth = modtbulk * mf.calcgrowth(
            bulk, MAOMs, avail, gmaxbMAOM, KSbulk * bulk
        )
       
        # print('GMAX', GMAX,
              # "\nrhizDOMgrowth", rhizDOMgrowth,
              # "\nrhizPOMgrowth", rhizPOMgrowth,
              # '\nrhizMAOMgrowth', rhizMAOMgrowth,      
              # '\nGMAXfungi', GMAXfungi,
              # "\nfungiDOMgrowth", fungiDOMgrowth, 
              # '\nfungiPOMgrowth', fungiPOMgrowth,               
              # '\nfungiMAOMgrowth', fungiMAOMgrowth) 
        
        # recalculate substrate derived C in bulk soil microbes and C pools
        DOM_sub_abs = DOM * DOM_sub  # recalculate because changesin calc.Rhizosphere
        POM_sub_abs = POM * POM_sub  # recalculate because changes in calc.Rhizosphere
        MAOMs_sub_abs = (
            MAOMs * MAOMs_sub
        )  # recalculate because changes in calc.Rhizosphere and calc.MAOM
        bulk_sub_abs = (
            bulk * bulk_sub
        )  # absolute substrate derived C in bulk [gC/m3]
       
        
        #calculate the overall change in bulk micrbial biomass, only if there was not enough for respiration this has become death
        dbulk = bulkDOMgrowth + bulkPOMgrowth + bulkMAOMgrowth - DEATHbulk * bulk - respRest

        

        #the consequent changes in the pools being eaten for growth or respired
        DOM += - bulkDOMgrowth + DEATHbulk * bulk + respRest - respDOMbulk  # add dead bulk to DOM and death from no C to resp
        POM += -bulkPOMgrowth - respPOM  # subtract what has been eaten from POM to grow and to respire
        MAOMs += -bulkMAOMgrowth - respMAOMs   # and MAOMs

        # update CN DOM
        DOM_N += - bulkDOMgrowth / CN_bulk + DEATHbulk * bulk / CN_bulk + respRest / CN_bulk

        CN_DOM = DOM / DOM_N  # recalculate CN DOM

        # if treatmentID == 5:
        #     print(treatmentID, "CN_DOM in the end of the day: ", CN_DOM)

        #    if (-drhiz>rhiz):
        #        print('mainLine307  rhiz, rhizPOMgrowth, POM, rhizMAOMgrowth, DEATH*rhiz, rRESPrhiz*rhiz', rhiz, rhizPOMgrowth, POM, rhizMAOMgrowth, DEATH*rhiz, rRESPrhiz*rhiz)

        DOM_sub_abs += (
            - bulkDOMgrowth * DOM_sub + DEATHbulk * bulk * bulk_sub+ respRest* bulk_sub - respDOMbulk * DOM_sub
        )  # add corresponding part of substrate derived C to DOM
        POM_sub_abs -= bulkPOMgrowth * POM_sub
        MAOMs_sub_abs -= bulkMAOMgrowth * MAOMs_sub
        
        bulk_sub_abs += (
            bulkDOMgrowth * DOM_sub
            + bulkMAOMgrowth * MAOMs_sub
            + bulkPOMgrowth * POM_sub
            - DEATHbulk * bulk * bulk_sub
            - rRESPbulk * bulk * bulk_sub
        )  # add the corresponding part of growth on MAOM as substrate derived C, subtract correspodning part of death and respiration

        baselineRespbulk = rRESPbulk * bulk - respRest
        baselineRespbulk_sub_abs = (
            baselineRespbulk * bulk_sub
        )  # what part of this respiration is substrate derived
        baselineRespbulk_sub = baselineRespbulk_sub_abs / baselineRespbulk
        bulk += dbulk

        
        # update relative substrate derived C proportions
        DOM_sub = DOM_sub_abs / DOM  # relative substrate derived C in DOM
        POM_sub = POM_sub_abs / POM  # relative substrate derived C in DOM
        MAOMs_sub = MAOMs_sub_abs / MAOMs  # relative substrate derived C in DOM
        
        bulk_sub = (
            bulk_sub_abs / bulk
        )  # update relative substrate derived C in bulk soil microbes
        # add up things
        MAOM = MAOMp + MAOMs
        # bact_total = rhiz*bact_rhiz_rel + bulk*bact_bulk_rel
        # baseline respiration without priming
        baselineResp = (
            baselineRespbulk + respDOM
        )  # of course this respDOM is higher if previous day DOM-feeding rhizeria grew more because of priming
        # all respiration
        resp = baselineResp + respPriming
        

        # calculate average substrate proportions
        MAOM_sub = MAOMp_sub * (MAOMp / MAOM) + MAOMs_sub * (
            MAOMs / MAOM
        )  # average substrate proportion in MAOM

        resp_sub = (
            baselineRespbulk_sub * (baselineRespbulk / resp)
            + respDOM_sub * (respDOM / resp)
            + respPriming_sub * (respPriming / resp)
        )
        # print ("fractions resp",baselineRespbulk / resp,respDOM / resp,respPriming / resp)
        # print ("fractions resp sub",baselineRespbulk_sub,respDOM_sub,respPriming_sub)
        
        
        respSubstrate = resp_sub * resp  # substrate derived respiration (absolute)
        respSoil = resp - respSubstrate  # soil-derived respiration (absolute)
        
        #soil-derived respiration from all sources except respPriming     
        #respDOM - respiration of DOM feeding rhizeria without priming being activ
        respSoilBaseline = respDOM * (1-respDOM_sub) + baselineRespbulk * (1-baselineRespbulk_sub) 
        # AllC = DOM + POM + MAOM + rhiz_total + fungi + resp - DOMadded
        # print("Potprim function line 635")
        #print("fractions resp: baselineRespbulk, respDOM, respPriming", baselineRespbulk/resp, respDOM/resp, respPriming/resp)
        #print("fractions resp_sub:", baselineRespbulk_sub, respDOM_sub, respPriming_sub)
        # print("")
        #soil and substrate derived C pools
        DOMSubstrate = DOM_sub  * DOM
        DOMSoil = DOM - DOMSubstrate
        POMSubstrate = POM_sub  * POM
        POMSoil = POM - POMSubstrate
        MAOMSubstrate = MAOM_sub  * MAOM
        MAOMSoil = MAOM - MAOMSubstrate
        
        #calculating fungi and bacteria back, using fixed FB ratios of rhizosphere and bulk soil
        bact_rhiz = rhiz/(FB_rhiz + 1) 
        fungi_rhiz = bact_rhiz * FB_rhiz
        bact_bulk = bulk/(FB_bulk + 1) 
        fungi_bulk = bact_bulk * FB_bulk
        bact = bact_rhiz + fungi_rhiz
        fungi = bact_bulk + fungi_bulk
        
        #calculating substrate derived proportion in bacteria and fungi               
        # bact_sub = 0
        # fungi_sub = 0        
        #first calculate in absolute units amount of substrate derived rhizosphere microbes
        rhiz_sub_abs = rhiz * rhiz_sub
        bulk_sub_abs = bulk * bulk_sub # same for bulk soil microbes
        
        # use the same equations as above when calculating the total four microbial pools, just use the substrate derived proportions instead
        bact_rhiz_sub_abs = rhiz_sub_abs/(FB_rhiz + 1) 
        fungi_rhiz_sub_abs = bact_rhiz_sub_abs * FB_rhiz
        bact_bulk_sub_abs = bulk_sub_abs/(FB_bulk + 1) 
        fungi_bulk_sub_abs = bact_bulk_sub_abs * FB_bulk
        bact_sub_abs = bact_rhiz_sub_abs + fungi_rhiz_sub_abs
        fungi_sub_abs = bact_bulk_sub_abs + fungi_bulk_sub_abs
        bact_sub = bact_sub_abs/bact
        fungi_sub = fungi_sub_abs/fungi
        
        #calculate soil and substrate derived bacteria and fungi
        bactSubstrate = bact_sub  * bact
        bactSoil = bact - bactSubstrate
        fungiSubstrate = fungi_sub  * fungi
        fungiSoil = fungi - fungiSubstrate
        
        #todo check when CN_MAOM changes, calculate CN MAOM
        CN_MAOM = (MAOMs * CN_MAOMs + MAOMp * CN_MAOMp)/MAOM
        
        # if treatmentID == 5 or treatmentID == 1 or treatmentID == 3:
        #     print(
        #         treatmentID,
        #         "respSubstrate: ",
        #         respSubstrate,
        #         "resp_sub: ",
        #         resp_sub,
        #         "resp: ",
        #         resp,
        #         "baselineResp: ",
        #         baselineResp,
        #     )

        # if treatmentID == 5 or treatmentID == 1 or treatmentID == 3:
        #     print(treatmentID, respPriming)

        if Plotting:  # save data for Plotting
            outtreatment.append(treatment)
            time_d.append(d)  # store days in an array for plotting
            outDOMadded.append(DOM_added / 0.8)  # change units from gC/m3 µgC/g soil
            outMAOM.append(MAOM / (0.8 * 1000))  # change units from gC/m3 mgC/g soil
            outMAOMp.append(MAOMp / (0.8 * 1000))  # change units from gC/m3 mgC/g soil)
            outMAOMs.append(MAOMs / (0.8 * 1000))  # change units from gC/m3 mgC/g soil)
            outPOM.append(POM / (0.8 * 1000))  # change units from gC/m3 mgC/g soil
            outDOM.append(DOM / 0.8)  # change units from gC/m3 µgC/g soil)
            outrhiz.append(rhiz / 0.8)  # change units from gC/m3 µgC/g soil  
            outbulk.append(bulk / 0.8)  # change units from gC/m3 µgC/g soil            
            # outrhiz.append(rhiz / 0.8)  # change units from gC/m3 µgC/g soil
            # outFungi.append(fungi / 0.8)  # change units from gC/m3 µgC/g soil
            outRespSubstrate.append(
                respSubstrate / (0.8 * 24)
            )  # change units from gC/m3/day
            outRespSoilBaseline.append(
                respSoilBaseline / (0.8 * 24)
            )  # change units from gC/m3/day
            outRespSoil.append(respSoil / (0.8 * 24))  # change units from gC/m3/day
            # substrate-derived %
            outrhiz_sub.append(rhiz_sub)
            outbulk_sub.append(bulk_sub)
            outDOM_sub.append(DOM_sub)
            outPOM_sub.append(POM_sub)
            outMAOM_sub.append(MAOM_sub)
            outMAOMs_sub.append(MAOMs_sub)
            outMAOMp_sub.append(MAOMp_sub)
            outResp_sub.append(resp_sub)

            # make different output depending on the type of run
            # make different output depending on the type of run
        if mode_ == "Sensitivity":
            results_df.loc[len(results_df)] = [
                treatment,
                d,
                DOM_added / 0.8,  # change units from gC/m3 µgC/g soil
                DOM / 0.8,  # change units from gC/m3 µgC/g soil
                rhiz / 0.8,  # change units from gC/m3 µgC/g soil
                bulk / 0.8,  # change units from gC/m3 µgC/g soil
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
                respSoil / (0.8 * 24),     # change units from gC/m3/day to µg CO2-C/g soil/h
                respSubstrate / (0.8 * 24),     # change units from gC/m3/day to µg CO2-C/g soil/h
                DOMSoil / (0.8 * 1000), # change units from gC/m3 mgC/g soil
                DOMSubstrate / (0.8 * 1000),
                POMSoil / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                POMSubstrate / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOMSoil / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOMSubstrate / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                CN_DOM,
                CN_MAOM,
                bactSoil / 0.8,# change units from gC/m3 µgC/g soil
                bactSubstrate / 0.8,# change units from gC/m3 µgC/g soil
                fungiSoil / 0.8, # change units from gC/m3 µgC/g soil
                fungiSubstrate / 0.8, # change units from gC/m3 µgC/g soil
                ]
         

        if mode_ == "Normal":  # for normal runs
            results_df.loc[len(results_df)] = [
                treatment,
                d,
                DOM_added / 0.8,  # change units from gC/m3 µgC/g soil
                respSoil / (0.8 * 24),     # change units from gC/m3/day to µg CO2-C/g soil/h
                respSubstrate / (0.8 * 24),     # change units from gC/m3/day to µg CO2-C/g soil/h
                DOMSoil / (0.8 * 1000), # change units from gC/m3 mgC/g soil
                DOMSubstrate / (0.8 * 1000),
                POMSoil / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                POMSubstrate / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOMSoil / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                MAOMSubstrate / (0.8 * 1000),  # change units from gC/m3 mgC/g soil
                CN_DOM,
                CN_MAOM,
                bactSoil / 0.8,# change units from gC/m3 µgC/g soil
                bactSubstrate / 0.8,# change units from gC/m3 µgC/g soil
                fungiSoil / 0.8, # change units from gC/m3 µgC/g soil
                fungiSubstrate / 0.8, # change units from gC/m3 µgC/g soil
               
            ]  

        ############# end of daily run of coreMAOM   #############

        # after the total run is completed (after numDays)

        ############# Plotting   #############
    if Plotting:  # transform data for Plotting to adjusted units matching the data
        # change units to easily understandable for the plot
        # outDOMadded2 = np.divide(outDOMadded, 0.8) # change units from gC/m3 µgC/g soil
        # outrhiz_total2 = np.divide(outrhiz_total, 0.8) # change units from gC/m3 µgC/g soil
        # outrhiz_DOM2 = np.divide(outrhiz_DOM, 0.8)
        # outrhiz2 = np.divide(outrhiz, 0.8)
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
            outrhiz,
            outbulk,
            outRespSubstrate,
            outRespSoil,
            outRespSoilBaseline,
            outPOM,
            outMAOMp,
            outMAOMs,
            treatment,
            path,
        ):  # plot in original KEYLINK units

            fig, ((p1, p2, p3), (p4, p5, p6)) = plt.subplots(
                nrows=2, ncols=3, figsize=(12, 8)
            )  # was 10, 12
            fig.suptitle(treatment, size=16)
            fig.tight_layout(pad=2.0)
            ps = (p1, p2, p3, p4, p5, p6)
            plt.subplots_adjust(bottom=0.2, hspace=0.6)

            # counter = count(0, 1)
            # columns = list(df)
            ps[0].set_title("DOM additions, µgC g-1 soil")
            ps[1].set_title("DOM, µgC g-1 soil")
            ps[2].set_title("POM, mgC g-1 soil")
            ps[3].set_title("Microbial biomass, µgC g-1 soil")
            ps[4].set_title("Respiration, µg C-CO2 g-1 soil h-1")
            # ps[4].set_title("SOM, mgC g-1 soil")
            ps[5].set_title("MAOM, mgC g-1 soil")

            p1.plot(time_d, outDOMadded, label="DOM")
            ps[0].legend(
                loc="upper left", bbox_to_anchor=(0, -0.15), shadow=True
            )  # loc='upper left',

            p2.plot(time_d, outDOM, label="DOM")
            ps[1].legend(
                loc="upper left", bbox_to_anchor=(0, -0.15), shadow=True
            )  # loc='upper left',

            p3.plot(time_d, outPOM, label="POM")
            ps[2].legend(
                loc="upper left", bbox_to_anchor=(0, -0.15), shadow=True
            )  # loc='upper left',

            p4.plot(time_d, outrhiz, label="rhizosphere microbes")
            p4.plot(time_d, outbulk, label="bulk soil microbes")
            ps[3].legend(
                loc="upper left", bbox_to_anchor=(0, -0.15), shadow=True
            )  # loc='bottom right',

            p5.plot(time_d, outRespSubstrate, label="substrate-derived")
            p5.plot(time_d, outRespSoil, label="soil-derived incl. priming")
            p5.plot(
                time_d,
                outRespSoilBaseline,
                label="soil-derived baseline",
            )
            ps[4].legend(
                loc="upper left", bbox_to_anchor=(0, -0.15), shadow=True
            )  # loc='bottom right',

            p6.plot(time_d, outMAOM, label="MAOM")
            p6.plot(time_d, outMAOMp, label="primary MAOM")
            p6.plot(time_d, outMAOMs, label="secondary MAOM")
            ps[5].legend(
                loc="upper left", bbox_to_anchor=(0, -0.15), shadow=True
            )  # loc='bottom right',
            
           
            plt.savefig(os.path.join(figures_path, "Dailyplot1_" + treatment + ".png"))
            plt.close()

        # plot substrate-derived proportions
        def Dailyplot2(
            outrhiz_sub,
            outbulk_sub,
            outDOM_sub,
            outPOM_sub,
            outMAOMs_sub,
            outMAOMp_sub,
            treatment,
            figures_path,
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
            ps[3].set_title("total respiration")

            p1.plot(time_d, outrhiz_sub, label="rhizosphere microbes")
            p1.plot(time_d, outbulk_sub, label="bulk soil microbes")
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

            p4.plot(time_d, outResp_sub, label="respiration")
            ps[3].legend(
                loc="upper left", bbox_to_anchor=(1, 1), shadow=True
            )  # loc='bottom right',

            plt.savefig(
                os.path.join(figures_path, "Dailyplot2_" + treatment + ".png"),
                bbox_inches="tight",
            )
            plt.close()

        # # check if there is a "figures" folder; if not, create one
        # try:
        #     os.makedirs("./output/figures")
        # except FileExistsError:
        #     # directory already exists
        #     pass
    
# plot in normal (not KEYLINK) units
        # after each run, make a plot
        
        #define the folder where figures should be saved as output/figures
        figures_path = os.path.join(path, "figures")
        
        Dailyplot1(
            outDOMadded,
            outDOM,
            outrhiz,
            outbulk,
            outRespSubstrate,
            outRespSoil,
            outRespSoilBaseline,
            outPOM,
            outMAOMp,
            outMAOMs,
            treatment,
            figures_path,
        )
        
        
        
        Dailyplot2(
            outrhiz_sub,
            outbulk_sub,
            outDOM_sub,
            outPOM_sub,
            outMAOMs_sub,
            outMAOMp_sub,
            treatment,
            figures_path,
        )

    ############# end of Plotting   #############

    # if (Bayesian):
    return results_df
