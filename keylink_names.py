"""
Names for the KEYLINK state vector: functional groups, organic matter pools and the daily flux
outputs. Use these instead of literal indices (B[10] -> B[SOM]).

State vector B (length NSTATE = 22), also the columns of KeylinkResults.pools:
    0-8   the nine functional groups (Flores et al. 2021, Fig. 2), all g C m-3
    9-12  litter, SOM, fine roots, cumulative CO2
    13-20 flux outputs of the current day (g C m-3 day-1), integrated alongside the pools and
          reset at the start of every day; 21 is a spare entry that stays 0
"""
from enum import IntEnum


class Pool(IntEnum):
    BACT = 0        # bacteria
    FUNGI = 1       # non-mycorrhizal fungi
    MYC = 2         # mycorrhizal fungi
    BVORES = 3      # bacterivores
    FVORES = 4      # fungivores
    SAP = 5         # detritivores (non-engineer saprotrophs)
    ENG = 6         # ecosystem engineers (engineer detritivores)
    HVORES = 7      # root herbivores
    PRED = 8        # predators
    LITTER = 9
    SOM = 10
    ROOTS = 11
    CO2 = 12        # cumulative respiration


class Flux(IntEnum):
    BACT_RESP = 13        # bacterial respiration
    FUN_RESP = 14         # fungal respiration
    MYC_RESP = 15         # mycorrhizal respiration
    BACT_GROWTH_SOM = 16  # bacterial growth on SOM
    BACT_GROWTH_LIT = 17  # bacterial growth on litter
    SOM_EATEN = 18        # SOM consumed by all consumers
    LIT_EATEN = 19        # litter consumed by all consumers (incl. faeces of detritivores/engineers)
    LIT_EATEN_ENG = 20    # litter consumed by engineers


# module level aliases, convenient in dense equation code
BACT, FUNGI, MYC, BVORES, FVORES, SAP, ENG, HVORES, PRED, LITTER, SOM, ROOTS, CO2 = Pool
BACT_RESP, FUN_RESP, MYC_RESP, BACT_GROWTH_SOM, BACT_GROWTH_LIT, SOM_EATEN, LIT_EATEN, LIT_EATEN_ENG = Flux

NGROUPS = 9           # functional groups
NPOOLS = 13           # groups + litter, SOM, roots, CO2
NSTATE = 22           # pools + 8 flux outputs + 1 spare

GROUPS = [p for p in Pool if p < NGROUPS]
GROUP_KEYS = ['bact', 'fungi', 'myc', 'bvores', 'fvores', 'sap', 'eng', 'hvores', 'pred']
POOL_KEYS = GROUP_KEYS + ['litter', 'SOM', 'roots', 'CO2']
FLUX_KEYS = ['bactResp', 'funResp', 'EMresp', 'bactGrowthSOM', 'bactGrowthLit', 'SOMeaten', 'LITeaten', 'LITeatenEng']

POOL_NAMES = ['Bacteria', 'Fungi', 'Mycorrhiza', 'bacterivores', 'fungivores', 'saprotrophs',
              'engineers', 'herbivores', 'predators', 'litter', 'SOM', 'roots', 'CO2']
POOL_COLOURS = ['blue', 'red', 'darkcyan', 'cyan', 'orange', 'purple', 'darkgreen', 'magenta',
                'black', 'brown', 'grey', 'chartreuse', 'yellow']

# pore size classes (Flores et al. 2021, 'Structural effects'), index into PV / PW arrays
INACC, BACTPORES, MICRO, MESO, MACRO = range(5)
PORE_KEYS = ['inaccessible', 'bacterial', 'micro', 'meso', 'macro']

# entries of the availability vector returned by keylink_functions.calcAvail
AVAIL_KEYS = ['SOMbact', 'SOMfungi', 'SOMeng', 'SOMsap', 'bbvores', 'ffvores', 'fvorespred',
              'bvorespred', 'hvorespred', 'sappred', 'engpred', 'SOMunavail']
