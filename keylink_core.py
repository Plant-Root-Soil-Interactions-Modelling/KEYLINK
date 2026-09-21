'''
KEYLINK core: configuration, daily time loop and results.

Created on 27.07.2016 (A. Schnepf - G. Deckmyn - G. Cauwenberg - O. Flores)

Usage
-----
    import keylink_core as core
    cfg     = core.KeylinkConfig.from_files()          # reads the KL_*.txt parameter files of the cwd
    model   = core.KeylinkModel(cfg)
    results = model.run(cfg.gmax)                      # or any other gmax vector (calibration)
    results.pools[:, core.SOM]                          # daily values, g C m-3

Importing this module has no side effects: nothing is read and no state is kept between runs.
All parameters live in a KeylinkConfig (a frozen dataclass; use cfg.replace(...) for variants),
every run allocates its own output arrays and returns them in a KeylinkResults object.
The model equations themselves are in keylink_functions (fCompSpecies is the ODE right-hand side).
'''
import dataclasses
from dataclasses import dataclass, field
import os

import numpy as np
from scipy.integrate import odeint

import keylink_functions as mf

# ---------------------------------------------------------------------------------------------
# pool indices of the state vector B / the daily output (Flores et al. 2021, Fig. 2)
BACT, FUNGI, MYC, BVORES, FVORES, SAP, ENG, HVORES, PRED = range(9)   # the nine functional groups
LITTER, SOM, ROOTS, CO2 = 9, 10, 11, 12                                # organic matter pools, cumulative CO2
# daily flux outputs carried along in the state vector (reset every day)
BACT_RESP, FUN_RESP, MYC_RESP, BACT_GROWTH_SOM, BACT_GROWTH_LIT, SOM_EATEN, LIT_EATEN, LIT_EATEN_ENG = range(13, 21)
NGROUPS = 9          # functional groups
NPOOLS = 13          # groups + litter, SOM, roots, CO2
NSTATE = 22          # pools + 8 flux outputs + 1 spare (always 0)

POOL_NAMES = np.array(['Bacteria', 'Fungi', 'Mycorrhiza', 'bacterivores', 'fungivores', 'saprotrophs',
                       'engineers', 'herbivores', 'predators', 'litter', 'SOM', 'roots', 'CO2'])
POOL_COLOURS = np.array(['blue', 'red', 'darkcyan', 'cyan', 'orange', 'purple', 'darkgreen', 'magenta',
                         'black', 'brown', 'grey', 'chartreuse', 'yellow'])

DAYS_IN_MONTH = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])


def import_pools(filename, directory='.'):
    """Load an array from a KL_*.txt parameter file (the '.txt' is added)."""
    return np.loadtxt(os.path.join(directory, filename + '.txt'), encoding='latin1')


def export_pools(filename, array, directory='.'):
    """Save an array in a readable format, compatible with R."""
    np.savetxt(os.path.join(directory, filename + '.txt'), array)


def _is_leap_year(year):
    return year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)


# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class KeylinkConfig:
    """All parameters and inputs of one KEYLINK set-up. Arrays indexed by functional group have
    length 9 in the order bact, fungi, myc, bacterivores, fungivores, detritivores, engineers,
    herbivores, predators. Build it with KeylinkConfig.from_files() and derive variants with
    cfg.replace(field=value, ...)."""
    # --- functional groups (KL_FaunalParams.txt)
    gmax: np.ndarray          # maximal growth rate, g C g C-1 day-1 (usually the calibrated vector)
    ks: np.ndarray            # half saturation constant, g C m-3
    death: np.ndarray         # turnover rate d, g C g C-1 day-1
    resp: np.ndarray          # respiration rate at Topt, g C g C-1 day-1
    faec: np.ndarray          # fraction of consumed C to faeces
    cn: np.ndarray            # C:N of the biomass
    rec: np.ndarray           # recalcitrance of the biomass (unused in the equations)
    mcn: np.ndarray           # sensitivity of gmax to the C:N of the food, pmC:N (Eqs. 21-22)
    mrec: np.ndarray          # sensitivity of gmax to litter recalcitrance, pmRec (Eqs. 23-24)
    tmin: np.ndarray          # temperature response (Eq. 17), degC
    topt: np.ndarray
    tmax: np.ndarray
    q10: np.ndarray           # Q10 of respiration
    # --- soil (KL_initSoil.txt)
    d: float                  # depth of the soil layer, m
    BD: float                 # bulk density, kg m-3 (not used by the equations)
    alpha: float              # van Genuchten parameters of the matrix
    n: float
    m: float
    Ksat: float
    pH: float
    litterCN0: float          # initial C:N of the litter pool
    SOMCN0: float             # initial C:N of SOM
    drainmax: float           # infiltration capacity of the layer below, mm day-1
    PVstruct: np.ndarray      # textural pore volume of the 5 classes (inaccessible, bact, micro, meso, macro), l m-3
    # --- engineers (KL_engineerParams.txt)
    ratioPVBeng: float        # burrow volume per engineer biomass, l g C-1
    fPVB: float               # fraction of burrow volume that is additional porosity
    tPVB: float               # daily turnover of burrows
    PVBmax: float             # maximal burrow volume, l m-3
    frag: float               # fragmentation: fraction of litter eaten by engineers that becomes SOM
    pfaec: np.ndarray         # sensitivity of faeces fraction to C:N (index 5 detritivores, 6 engineers)
    bioturbRate: float        # % of SOM mixed downwards per g C engineer per day
    moveRate: float           # % of litter moved downwards per g C engineer per day
    # --- run parameters (KL_runparams.txt)
    tStop: int                # number of simulated days
    initWater: float          # initial water content, % of pore volume
    Nmin0: float              # initial mineral N, g N m-3
    rrg: float                # daily root growth, g C m-3 day-1
    rootTO: float             # root turnover, day-1
    inputLit: float           # daily litter input, g C m-3 day-1
    CNlit: float              # C:N of the litter input
    recLit: float             # recalcitrance of the litter input, %
    CtoMyc: float             # plant C input to mycorrhiza, g C m-3 day-1
    NmyctoPlant: float        # fraction of plant N uptake delivered by mycorrhiza
    ee: float                 # effective / potential evapotranspiration
    # --- climate (KL_climateParams.txt, daily climate record)
    TEMP: np.ndarray          # monthly mean temperature, degC
    SUNH: np.ndarray          # monthly sunshine hours
    climatefile: str          # daily record: DOY year month DOM P Temperature (one line per simulated day)
    # --- initial pools (KL_initC_Pools.txt), g C m-3: 9 groups, litter, SOM, roots, CO2
    B0: np.ndarray
    # --- density dependent turnover (KL_densityParams.txt), death = d*B*(B/Bref)**(beta-1)
    density_dependent_turnover: bool = True
    beta: np.ndarray = None   # exponent per group (1 = first order)
    bref: np.ndarray = None   # reference biomass per group, g C m-3 (<= 0: initial biomass)

    @classmethod
    def from_files(cls, directory='.', climatefile='PrecipKMIBrass.txt', density_dependent_turnover=True,
                   tStop=None):
        """Read the KL_*.txt parameter files from `directory`. `tStop` overrides the run length of
        KL_runparams.txt; the climate file must have at least tStop lines after its header."""
        P = import_pools('KL_FaunalParams', directory)
        (gmax, ks, death, resp, faec, cn, rec, mcn, mrec, tmin, topt, tmax, q10) = P
        soil = import_pools('KL_initSoil', directory)
        eng = import_pools('KL_engineerParams', directory)
        run = import_pools('KL_runparams', directory)
        TEMP, SUNH = import_pools('KL_climateParams', directory)
        B0 = import_pools('KL_initC_Pools', directory)
        pfaec = np.zeros(7)
        pfaec[5], pfaec[6] = eng[5], eng[6]
        beta = bref = None
        if density_dependent_turnover:
            beta, bref = import_pools('KL_densityParams', directory)
        cfg = cls(gmax=gmax, ks=ks, death=death, resp=resp, faec=faec, cn=cn, rec=rec, mcn=mcn, mrec=mrec,
                  tmin=tmin, topt=topt, tmax=tmax, q10=q10,
                  d=soil[0], BD=soil[1], alpha=soil[2], n=soil[3], m=soil[4], Ksat=soil[5], pH=soil[6],
                  litterCN0=soil[7], SOMCN0=soil[8], drainmax=soil[9], PVstruct=np.array(soil[10:15]),
                  ratioPVBeng=eng[0], fPVB=eng[1], tPVB=eng[2], PVBmax=eng[3], frag=eng[4], pfaec=pfaec,
                  bioturbRate=eng[7], moveRate=eng[8],
                  tStop=int(run[0] if tStop is None else tStop), initWater=run[1], Nmin0=run[2], rrg=run[3],
                  rootTO=run[4], inputLit=run[5], CNlit=run[6], recLit=run[7], CtoMyc=run[8],
                  NmyctoPlant=run[9], ee=run[10],
                  TEMP=TEMP, SUNH=SUNH, climatefile=os.path.join(directory, climatefile) if directory != '.' else climatefile,
                  B0=np.array(B0[:NPOOLS], dtype=float),
                  density_dependent_turnover=density_dependent_turnover, beta=beta, bref=bref)
        cfg.validate()
        return cfg

    def replace(self, **changes):
        """A copy of the configuration with some fields changed (validated)."""
        cfg = dataclasses.replace(self, **changes)
        cfg.validate()
        return cfg

    def validate(self):
        for name in ('gmax', 'ks', 'death', 'resp', 'faec', 'cn', 'mcn', 'mrec', 'tmin', 'topt', 'tmax', 'q10'):
            if np.shape(getattr(self, name)) != (NGROUPS,):
                raise ValueError('%s must have %d values (one per functional group)' % (name, NGROUPS))
        if np.shape(self.PVstruct) != (5,) or np.any(self.PVstruct <= 0):
            raise ValueError('PVstruct must hold 5 positive pore volumes')
        if np.shape(self.B0) != (NPOOLS,):
            raise ValueError('B0 must hold %d initial pools' % NPOOLS)
        if self.tStop < 1:
            raise ValueError('tStop must be at least 1 day')
        if self.density_dependent_turnover:
            if np.shape(self.beta) != (NGROUPS,) or np.shape(self.bref) != (NGROUPS,):
                raise ValueError('beta and bref must have %d values' % NGROUPS)
            if np.any(self.beta < 1):
                raise ValueError('beta must be >= 1')

    @property
    def bref_effective(self):
        """Reference biomass for the density dependent turnover: bref, or the initial biomass where bref <= 0."""
        if not self.density_dependent_turnover:
            return None
        return np.where(self.bref > 0, self.bref, self.B0[:NGROUPS])

    @property
    def beta_effective(self):
        return self.beta if self.density_dependent_turnover else None


# ---------------------------------------------------------------------------------------------
@dataclass
class KeylinkResults:
    """Daily output of one run. Rows are days; all pools in g C m-3, water and pore volumes in l m-3,
    nitrogen in g N m-3."""
    gmax: np.ndarray                 # the gmax vector the run was made with
    pools: np.ndarray                # (tStop, 22): 13 pools (see BACT .. CO2) and the daily fluxes 13-20
    soil_water: np.ndarray           # (tStop,) total soil water
    porosity: np.ndarray             # (tStop,) total pore volume
    pores: np.ndarray                # (tStop, 5) pore volume per size class
    Nmin: np.ndarray                 # (tStop,) mineral N at the end of the day
    Ntot: np.ndarray                 # (tStop,) total N of the layer at the end of the day
    SOMCN: np.ndarray                # (tStop,) C:N of SOM
    litterCN: np.ndarray             # (tStop,) C:N of the litter pool
    mineralisation: np.ndarray       # (tStop,) net N mineralisation by bacteria of the day
    B0: np.ndarray                   # initial pools
    Nmin0: float = 0.0               # initial mineral N
    Ntot0: float = 0.0               # initial total N
    errors: list = field(default_factory=list)   # messages about numerical problems (e.g. negative SOM N)

    @property
    def ndays(self):
        return self.pools.shape[0]


# ---------------------------------------------------------------------------------------------
class KeylinkModel:
    """The KEYLINK soil model for one configuration. `run()` integrates the daily loop of
    Flores et al. (2021) ('Calculation order') and returns a KeylinkResults; it can be called
    repeatedly (e.g. in a calibration), every call starts from the initial state of the config."""

    def __init__(self, config):
        self.cfg = cfg = config
        cfg.validate()
        # heat index and exponent of the Thornthwaite equation (from the monthly mean temperatures)
        HI = float(np.sum((cfg.TEMP/5.0)**1.514))
        self.HI = HI
        self.alfa = 0.000000675*HI**3 - 0.0000771*HI**2 + 0.01792*HI + 0.49239
        # initial nitrogen: biota (fixed C:N) + mineral + SOM + litter
        B0 = cfg.B0
        self.Ntot0 = float(np.sum(B0[:NGROUPS]/cfg.cn) + cfg.Nmin0 + B0[SOM]/cfg.SOMCN0 + B0[LITTER]/cfg.litterCN0)
        self.iniSOM = float(B0[SOM])                      # SOM used for the inaccessible-pore fraction
        self.PW0 = cfg.initWater/100.0*cfg.PVstruct        # initial water per pore class, l m-3

    # -----------------------------------------------------------------------------------------
    def _read_climate(self):
        """Daily (year, month index, precipitation, temperature) for the tStop days of the run."""
        cfg = self.cfg
        data = np.loadtxt(cfg.climatefile, skiprows=1, encoding='latin1')
        if data.ndim == 1:
            data = data[None, :]
        if data.shape[0] < cfg.tStop:
            raise ValueError('climate file %s has %d days, run needs %d' % (cfg.climatefile, data.shape[0], cfg.tStop))
        data = data[:cfg.tStop]
        return data[:, 1].astype(int), data[:, 2].astype(int) - 1, data[:, 4], data[:, 5]

    # -----------------------------------------------------------------------------------------
    def run(self, gmax=None):
        """Run the model for cfg.tStop days with the given gmax vector (default cfg.gmax)."""
        cfg = self.cfg
        GMAX = np.array(cfg.gmax if gmax is None else gmax, dtype=float)
        if GMAX.shape != (NGROUPS,):
            raise ValueError('gmax must have %d values' % NGROUPS)
        T = cfg.tStop
        years, months, precips, temps = self._read_climate()

        # per-run state
        B = np.zeros(NSTATE); B[:NPOOLS] = cfg.B0
        PW = self.PW0.copy()
        pv = cfg.PVstruct.copy()
        Nmin, Ntot = cfg.Nmin0, self.Ntot0
        litterCN, SOMCN = cfg.litterCN0, cfg.SOMCN0
        modt = np.zeros(NGROUPS); rRESP = np.zeros(NGROUPS)
        BETA, BREF = cfg.beta_effective, cfg.bref_effective

        res = KeylinkResults(gmax=GMAX.copy(), pools=np.zeros((T, NSTATE)), soil_water=np.zeros(T),
                             porosity=np.zeros(T), pores=np.zeros((T, 5)), Nmin=np.zeros(T), Ntot=np.zeros(T),
                             SOMCN=np.zeros(T), litterCN=np.zeros(T), mineralisation=np.zeros(T),
                             B0=cfg.B0.copy(), Nmin0=cfg.Nmin0, Ntot0=self.Ntot0)

        for i in range(T):
            # a) pore size distribution (aggregation, burrows) ------------------------------------
            ag = mf.calcAg(B[FUNGI], B[MYC], B[SOM])
            pv = mf.calcPVD(cfg.PVstruct, pv, ag, cfg.ratioPVBeng, cfg.fPVB, cfg.tPVB, cfg.PVBmax, cfg.d, B)
            res.pores[i] = pv

            # climate of the day ------------------------------------------------------------------
            nd = DAYS_IN_MONTH.copy()
            if _is_leap_year(years[i]):
                nd[1] += 1
            precip, temp = precips[i], temps[i]
            for j in range(NGROUPS):
                modt[j] = mf.calcmodt(temp, cfg.topt[j], cfg.tmin[j], cfg.tmax[j])
                rRESP[j] = mf.calcresp(temp, cfg.topt[j], cfg.resp[j], cfg.q10[j])
            pet = mf.PET(temp, nd[months[i]], cfg.SUNH[months[i]], self.HI, self.alfa)

            # b) water ----------------------------------------------------------------------------
            SatW = sum(pv[0:4])/(cfg.d - pv[4]/sum(pv[0:5]))   # saturated water content (excluding macropores)
            PW, drainage, runoff = mf.calcPW(pv, precip, PW, cfg.drainmax, cfg.d, pv[0]/cfg.d, SatW,
                                             cfg.alpha, cfg.n, cfg.m, cfg.Ksat)
            res.soil_water[i] = sum(PW)
            res.porosity[i] = sum(pv)

            # d) availability of each pool to each consumer --------------------------------------
            avail = mf.calcAvail(pv, PW, self.iniSOM)

            # litter input (plant litter and root turnover) and its nitrogen ----------------------
            lit_in = cfg.inputLit + mf.rootTurnover(B[ROOTS], cfg.rootTO)
            litterCN = (B[LITTER] + lit_in)/(B[LITTER]/litterCN + lit_in/cfg.CNlit)
            B[LITTER] += lit_in
            Ntot += lit_in/cfg.CNlit
            B[ROOTS] += mf.rootgrowth(cfg.rrg) - mf.rootTurnover(B[ROOTS], cfg.rootTO)
            # (the plant C input to mycorrhiza, CtoMyc, is added inside the ODE)

            # plant N uptake: via mycorrhiza first (mineral N first, the rest from SOM), then direct
            CNsoil = (B[LITTER] + B[SOM])/(B[LITTER]/litterCN + B[SOM]/SOMCN)
            NviaMyc = mf.mycNtoPlant(cfg.NmyctoPlant, cfg.inputLit, litterCN, cfg.CtoMyc, CNsoil)
            Ntot -= NviaMyc
            Nmin = Nmin - NviaMyc if NviaMyc < Nmin else 0.0
            upt = min(mf.plantNuptake(litterCN, cfg.inputLit, cfg.NmyctoPlant), Nmin)
            Ntot -= upt
            Nmin -= upt

            # e)-g) growth modifiers and the daily ODE of all C pools ------------------------------
            B[13:NSTATE] = 0   # flux outputs of this day
            day = odeint(mf.fCompSpecies, B, [i, i+1],
                         args=(avail, modt, GMAX, litterCN, SOMCN, mf, cfg.cn, cfg.mcn, cfg.mrec, cfg.pH,
                               cfg.recLit, cfg.faec, cfg.pfaec, rRESP, cfg.ks, cfg.death, cfg.CtoMyc, BETA, BREF))
            B = day[1, :].copy()
            res.pools[i] = B

            # evapotranspiration at the end of the day (otherwise the soil is always dry)
            PW = mf.wl(PW, cfg.ee*pet)

            # h) nitrogen: bacterial mineralisation, SOM C:N as the residual of the N budget ------
            mineralisation = (B[BACT_GROWTH_SOM]/SOMCN + B[BACT_GROWTH_LIT]/litterCN
                              - (B[BACT_GROWTH_SOM] + B[BACT_GROWTH_LIT] - B[BACT_RESP])/cfg.cn[BACT])
            Nmin = max(0.0, Nmin + mineralisation)   # deficit is taken from SOM (the residual below)
            Nfauna = np.sum(B[:NGROUPS]/cfg.cn)
            NSOM = Ntot - Nfauna - Nmin - B[LITTER]/litterCN
            if NSOM <= 0:
                res.errors.append('day %d: SOM nitrogen <= 0 (%.3f g N m-3)' % (i, NSOM))
                NSOM = abs(NSOM)
            SOMCN = B[SOM]/NSOM

            # i2) bioturbation: SOM and litter moved out of the layer by the engineers ------------
            SOMdown = mf.calcBioturb(B[ENG], cfg.bioturbRate, B[SOM])
            B[SOM] -= SOMdown
            Ntot -= SOMdown/SOMCN
            Litdown = mf.calcLittermove(B[ENG], cfg.moveRate, B[LITTER])
            B[LITTER] -= Litdown
            Ntot -= Litdown/litterCN

            res.Nmin[i], res.Ntot[i], res.SOMCN[i], res.litterCN[i], res.mineralisation[i] = \
                Nmin, Ntot, SOMCN, litterCN, mineralisation

            # fragmentation of litter eaten by engineers into SOM
            # NOTE: index 21 of the state is the spare (always 0) entry, LIT_EATEN_ENG is index 20;
            # kept as in the original code so this refactoring does not change results.
            B[LITTER] -= cfg.frag*res.pools[i, 21]
            B[SOM] += cfg.frag*res.pools[i, 21]

            for s in range(SOM+1):   # security: avoid negative biomasses
                if B[s] <= 0:
                    B[s] = 0.001
        return res


# ---------------------------------------------------------------------------------------------
def show_plot(results, gr=(4000, 2000, 500, 50, 0.1)):
    """Plot soil water / porosity and the C pools of a run in four biomass ranges (gr)."""
    import matplotlib.pyplot as plt
    soln, tStop = results.pools, results.ndays
    t = np.arange(0., tStop)

    plt.clf()
    plt.subplot(2, 1, 1)
    plt.plot(t, results.soil_water, label="SW")
    plt.plot(t, results.porosity, label="porosity")
    plt.ylabel('Water, l/m2')
    plt.ylim(0, 400)
    plt.legend(loc=(1.01, 0), shadow=True)
    plt.subplot(2, 1, 2)
    for p in range(NPOOLS):   # pools with maximum above gr[1]
        if max(soln[:, p]) > gr[1]:
            plt.plot(t, soln[:, p], label=POOL_NAMES[p], c=POOL_COLOURS[p])
    plt.xlabel('Time,days'); plt.ylabel('Biomass, gC/m3'); plt.ylim(0, gr[0])
    plt.legend(loc=(1.01, 0), shadow=True)
    plt.show()

    def _range_page(upper, lower):
        """pools with maximum or mean biomass between gr[lower] and gr[upper] (upper < lower index)"""
        for p in range(NPOOLS):
            maxvalue = max(soln[:, p]); meanvalue = sum(soln[:, p])/tStop
            if (gr[lower] < maxvalue < gr[upper]) or (gr[lower] < meanvalue < gr[upper]):
                plt.plot(t, soln[:, p], label=POOL_NAMES[p], c=POOL_COLOURS[p])
        plt.ylabel('Biomass, gC/m3'); plt.ylim(0, gr[upper])
        plt.legend(loc=(1.01, 0), shadow=True)

    plt.clf()
    plt.subplot(2, 1, 1); _range_page(1, 2)
    plt.subplot(2, 1, 2); _range_page(2, 3); plt.xlabel('Time,days')
    plt.show()

    plt.clf()
    plt.subplot(2, 1, 1); _range_page(3, 4)
    plt.subplot(2, 1, 2)
    te = tStop - int(round(tStop/5))
    for p in range(NPOOLS):   # mean of the last 20 % of days below gr[4] (assumed local extinction)
        if sum(soln[te:tStop, p])/(tStop-te) < gr[4]:
            plt.plot(t, soln[:, p], label=POOL_NAMES[p], c=POOL_COLOURS[p])
    plt.xlabel('Time,days'); plt.ylabel('Biomass, gC/m3'); plt.ylim(0, gr[4])
    plt.legend(loc=(1.01, 0), shadow=True)
    plt.show()
