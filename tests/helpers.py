"""
Shared helpers for the KEYLINK tests: a short climate fixture, a model runner that starts from a
clean module state, and a reference state for the daily derivative functions.
"""
import os
import numpy as np

from tests import MODEL_DIR
import keylink_functions as mf
import keylink_core as core

# two years (1999-2000, 731 days) of the Brasschaat record, instead of the full 10-year file
FIXTURE_CLIMATE = os.path.join(MODEL_DIR, 'tests', 'fixtures', 'climate_1999_2000.txt')
FIXTURE_DAYS = 731

# pool indices: 0=bact, 1=fungi, 2=myc, 3=bvores, 4=fvores, 5=sap (detritivores), 6=eng,
# 7=hvores, 8=pred, 9=litter, 10=SOM, 11=roots, 12=CO2; flux outputs: 13=bactResp, 14=funResp,
# 15=EMresp, 16=bactGrowthSOM, 17=bactGrowthLit, 18=SOMeaten, 19=LITeaten, 20=LITeatenEng, 21=0
NPOOL = 22
POOLS = slice(0, 13)
SOMINI = 11470.


def make_config(ndays=FIXTURE_DAYS, density_dependent=True, **changes):
    """KeylinkConfig from the parameter files of the model directory, on the climate fixture."""
    if ndays > FIXTURE_DAYS:
        raise ValueError('fixture climate file has only %d days' % FIXTURE_DAYS)
    cfg = core.KeylinkConfig.from_files('.', climatefile=FIXTURE_CLIMATE,
                                        density_dependent_turnover=density_dependent, tStop=ndays)
    return cfg.replace(**changes) if changes else cfg


def run_model(ndays=FIXTURE_DAYS, density_dependent=True, gmax=None, **changes):
    """Run KeylinkModel on the fixture climate and return the KeylinkResults
    (results.B0 = initial pools, results.pools = daily output)."""
    cfg = make_config(ndays, density_dependent, **changes)
    return core.KeylinkModel(cfg).run(cfg.gmax if gmax is None else gmax)


def base_inputs():
    """A plausible, fully 'active' state (all biomasses > 0, all availabilities > 0, no temperature
    limitation) built from the parameter files of the repository, for calling the derivative
    functions directly."""
    P = np.loadtxt('KL_FaunalParams.txt', encoding='latin1')
    GMAX, KS, DEATH, RESP, FAEC, CN, REC, MCN, MREC = P[:9]
    B = np.zeros(NPOOL)
    B[:12] = [15.1, 15.1, 160., 0.1, 0.8, 0.6, 0.2, 0.2, 0.4, 2680., SOMINI, 320.]
    SOMunavail = 45./325.*SOMINI
    avail = np.array([0.30, 0.25, 0.30, 0.20, 0.25, 0.30, 0.20, 0.20, 0.20, 0.20, 1.0, SOMunavail])
    kw = dict(avail=avail, modt=np.ones(9), GMAX=GMAX.copy(), litterCN=38.6, SOMCN=18., mf=mf,
              CN=CN.copy(), MCN=MCN.copy(), MREC=MREC.copy(), pH=3.9, recLit=40., FAEC=FAEC,
              pfaec=np.array([0, 0, 0, 0, 0, 0.3, 0.2]), rRESP=RESP.copy(), KS=KS, DEATH=DEATH,
              CtoMyc=0.54)
    return B, kw


def rhs(B, kw):
    """Daily derivatives from fCompSpecies as an array."""
    return np.array(mf.fCompSpecies(B.copy(), 0., **kw), dtype=float)


def growth_of(B, kw, consumer, source_idx, avail_idx, gmax, SOM=False):
    """Independent re-computation of one Monod consumption term (Flores et al. 2021, Eq. 9)."""
    S = B[source_idx] - kw['avail'][11] if SOM else B[source_idx]
    return kw['modt'][consumer]*mf.calcgrowth(B[consumer], S, kw['avail'][avail_idx], gmax, kw['KS'][consumer])


def gmax_on(kw, group, source):
    """gmax of bacteria(0)/fungi(1)/mycorrhiza(2) on 'litter' or 'SOM' incl. C:N/rec/pH modifiers."""
    pH_id = 1 if group == 0 else 2
    if source == 'litter':
        return mf.calcgmaxmod(kw['CN'][group], kw['litterCN'], kw['MCN'][group], kw['recLit'],
                              kw['MREC'][group], kw['pH'], pH_id)*kw['GMAX'][group]
    return mf.calcgmaxmod(kw['CN'][group], kw['SOMCN'], kw['MCN'][group], 0.0,
                          kw['MREC'][group], kw['pH'], pH_id)*kw['GMAX'][group]
