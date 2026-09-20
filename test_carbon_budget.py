"""
Whole-model carbon budget test for keylink_core.KeylinkModel.

The carbon held in the 13 modelled pools (9 biota, litter, SOM, roots and the
cumulative CO2 pool) may only change through fluxes across the model boundary:

  inputs  : daily litter fall (inputLit), fine-root growth (rrg), plant C to
            mycorrhiza (CtoMyc, added inside the ODE), root exudation to SOM
            (mf.exudation(), inside the ODE), and the 'security' clamp that
            lifts non-positive pools to 0.001 g C m-3
  outputs : SOM bioturbated out of the layer (calcBioturb) and litter moved
            down by engineers (calcLittermove)

Everything else (feeding, faeces, necromass, respiration -> CO2, root turnover
-> litter, fragmentation) is internal and must cancel in the sum. The test
rebuilds the expected budget day by day from the model output and checks that
the pools close to solver precision. It fails if CtoMyc is added twice (once in
the ODE and once more in the daily loop), which was the case before.

Run with:  python -m unittest test_carbon_budget   (or pytest)
"""
import os
import importlib
import unittest
import numpy as np

os.environ.setdefault('MPLBACKEND', 'Agg')
import keylink_functions as mf
import keylink_core as core

NDAYS = 200          # long enough to include wet/dry periods and pool clamps
POOLS = slice(0, 13)  # 9 biota + litter + SOM + roots + CO2


def run_model(ndays):
    """Fresh import (initial state is mutated by a run), shortened run length."""
    importlib.reload(core)
    core.tStop = ndays
    core.psoln = np.zeros((ndays, core.NrGroups+14))
    core.PVt = np.zeros(ndays)
    core.PWt = np.zeros(ndays)
    gmax = core.import_pools('KL_FaunalParams')[0]
    B0 = core.Bini.copy()
    C, PW, PV, pores = core.KeylinkModel(gmax.copy())
    return B0, C.copy()


def between_days(B_after_ode):
    """Replays the explicit (non-ODE) part of the daily loop that follows the ODE
    of day i and precedes the ODE of day i+1. Returns the state handed to the
    next ODE call plus the boundary fluxes involved."""
    B = B_after_ode.copy()
    somdown = mf.calcBioturb(B[6], core.bioturbRate, B[10]); B[10] -= somdown
    litdown = mf.calcLittermove(B[6], core.moveRate, B[9]);  B[9] -= litdown
    # fragmentation (litter -> SOM) is internal: no effect on the sum
    clamp = 0.0
    for s in range(11):
        if B[s] <= 0:
            clamp += 0.001 - B[s]
            B[s] = 0.001
    return start_of_day(B), somdown, litdown, clamp


def start_of_day(B):
    """Litter input, root turnover (internal) and root growth applied before the ODE."""
    B = B.copy()
    rt = mf.rootTurnover(B[11], core.rootTO)
    B[9] += core.inputLit + rt
    B[11] += core.rrg - rt
    return B


class TestCarbonBudget(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.B0, cls.C = run_model(NDAYS)

    def test_run_is_meaningful(self):
        """Sanity: pools are finite and positive, CO2 accumulates, engineers are
        present so that the bioturbation/litter-moving outputs are non-zero."""
        C = self.C
        self.assertTrue(np.all(np.isfinite(C[:, POOLS])))
        self.assertTrue(np.all(C[:, :12] > 0))
        self.assertGreater(C[-1, 12], 0)
        self.assertGreater(mf.calcBioturb(C[:, 6], core.bioturbRate, C[:, 10]).sum(), 0)

    def test_daily_ode_conserves_carbon(self):
        """Over each ODE call, the sum of the 13 pools increases by exactly
        exudation + CtoMyc (the only inputs inside the ODE)."""
        C = self.C
        pre = start_of_day(self.B0)
        residuals = []
        for i in range(NDAYS):
            d_sum = C[i, POOLS].sum() - pre[POOLS].sum()
            residuals.append(d_sum - (mf.exudation() + core.CtoMyc))
            pre, _, _, _ = between_days(C[i])
        residuals = np.abs(residuals)
        self.assertLess(residuals.max(), 1e-8, "max daily residual %.3e g C m-3" % residuals.max())

    def test_whole_run_budget_closes(self):
        """sum(pools_end) == sum(pools_start) + inputs - outputs over the whole run.
        Before the CtoMyc fix this failed by exactly CtoMyc * NDAYS."""
        C = self.C
        inputs = outputs = 0.0
        pre = start_of_day(self.B0)
        inputs += core.inputLit + core.rrg           # day 0 explicit inputs
        for i in range(NDAYS):
            inputs += mf.exudation() + core.CtoMyc   # inputs inside the ODE of day i
            if i < NDAYS-1:
                pre, somdown, litdown, clamp = between_days(C[i])
                outputs += somdown + litdown
                inputs += clamp + core.inputLit + core.rrg
        expected_end = self.B0[POOLS].sum() + inputs - outputs
        closure = C[-1, POOLS].sum() - expected_end
        self.assertAlmostEqual(closure, 0.0, places=6,
                               msg="budget not closed: %.6f g C m-3 (CtoMyc*NDAYS = %.2f)" % (closure, core.CtoMyc*NDAYS))


if __name__ == '__main__':
    unittest.main()
