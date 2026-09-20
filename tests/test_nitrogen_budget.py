"""
Nitrogen bookkeeping test for keylink_core.KeylinkModel.

Design of the N bookkeeping in KEYLINK (one layer):
  * biota have a fixed C:N, so their N is B/CN;
  * litter C:N changes only by mixing with the daily litter input (CNlit);
  * mineral N (Nmin) changes by plant uptake (direct and via mycorrhiza) and by the
    net mineralisation / immobilisation of the bacteria;
  * SOM N is the residual of the total N budget, SOMCN = SOM / NSOM;
  * total N (Ntot) changes only by the N fluxes crossing the model boundary:
      + litter fall and root turnover (at CNlit)
      - N delivered to the plant by mycorrhiza, - direct plant N uptake
      - SOM bioturbated out of the layer, - litter moved out by engineers.

The tests replay these rules from the daily model output and the N diagnostics
(Nmint, Ntott, SOMCNt, litterCNt, mineralisationt) and check that
  (1) the initial total N really contains the SOM and litter N,
  (2) on every day   Ntot == N_biota + Nmin + N_litter + N_SOM   (closure),
  (3) Ntot changes only by the boundary fluxes listed above,
  (4) the bacterial mineralisation uses the DAILY growth/respiration fluxes with the
      C:N of the right source (SOM growth / SOM C:N, litter growth / litter C:N),
  (5) the flux outputs (columns 13-20) are daily values, not cumulative sums,
  (6) the SOM C:N stays physically plausible (no sign flips / 'error NSOM' hacks).

Run with:  python run_tests.py   (all tests) from the model directory.
"""
import unittest
import numpy as np

from tests.helpers import mf, core, run_model

NDAYS = 365


class TestNitrogenBudget(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.B0, cls.C = run_model(NDAYS)
        cls.CN = core.CN

    # ------------------------------------------------------------------ helpers
    def state_before_ode(self, i):
        """Pools handed to the ODE of day i and the explicit N changes made before it.
        Returns (B_pre, litterCN_i, SOMCN_used, Nmin_pre, Ntot_pre, dNtot_explicit)."""
        C = self.C
        if i == 0:
            B = self.B0.copy(); litterCN = core.litterCNini; SOMCN = core.SOMCNini
            Nmin = core.Nminini; Ntot = core.Ntotini
        else:
            B = C[i-1].copy()
            B[10] -= mf.calcBioturb(B[6], core.bioturbRate, B[10])
            B[9] -= mf.calcLittermove(B[6], core.moveRate, B[9])
            for s in range(11):
                if B[s] <= 0:
                    B[s] = 0.001
            litterCN = core.litterCNt[i-1]; SOMCN = core.SOMCNt[i-1]
            Nmin = core.Nmint[i-1]; Ntot = core.Ntott[i-1]
        lit_in = core.inputLit + mf.rootTurnover(B[11], core.rootTO)
        litterCN = (B[9]+lit_in)/(B[9]/litterCN + lit_in/core.CNlit)
        B[9] += lit_in
        Ntot += lit_in/core.CNlit
        B[11] += core.rrg - mf.rootTurnover(B[11], core.rootTO)
        CNsoil = (B[9]+B[10])/(B[9]/litterCN + B[10]/SOMCN)
        NviaMyc = mf.mycNtoPlant(core.NmyctoPlant, core.inputLit, litterCN, core.CtoMyc, CNsoil)
        Ntot -= NviaMyc
        Nmin = Nmin - NviaMyc if NviaMyc < Nmin else 0.0
        upt = min(mf.plantNuptake(litterCN, core.inputLit, core.NmyctoPlant), Nmin)
        Ntot -= upt; Nmin -= upt
        return B, litterCN, SOMCN, Nmin, Ntot

    def end_of_day_pools(self, i):
        B = self.C[i].copy()
        somdown = mf.calcBioturb(B[6], core.bioturbRate, B[10]); B[10] -= somdown
        litdown = mf.calcLittermove(B[6], core.moveRate, B[9]);  B[9] -= litdown
        return B, somdown, litdown

    # -------------------------------------------------------------------- tests
    def test_initial_total_N_includes_SOM_and_litter(self):
        B0 = self.B0
        expected = (B0[:9]/self.CN).sum() + core.Nminini + B0[10]/core.SOMCNini + B0[9]/core.litterCNini
        self.assertAlmostEqual(core.Ntotini, expected, places=9)
        self.assertGreater(core.Ntotini, B0[10]/core.SOMCNini)   # SOM N (~637) is really in there
        self.assertEqual(core.iniSOM, B0[10])                      # and iniSOM is the SOM pool

    def test_daily_closure_total_equals_sum_of_pools(self):
        """Ntot == N_biota + Nmin + N_litter + N_SOM at the end of every day."""
        for i in range(NDAYS):
            B, _, _ = self.end_of_day_pools(i)
            Nsum = (B[:9]/self.CN).sum() + core.Nmint[i] + B[9]/core.litterCNt[i] + B[10]/core.SOMCNt[i]
            self.assertAlmostEqual(Nsum, core.Ntott[i], places=8, msg="day %d" % i)

    def test_total_N_changes_only_by_boundary_fluxes(self):
        """Replay: Ntot(i) = Ntot(i-1) + litter N in - N to plant (myc + direct) - N exported by engineers."""
        exported = 0.0
        for i in range(NDAYS):
            _, litterCN, SOMCN_used, _, Ntot_pre = self.state_before_ode(i)
            _, somdown, litdown = self.end_of_day_pools(i)
            expected = Ntot_pre - somdown/core.SOMCNt[i] - litdown/litterCN
            self.assertAlmostEqual(core.Ntott[i], expected, places=9, msg="day %d" % i)
        exported += somdown/core.SOMCNt[i] + litdown/litterCN
        # all four boundary routes are actually active in this run (so the replay is a real check)
        self.assertGreater(exported, 0)                                        # engineer export
        self.assertGreater(mf.plantNuptake(core.litterCNt[0], core.inputLit, core.NmyctoPlant), 0)
        self.assertGreater(core.inputLit/core.CNlit, 0)                        # litter N input
        self.assertNotAlmostEqual(core.Ntott[-1], core.Ntotini, places=3)      # Ntot does change

    def test_mycorrhizal_N_transfer_uses_plant_demand_from_litter_input(self):
        """mycNtoPlant must be evaluated with inputLit as plant-demand proxy (as plantNuptake is),
        not with CtoMyc: check the day-0 Ntot against both variants."""
        B = self.B0.copy()
        lit_in = core.inputLit + mf.rootTurnover(B[11], core.rootTO)
        litterCN = (B[9]+lit_in)/(B[9]/core.litterCNini + lit_in/core.CNlit)
        B[9] += lit_in
        CNsoil = (B[9]+B[10])/(B[9]/litterCN + B[10]/core.SOMCNini)
        right = mf.mycNtoPlant(core.NmyctoPlant, core.inputLit, litterCN, core.CtoMyc, CNsoil)
        wrong = mf.mycNtoPlant(core.NmyctoPlant, core.CtoMyc, litterCN, core.CtoMyc, CNsoil)
        self.assertNotAlmostEqual(right, wrong, places=6)
        _, somdown, litdown = self.end_of_day_pools(0)
        Nmin0 = core.Nminini - right
        upt = min(mf.plantNuptake(litterCN, core.inputLit, core.NmyctoPlant), Nmin0)
        self.assertAlmostEqual(core.Ntott[0], core.Ntotini + lit_in/core.CNlit - right - upt
                               - somdown/core.SOMCNt[0] - litdown/litterCN, places=9)

    def test_bacterial_mineralisation_uses_daily_fluxes_and_right_CN(self):
        """mineralisation = G_SOM/CN_SOM + G_lit/CN_lit - (G_SOM + G_lit - R)/CN_bact with the
        fluxes of that day (columns 16, 17, 13) and Nmin(i) = max(0, Nmin_pre + mineralisation)."""
        C = self.C
        n_checked = 0
        for i in range(NDAYS):
            _, litterCN, SOMCN_used, Nmin_pre, _ = self.state_before_ode(i)
            G_som, G_lit, R = C[i, 16], C[i, 17], C[i, 13]
            expected = G_som/SOMCN_used + G_lit/litterCN - (G_som + G_lit - R)/self.CN[0]
            self.assertAlmostEqual(core.mineralisationt[i], expected, places=10, msg="day %d" % i)
            self.assertAlmostEqual(core.Nmint[i], max(0.0, Nmin_pre + expected), places=10, msg="day %d" % i)
            if G_som > 1e-6 and G_lit > 1e-6:
                # with the C:N swapped the value must differ (guards the index fix)
                swapped = G_som/litterCN + G_lit/SOMCN_used - (G_som + G_lit - R)/self.CN[0]
                self.assertNotAlmostEqual(expected, swapped, places=8)
                n_checked += 1
        self.assertGreater(n_checked, 10)

    def test_flux_outputs_are_daily_not_cumulative(self):
        """Columns 13-20 hold the flux of the day: they must be bounded by what the pool can
        produce in one day and must not grow monotonically like a running sum."""
        C = self.C
        resp_rate_max = core.RESP[0]*core.Q10[0]**((core.T_MAX[0]-core.T_OPT[0])/10)
        for i in range(1, NDAYS):
            bact_max = max(C[i-1, 0], C[i, 0])*3    # generous bound on within-day biomass
            self.assertLessEqual(C[i, 13], resp_rate_max*bact_max + 1e-9, "day %d" % i)
        running = np.cumsum(C[:, 13])
        self.assertLess(C[-1, 13], 0.5*running[-1])
        self.assertFalse(np.all(np.diff(C[:, 13]) >= 0))

    def test_SOM_CN_stays_plausible(self):
        self.assertTrue(np.all(core.SOMCNt > 5) and np.all(core.SOMCNt < 60),
                        "SOM C:N out of range: %.2f .. %.2f" % (core.SOMCNt.min(), core.SOMCNt.max()))
        self.assertLess(abs(core.SOMCNt[0]-core.SOMCNini), 0.5)
        self.assertTrue(np.all(core.Nmint >= 0))
        self.assertTrue(np.all(core.litterCNt > 0))


if __name__ == '__main__':
    unittest.main()
