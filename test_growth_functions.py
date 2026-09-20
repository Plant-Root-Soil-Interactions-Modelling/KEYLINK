"""
Regression tests for the daily C-pool derivative functions
    keylink_functions.fCompSpecies   (used when CompetingSpecies == True)
    keylink_core.fEight              (used when CompetingSpecies == False)

They guard against two bugs found in both functions:
  * the summed flux outputs SOMeaten / LITeaten silently dropped the fungal,
    mycorrhizal, detritivore and engineer terms (continuation lines were not
    enclosed in parentheses, so only the bacterial term was assigned);
  * mycorrhiza grew with the fungal gmax / C:N / recalcitrance parameters
    (GMAX[1], CN[1], MCN[1], MREC[1]); its own parameters (index 2) were unused.

Run with:  python -m unittest test_growth_functions   (or pytest)
"""
import os
import unittest
import numpy as np

os.environ.setdefault('MPLBACKEND', 'Agg')  # keylink_core imports pyplot
import keylink_functions as mf
import keylink_core as core

# pool indices: 0=bact, 1=fungi, 2=myc, 3=bvores, 4=fvores, 5=sap (detritivores),
# 6=eng, 7=hvores, 8=pred, 9=litter, 10=SOM, 11=roots, 12=CO2
# output indices: 13=bactResp, 14=funResp, 15=EMresp, 16=bactGrowthSOM,
# 17=bactGrowthLit, 18=SOMeaten, 19=LITeaten, 20=LITeatenEng, 21=0
NPOOL = 22
SOMINI = 11470.


def base_inputs():
    """A plausible, fully 'active' state (all biomasses > 0, all availabilities > 0,
    no temperature limitation) built from the parameter files of the repository."""
    P = np.loadtxt('KL_FaunalParams.txt', encoding='latin1')
    GMAX, KS, DEATH, RESP, FAEC, CN, REC, MCN, MREC = P[:9]
    B = np.zeros(NPOOL)
    B[:12] = [15.1, 15.1, 160., 0.1, 0.8, 0.6, 0.2, 0.2, 0.4, 2680., SOMINI, 320.]
    SOMunavail = 45./325.*SOMINI
    avail = np.array([0.30, 0.25, 0.30, 0.20, 0.25, 0.30, 0.20, 0.20, 0.20, 0.20, 1.0, SOMunavail])
    modt = np.ones(9)
    kw = dict(avail=avail, modt=modt, GMAX=GMAX.copy(), litterCN=38.6, SOMCN=18., mf=mf,
              CN=CN.copy(), MCN=MCN.copy(), MREC=MREC.copy(), pH=3.9, recLit=40., FAEC=FAEC,
              pfaec=np.array([0, 0, 0, 0, 0, 0.3, 0.2]), rRESP=RESP.copy(), KS=KS, DEATH=DEATH,
              CtoMyc=0.54)
    return B, kw


def rhs(B, kw):
    return np.array(mf.fCompSpecies(B.copy(), 0., **kw), dtype=float)


def growth_of(B, kw, consumer, source_idx, avail_idx, gmax, SOM=False):
    """Independent re-computation of one Monod consumption term (Eq. 9)."""
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


class TestSummedConsumptionOutputs(unittest.TestCase):
    """SOMeaten (idx 18) and LITeaten (idx 19) must include ALL consumers."""

    def setUp(self):
        self.B, self.kw = base_inputs()
        self.out = rhs(self.B, self.kw)

    def expected_terms(self):
        B, kw = self.B, self.kw
        gEng = min(mf.calcgmaxEng(kw['GMAX'][6], kw['pH']), kw['GMAX'][6])
        fSAP = min(1, mf.calcFaec(kw['GMAX'][5], kw['FAEC'][5], kw['pfaec'][5], kw['litterCN'], kw['CN'][5], kw['rRESP'][5]))
        fEng = min(1, mf.calcFaec(gEng, kw['FAEC'][6], kw['pfaec'][6], kw['litterCN'], kw['CN'][6], kw['rRESP'][6]))
        som = {'bact': growth_of(B, kw, 0, 10, 0, gmax_on(kw, 0, 'SOM'), SOM=True),
               'fungi': growth_of(B, kw, 1, 10, 1, gmax_on(kw, 1, 'SOM'), SOM=True),
               'myc': growth_of(B, kw, 2, 10, 1, gmax_on(kw, 2, 'SOM'), SOM=True),
               'sap': growth_of(B, kw, 5, 10, 3, kw['GMAX'][5], SOM=True),
               'eng': growth_of(B, kw, 6, 10, 2, gEng, SOM=True)}
        lit = {'bact': growth_of(B, kw, 0, 9, 0, gmax_on(kw, 0, 'litter')),
               'fungi': growth_of(B, kw, 1, 9, 0, gmax_on(kw, 1, 'litter')),
               'myc': growth_of(B, kw, 2, 9, 0, gmax_on(kw, 2, 'litter')),
               'sap': (1+fSAP)*growth_of(B, kw, 5, 9, 0, kw['GMAX'][5]),
               'eng': (1+fEng)*growth_of(B, kw, 6, 9, 0, gEng)}
        return som, lit

    def test_all_terms_are_nonzero_in_this_state(self):
        som, lit = self.expected_terms()
        for d in (som, lit):
            for k, v in d.items():
                self.assertGreater(v, 0, k)

    def test_SOMeaten_is_sum_over_all_consumers(self):
        som, _ = self.expected_terms()
        self.assertAlmostEqual(self.out[18], sum(som.values()), places=10)
        # and not just the bacterial share (the old behaviour)
        self.assertGreater(self.out[18], som['bact']*1.5)

    def test_LITeaten_is_sum_over_all_consumers(self):
        _, lit = self.expected_terms()
        self.assertAlmostEqual(self.out[19], sum(lit.values()), places=10)
        self.assertGreater(self.out[19], lit['bact']*1.5)
        # engineer-only litter consumption output (idx 20) is the engineer term
        self.assertAlmostEqual(self.out[20], lit['eng'], places=10)

    def test_LITeaten_consistent_with_litter_derivative(self):
        """d(litter)/dt = -LITeaten + necromass of detritivores, engineers, herbivores, predators."""
        B, kw = self.B, self.kw
        deaths = sum(kw['DEATH'][i]*B[i] for i in (5, 6, 7, 8))
        self.assertAlmostEqual(self.out[9], -self.out[19] + deaths, places=10)

    def test_SOMeaten_consistent_with_SOM_derivative(self):
        """d(SOM)/dt = exudation - SOMeaten + faeces (pred, sap, eng on litter, hvores)
        + necromass of bacteria, fungi, myc, bacterivores, fungivores."""
        B, kw = self.B, self.kw
        _, lit = self.expected_terms()
        gEng = min(mf.calcgmaxEng(kw['GMAX'][6], kw['pH']), kw['GMAX'][6])
        fSAP = min(1, mf.calcFaec(kw['GMAX'][5], kw['FAEC'][5], kw['pfaec'][5], kw['litterCN'], kw['CN'][5], kw['rRESP'][5]))
        fEng = min(1, mf.calcFaec(gEng, kw['FAEC'][6], kw['pfaec'][6], kw['litterCN'], kw['CN'][6], kw['rRESP'][6]))
        faeces = (kw['FAEC'][8]*sum(growth_of(B, kw, 8, prey, av, kw['GMAX'][8])
                                    for prey, av in ((3, 7), (4, 6), (5, 9), (6, 10), (7, 8)))
                  + fSAP*growth_of(B, kw, 5, 9, 0, kw['GMAX'][5])
                  + fEng*growth_of(B, kw, 6, 9, 0, gEng)
                  + kw['FAEC'][7]*kw['modt'][7]*mf.calcgrowth(B[7], B[11], 1, kw['GMAX'][7], kw['KS'][7]))
        deaths = sum(kw['DEATH'][i]*B[i] for i in range(5))
        self.assertAlmostEqual(self.out[10], mf.exudation() - self.out[18] + faeces + deaths, places=10)


class TestMycorrhizaParameters(unittest.TestCase):
    """Mycorrhiza must use GMAX[2], CN[2], MCN[2], MREC[2] - not the fungal values."""

    def setUp(self):
        self.B, self.kw = base_inputs()
        self.ref = rhs(self.B, self.kw)

    def _with(self, **changes):
        kw = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in self.kw.items()}
        for key, (idx, val) in changes.items():
            kw[key][idx] = val
        return rhs(self.B, kw), kw

    def test_myc_growth_responds_to_its_own_gmax(self):
        out, _ = self._with(GMAX=(2, self.kw['GMAX'][2]*2))
        self.assertGreater(out[2], self.ref[2])
        # nothing but myc, litter, SOM and the summed outputs may change
        unchanged = [i for i in range(13) if i not in (2, 9, 10)]
        np.testing.assert_allclose(out[unchanged], self.ref[unchanged])
        # what myc gains is exactly what litter and SOM lose
        self.assertAlmostEqual(out[2]-self.ref[2], -(out[9]-self.ref[9]) - (out[10]-self.ref[10]), places=10)

    def test_myc_growth_does_not_respond_to_fungal_gmax(self):
        out, _ = self._with(GMAX=(1, self.kw['GMAX'][1]*2))
        self.assertAlmostEqual(out[2], self.ref[2], places=12)
        self.assertGreater(out[1], self.ref[1])

    def test_myc_uses_own_CN_and_recalcitrance_sensitivity(self):
        for key in ('MCN', 'MREC', 'CN'):
            own, _ = self._with(**{key: (2, self.kw[key][2]*0.5)})
            self.assertNotAlmostEqual(own[2], self.ref[2], places=8, msg=key + '[2] has no effect on myc')
            fung, _ = self._with(**{key: (1, self.kw[key][1]*0.5)})
            self.assertAlmostEqual(fung[2], self.ref[2], places=12, msg=key + '[1] wrongly affects myc')

    def test_myc_gmax_zero_leaves_only_input_minus_losses(self):
        out, kw = self._with(GMAX=(2, 0.0))
        B = self.B
        predation = kw['modt'][4]*mf.calcgrowth(B[4], B[2], kw['avail'][5], kw['GMAX'][4], kw['KS'][4])
        expected = kw['CtoMyc'] - kw['DEATH'][2]*B[2] - kw['rRESP'][2]*B[2] - predation
        self.assertAlmostEqual(out[2], expected, places=10)


class TestDetritivoresAndEngineers(unittest.TestCase):

    def setUp(self):
        self.B, self.kw = base_inputs()

    def test_engineer_growth_follows_pH_dependent_gmax(self):
        """Below pH 3 engineers cannot grow (Eq. 20): derivative = -death - resp - predation,
        and their litter/SOM consumption vanishes from the summed outputs."""
        kw = dict(self.kw); kw['pH'] = 2.5
        out = rhs(self.B, kw)
        B = self.B
        predation = (1+kw['FAEC'][8])*kw['modt'][8]*mf.calcgrowth(B[8], B[6], kw['avail'][10], kw['GMAX'][8], kw['KS'][8])
        self.assertAlmostEqual(out[6], -kw['DEATH'][6]*B[6] - kw['rRESP'][6]*B[6] - predation, places=10)
        self.assertEqual(out[20], 0.0)
        ref = rhs(self.B, self.kw)
        self.assertLess(out[19], ref[19])
        self.assertLess(out[18], ref[18])

    def test_detritivore_growth_uses_unmodified_gmax_and_own_availability(self):
        B, kw = self.B, self.kw
        out = rhs(B, kw)
        fSAP = min(1, mf.calcFaec(kw['GMAX'][5], kw['FAEC'][5], kw['pfaec'][5], kw['litterCN'], kw['CN'][5], kw['rRESP'][5]))
        g_lit = growth_of(B, kw, 5, 9, 0, kw['GMAX'][5])
        g_som = growth_of(B, kw, 5, 10, 3, kw['GMAX'][5], SOM=True)
        predation = (1+kw['FAEC'][8])*growth_of(B, kw, 8, 5, 9, kw['GMAX'][8])
        self.assertAlmostEqual(out[5], g_lit + g_som - predation - kw['DEATH'][5]*B[5] - kw['rRESP'][5]*B[5], places=10)
        # litter loses growth + faeces, SOM gets the faeces back
        self.assertGreater(fSAP, 0)
        kw0 = dict(kw); kw0 = {k: (v.copy() if isinstance(v, np.ndarray) else v) for k, v in kw.items()}
        kw0['GMAX'][5] = 0.0
        out0 = rhs(B, kw0)
        self.assertAlmostEqual(out0[9]-out[9], (1+fSAP)*g_lit, places=10)


class TestMassBalanceAndEquivalence(unittest.TestCase):

    def test_carbon_is_conserved(self):
        """Sum of all pool derivatives (biota, litter, SOM, roots, CO2) equals the
        external inputs (root exudation + plant C to mycorrhiza)."""
        B, kw = base_inputs()
        out = rhs(B, kw)
        self.assertAlmostEqual(out[:13].sum(), mf.exudation() + kw['CtoMyc'], places=9)

    def test_fEight_equals_fCompSpecies(self):
        """Both right-hand sides must give identical derivatives for identical inputs
        (fEight reads parameters from keylink_core globals, fCompSpecies gets them as arguments)."""
        B, kw = base_inputs()
        kw.update(CN=core.CN, MCN=core.MCN, MREC=core.MREC, pH=core.pH, recLit=core.recLit,
                  FAEC=core.FAEC, pfaec=core.pfaec, KS=core.KS, DEATH=core.DEATH, CtoMyc=core.CtoMyc)
        core.rRESP[:] = kw['rRESP']
        a = rhs(B, kw)
        b = np.array(core.fEight(B.copy(), 0., kw['avail'], kw['modt'], kw['GMAX'], kw['litterCN'], kw['SOMCN']), dtype=float)
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)
        self.assertGreater(b[18], 0); self.assertGreater(b[19], 0)

    def test_fEight_myc_uses_own_gmax(self):
        B, kw = base_inputs()
        core.rRESP[:] = kw['rRESP']
        def f(g):
            return np.array(core.fEight(B.copy(), 0., kw['avail'], kw['modt'], g, kw['litterCN'], kw['SOMCN']), dtype=float)
        g = kw['GMAX'].copy(); ref = f(g)
        g2 = g.copy(); g2[2] *= 2; self.assertGreater(f(g2)[2], ref[2])
        g1 = g.copy(); g1[1] *= 2; self.assertAlmostEqual(f(g1)[2], ref[2], places=12)


if __name__ == '__main__':
    unittest.main()
