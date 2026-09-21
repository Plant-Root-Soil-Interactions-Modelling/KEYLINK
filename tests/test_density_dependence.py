"""
Tests for the density-dependent turnover add-on (calcDeath, KeylinkConfig.density_dependent_turnover)
and for the stability of the microbial pools in whole-model runs.

Background: the growth equation (Flores et al. 2021, Eq. 9) is exponential in the consumer
biomass and, with first-order turnover (Eq. 12), the only brakes are substrate depletion and
predation. This produced boom-bust dynamics (microbial peaks far above plausible values
followed by decay to ~1e-60 g C m-3). Flores et al. (2021) suggest density-dependent microbial
turnover (Georgiou et al. 2017) as an add-on; it is implemented as
    Dt = d * B * (B/Bref)**(beta-1)
with beta and Bref per functional group in KL_densityParams.txt (beta = 1 -> Eq. 12).

Run with:  python run_tests.py   (all tests) from the model directory.
"""
import unittest
import numpy as np

from tests.helpers import core, run_model, make_config, base_inputs, rhs, FIXTURE_DAYS, foodweb, plant

NDAYS = FIXTURE_DAYS   # 2 years


class TestCalcDeath(unittest.TestCase):

    def test_first_order_when_beta_is_one(self):
        for B in (0.0, 0.3, 15.1, 5000.):
            self.assertEqual(foodweb.calcDeath(B, 0.05, 1.0, 15.1), 0.05*B)
            self.assertEqual(foodweb.calcDeath(B, 0.05), 0.05*B)

    def test_equals_d_at_reference_biomass_and_scales_with_density(self):
        d, Bref = 0.05, 15.1
        for beta in (1.5, 2.0):
            self.assertAlmostEqual(foodweb.calcDeath(Bref, d, beta, Bref), d*Bref)
            # specific turnover at 10x the reference biomass is 10**(beta-1) times d
            self.assertAlmostEqual(foodweb.calcDeath(10*Bref, d, beta, Bref)/(10*Bref), d*10**(beta-1))
            # and at 0.1x it is reduced accordingly (slower decline of small populations)
            self.assertAlmostEqual(foodweb.calcDeath(0.1*Bref, d, beta, Bref)/(0.1*Bref), d*0.1**(beta-1))

    def test_robust_for_non_positive_biomass(self):
        self.assertEqual(foodweb.calcDeath(0.0, 0.05, 2.0, 15.1), 0.0)
        self.assertTrue(np.isfinite(foodweb.calcDeath(-1e-9, 0.05, 1.7, 15.1)))


class TestInDerivativeFunctions(unittest.TestCase):

    def setUp(self):
        self.B, self.kw = base_inputs()
        self.beta = np.full(9, 2.0)
        self.bref = self.B[:9].copy()

    def test_beta_one_reproduces_first_order_model(self):
        a = rhs(self.B, self.kw)
        b = rhs(self.B, dict(self.kw, BETA=np.ones(9), BREF=self.bref))
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)

    def test_at_reference_biomass_identical_to_first_order(self):
        """With B == Bref the density-dependent turnover equals d*B: same derivatives."""
        a = rhs(self.B, self.kw)
        b = rhs(self.B, dict(self.kw, BETA=self.beta, BREF=self.bref))
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-12)

    def test_turnover_increases_above_reference_and_carbon_is_conserved(self):
        B = self.B.copy(); B[2] *= 10      # mycorrhiza at 10x its reference biomass
        kw_dd = dict(self.kw, BETA=self.beta, BREF=self.bref)
        a = rhs(B, self.kw); b = rhs(B, kw_dd)
        extra_death = foodweb.calcDeath(B[2], self.kw['DEATH'][2], 2.0, self.bref[2]) - self.kw['DEATH'][2]*B[2]
        self.assertGreater(extra_death, 0)
        self.assertAlmostEqual(a[2]-b[2], extra_death, places=10)      # myc loses it
        self.assertAlmostEqual(b[10]-a[10], extra_death, places=10)    # SOM (microbial necromass) gains it
        self.assertAlmostEqual(b[:13].sum(), plant.exudation()+kw_dd['CtoMyc'], places=9)  # mass balance


class TestWholeModelStability(unittest.TestCase):
    """Default parameter files, 2-year run on the climate fixture. Bacteria and mycorrhiza must stay within an order of
    magnitude of their observed (initial) biomass; with first-order turnover they do not."""

    @classmethod
    def setUpClass(cls):
        res = run_model(NDAYS, density_dependent=True)
        cls.B0, cls.C = res.B0, res.pools
        cls.C1 = run_model(NDAYS, density_dependent=False).pools

    def test_parameter_file_is_read(self):
        cfg = make_config()
        self.assertTrue(cfg.density_dependent_turnover)
        self.assertEqual(len(cfg.beta), 9)
        np.testing.assert_allclose(cfg.bref_effective, cfg.B0[:9])   # Bref <= 0 in the file -> initial biomass
        off = make_config(density_dependent=False)
        self.assertIsNone(off.beta_effective); self.assertIsNone(off.bref_effective)

    def test_microbial_pools_stay_within_plausible_range(self):
        C, B0 = self.C, self.B0
        for i, n in [(0, 'bacteria'), (2, 'mycorrhiza')]:
            ratio = C[:, i]/B0[i]
            self.assertLess(ratio.max(), 5, "%s peaks at %.1f x initial" % (n, ratio.max()))
            self.assertGreater(ratio.min(), 0.05, "%s drops to %.3f x initial" % (n, ratio.min()))
        # microbial C (bacteria + fungi + mycorrhiza) never exceeds 5 % of SOM C
        self.assertLess((C[:, :3].sum(1)/C[:, 10]).max(), 0.05)
        # no day-to-day explosion: biomass changes by less than 20 % per day
        for i in (0, 2):
            f = C[1:, i]/C[:-1, i]
            self.assertLess(f.max(), 1.2); self.assertGreater(f.min(), 0.8)

    def test_SOM_and_litter_remain_in_range(self):
        C = self.C
        self.assertGreater(C[:, 10].min(), 0.7*self.B0[10]); self.assertLess(C[:, 10].max(), 1.3*self.B0[10])
        self.assertGreater(C[:, 9].min(), 0.3*self.B0[9]);  self.assertLess(C[:, 9].max(), 1.5*self.B0[9])

    def test_first_order_turnover_is_far_less_stable(self):
        """Documents why the add-on is on by default: with beta = 1 the same run lets populations      
        run away or collapse over many orders of magnitude, instead of cycling around their
        observed biomass."""
        C, C1, B0 = self.C, self.C1, self.B0

        def orders(D, i):
            return np.log10(D[:, i].max()/max(D[:, i].min(), 1e-300))
            
        for i, n in [(0, 'bacteria'), (5, 'detritivores')]:
            self.assertGreater(orders(C1, i), orders(C, i) + 3,
                               '%s: %.1f orders with beta=1 vs %.1f with density dependence'
                               % (n, orders(C1, i), orders(C, i)))
            self.assertLess(orders(C, i), 2)
        # and some group leaves any plausible range entirely with first-order turnover
        self.assertGreater((C1[:, :9].max(0)/B0[:9]).max(), 10*(C[:, :9].max(0)/B0[:9]).max()/3)
        self.assertGreater((C1[:, :9].min(0)/B0[:9]).min()*1e6, 0)   # ... while others collapse 
        self.assertLess((C1[:, :9].min(0)/B0[:9]).min(), 1e-6)


if __name__ == '__main__':
    unittest.main()
