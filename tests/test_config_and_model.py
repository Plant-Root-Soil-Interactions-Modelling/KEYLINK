"""
Tests of the model API: KeylinkConfig (explicit, validated, immutable configuration),
KeylinkModel (no shared state between runs) and KeylinkResults, plus the matrix infiltration
capacity used by calcPW.

Run with:  python run_tests.py   (all tests) from the model directory.
"""
import unittest
import numpy as np

from tests.helpers import mf, core, make_config, FIXTURE_CLIMATE, FIXTURE_DAYS


class TestKeylinkConfig(unittest.TestCase):

    def test_from_files_reads_all_parameter_files(self):
        cfg = make_config()
        P = np.loadtxt('KL_FaunalParams.txt', encoding='latin1')
        np.testing.assert_array_equal(cfg.gmax, P[0]); np.testing.assert_array_equal(cfg.q10, P[12])
        soil = np.loadtxt('KL_initSoil.txt', encoding='latin1')
        self.assertEqual(cfg.m, soil[4]); self.assertEqual(cfg.pH, soil[6])
        np.testing.assert_array_equal(cfg.PVstruct, soil[10:15])
        np.testing.assert_array_equal(cfg.B0, np.loadtxt('KL_initC_Pools.txt', encoding='latin1')[:13])
        self.assertEqual(cfg.tStop, FIXTURE_DAYS); self.assertEqual(cfg.climatefile, FIXTURE_CLIMATE)
        self.assertEqual(cfg.pfaec[5], np.loadtxt('KL_engineerParams.txt', encoding='latin1')[5])

    def test_default_run_length_comes_from_runparams(self):
        cfg = core.KeylinkConfig.from_files('.')
        self.assertEqual(cfg.tStop, int(np.loadtxt('KL_runparams.txt', encoding='latin1')[0]))
        self.assertEqual(cfg.climatefile, 'PrecipKMIBrass.txt')

    def test_config_is_immutable_and_replace_gives_validated_copy(self):
        cfg = make_config()
        with self.assertRaises(Exception):
            cfg.pH = 5.0
        cfg2 = cfg.replace(pH=5.9, tStop=10)
        self.assertEqual(cfg2.pH, 5.9); self.assertEqual(cfg2.tStop, 10)
        self.assertEqual(cfg.pH, 3.9); self.assertEqual(cfg.tStop, FIXTURE_DAYS)
        with self.assertRaises(ValueError):
            cfg.replace(gmax=np.ones(8))
        with self.assertRaises(ValueError):
            cfg.replace(tStop=0)
        with self.assertRaises(ValueError):
            cfg.replace(beta=np.full(9, 0.5))

    def test_density_dependence_switch(self):
        on = make_config(); off = make_config(density_dependent=False)
        self.assertIsNotNone(on.beta_effective); self.assertIsNone(off.beta_effective)
        np.testing.assert_allclose(on.bref_effective, on.B0[:9])
        explicit = on.replace(bref=np.full(9, 2.0))
        np.testing.assert_allclose(explicit.bref_effective, 2.0)


class TestKeylinkModel(unittest.TestCase):

    def test_import_has_no_side_effects(self):
        """Importing keylink_core must not read files or create module-level model state."""
        for name in ('B', 'Bini', 'PW', 'PWini', 'psoln', 'tStop', 'Ntot', 'Nmin', 'GMAX', 'CN', 'fEight'):
            self.assertFalse(hasattr(core, name), 'module-level state %s still exists' % name)

    def test_repeated_runs_are_identical_and_independent(self):
        cfg = make_config(60)
        model = core.KeylinkModel(cfg)
        r1 = model.run(); r2 = model.run()
        np.testing.assert_array_equal(r1.pools, r2.pools)
        np.testing.assert_array_equal(r1.soil_water, r2.soil_water)
        np.testing.assert_array_equal(r1.Ntot, r2.Ntot)
        self.assertIsNot(r1.pools, r2.pools)
        r1.pools[:] = -1                                   # mutating one result does not touch the other
        self.assertFalse(np.any(r2.pools < 0))
        np.testing.assert_array_equal(cfg.B0, np.loadtxt('KL_initC_Pools.txt', encoding='latin1')[:13])  # config untouched

    def test_different_gmax_give_different_results_and_are_recorded(self):
        model = core.KeylinkModel(make_config(30))
        g = model.cfg.gmax.copy(); g2 = g.copy(); g2[0] *= 1.5
        r1 = model.run(g); r2 = model.run(g2)
        np.testing.assert_array_equal(r2.gmax, g2)
        self.assertGreater(np.abs(r1.pools[:, 0]-r2.pools[:, 0]).max(), 0)
        with self.assertRaises(ValueError):
            model.run(np.ones(5))

    def test_results_shapes_and_initial_day(self):
        cfg = make_config(20)
        r = core.KeylinkModel(cfg).run()
        self.assertEqual(r.ndays, 20)
        self.assertEqual(r.pools.shape, (20, 22)); self.assertEqual(r.pores.shape, (20, 5))
        for name in ('soil_water', 'porosity', 'Nmin', 'Ntot', 'SOMCN', 'litterCN', 'mineralisation'):
            self.assertEqual(getattr(r, name).shape, (20,))
        np.testing.assert_array_equal(r.B0, cfg.B0)
        self.assertAlmostEqual(r.Ntot0, (cfg.B0[:9]/cfg.cn).sum() + cfg.Nmin0 + cfg.B0[10]/cfg.SOMCN0 + cfg.B0[9]/cfg.litterCN0)
        self.assertEqual(r.errors, [])
        self.assertEqual(r.pools[0, 21], 0.0)             # spare entry of the state vector

    def test_climate_file_too_short_is_an_error(self):
        cfg = make_config().replace(tStop=FIXTURE_DAYS+1)
        with self.assertRaises(ValueError):
            core.KeylinkModel(cfg).run()

    def test_van_genuchten_m_is_the_parameter_file_value(self):
        """Regression: the old module overwrote the van Genuchten m with a loop variable (m = 11)."""
        cfg = make_config()
        self.assertEqual(cfg.m, 0.3)


class TestMatrixInfiltration(unittest.TestCase):
    """calcPW uses a constant maximal matrix infiltration rate (Ksat), as prescribed by Flores
    et al. (2021). With the former unsaturated conductivity K(theta) a dry soil could never rewet."""

    def test_dry_matrix_takes_up_rain(self):
        PV = np.array([45., 37., 37., 200., 6.])
        PW = np.array([45., 0., 0., 0., 0.])            # matrix completely dry
        PWn, drain, runoff = mf.calcPW(PV, 8.0, PW.copy(), 7., 1., 45., 319., 1.2, 1.7, 0.3, 10.)
        self.assertGreater(PWn[1]+PWn[2]+PWn[3], 0.1)    # some of the 8 mm entered the matrix
        self.assertEqual(runoff, 0.0)


if __name__ == '__main__':
    unittest.main()
