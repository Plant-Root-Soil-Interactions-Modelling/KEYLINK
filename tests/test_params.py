"""
Tests of the schema-driven parameter file reader (keylink_params): every line of a KL_*.txt file
is identified by the keyword in its comment and checked against its valid range, so that a
missing, reordered, mislabelled or out-of-range line is reported instead of silently shifting the
following values.

Run with:  python run_tests.py   (all tests) from the model directory.
"""
import os
import shutil
import tempfile
import unittest
import numpy as np

from tests.helpers import core, make_config
import keylink_params as kp
import keylink_names as kn


def copy_params(tmp):
    for f in os.listdir('.'):
        if f.startswith('KL_') and f.endswith('.txt'):
            shutil.copy(f, tmp)


def edit(tmp, name, fn):
    """Rewrite <tmp>/<name>.txt with fn(list of lines) -> list of lines."""
    path = os.path.join(tmp, name + '.txt')
    lines = open(path, encoding='latin1').read().splitlines()
    open(path, 'w', encoding='latin1').write('\n'.join(fn(lines)) + '\n')


class TestSchema(unittest.TestCase):

    def test_all_repository_files_load(self):
        for name in kp.SCHEMA:
            d = kp.load_parameter_file(name)
            self.assertEqual(list(d), [e.key for e in kp.SCHEMA[name]['entries']], name)

    def test_group_wise_files_return_one_array_per_parameter(self):
        d = kp.load_parameter_file('KL_FaunalParams')
        for k, v in d.items():
            self.assertEqual(v.shape, (kn.NGROUPS,), k)
        self.assertEqual(d['gmax'][kn.BACT], 1.24); self.assertEqual(d['ks'][kn.PRED], 2.0)
        self.assertEqual(d['q10'][kn.ENG], 2.5)

    def test_scalar_files_return_floats_and_names(self):
        soil = kp.load_parameter_file('KL_initSoil')
        self.assertEqual(soil['m'], 0.3); self.assertEqual(soil['Ksat'], 10.0); self.assertEqual(soil['PV_meso'], 200.0)
        run = kp.load_parameter_file('KL_runparams')
        self.assertEqual(run['tStop'], 3653.0); self.assertEqual(run['CtoMyc'], 0.54)
        pools = kp.load_parameter_file('KL_initC_Pools')
        self.assertEqual(pools['SOM'], 11470.0); self.assertEqual(pools['pred'], 0.4)
        clim = kp.load_parameter_file('KL_climateParams')
        self.assertEqual(clim['TEMP'].shape, (12,)); self.assertEqual(clim['SUNH'][6], 201.0)

    def test_load_all_and_describe(self):
        d = kp.load_all()
        self.assertEqual(set(d), {'FaunalParams', 'initSoil', 'engineerParams', 'runparams', 'initC_Pools',
                                  'climateParams', 'densityParams'})
        text = kp.describe('KL_initSoil')
        self.assertIn('Ksat', text); self.assertIn('[0, 1]', text)   # m has range [0, 1]


class TestErrorsAreReported(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        copy_params(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def assertRaisesMessage(self, fragment, fn):
        with self.assertRaises(kp.ParameterFileError) as ctx:
            fn()
        self.assertIn(fragment, str(ctx.exception))

    def test_missing_file(self):
        os.remove(os.path.join(self.tmp, 'KL_initSoil.txt'))
        self.assertRaisesMessage('not found', lambda: kp.load_parameter_file('KL_initSoil', self.tmp))

    def test_missing_line_is_detected(self):
        edit(self.tmp, 'KL_initSoil', lambda L: L[:5] + L[6:])      # drop Ksat
        self.assertRaisesMessage('expected 15 parameter lines', lambda: kp.load_parameter_file('KL_initSoil', self.tmp))

    def test_reordered_lines_are_detected(self):
        def swap(L):
            L = list(L); L[6], L[7] = L[7], L[6]; return L            # pH <-> litter CN
        edit(self.tmp, 'KL_initSoil', swap)
        self.assertRaisesMessage("should be 'pH'", lambda: kp.load_parameter_file('KL_initSoil', self.tmp))

    def test_mislabelled_line_is_detected(self):
        edit(self.tmp, 'KL_runparams', lambda L: [L[0]] + ['0.9\t#something else'] + L[2:])
        self.assertRaisesMessage("should be 'initWater'", lambda: kp.load_parameter_file('KL_runparams', self.tmp))

    def test_out_of_range_value_is_detected(self):
        edit(self.tmp, 'KL_initSoil', lambda L: L[:4] + ['11\t#m, van Genuchten'] + L[5:])   # the old m = 11 accident
        self.assertRaisesMessage('outside the valid range [0, 1]', lambda: kp.load_parameter_file('KL_initSoil', self.tmp))
        edit(self.tmp, 'KL_initSoil', lambda L: L[:4] + ['0.3\t#m, van Genuchten'] + L[5:])
        edit(self.tmp, 'KL_FaunalParams', lambda L: L[:5] + ['0 0 0 0 0 1.2 0.3 0.3 0.3\t# faec, fraction 0-1'] + L[6:])
        self.assertRaisesMessage('faec', lambda: kp.load_parameter_file('KL_FaunalParams', self.tmp))

    def test_wrong_group_order_in_header_is_detected(self):
        edit(self.tmp, 'KL_FaunalParams', lambda L: ['#fung\tbact\tEM\tbacVor\tfunVor\tdetrVor\tengin\therbVor\tpred'] + L[1:])
        self.assertRaisesMessage('group order', lambda: kp.load_parameter_file('KL_FaunalParams', self.tmp))

    def test_wrong_number_of_columns_is_detected(self):
        edit(self.tmp, 'KL_FaunalParams', lambda L: [L[0], '1.24 0.6 0.44 1.4 0.8 0.178 0.109 0.135\t# gmax'] + L[2:])
        self.assertRaisesMessage('expected 9 value(s)', lambda: kp.load_parameter_file('KL_FaunalParams', self.tmp))

    def test_non_numeric_token_is_detected(self):
        edit(self.tmp, 'KL_runparams', lambda L: ['abc\t#number of days of run'] + L[1:])
        self.assertRaisesMessage('cannot read number', lambda: kp.load_parameter_file('KL_runparams', self.tmp))

    def test_errors_propagate_through_config(self):
        edit(self.tmp, 'KL_engineerParams', lambda L: L[:-1])
        with self.assertRaises(kp.ParameterFileError):
            core.KeylinkConfig.from_files(self.tmp)


class TestConfigValidationUsesSchema(unittest.TestCase):

    def test_replace_checks_ranges(self):
        cfg = make_config()
        for field, value in (('m', 11.0), ('pH', 15.0), ('ee', 1.5), ('fPVB', -0.1), ('NmyctoPlant', 2.0)):
            with self.assertRaises(ValueError, msg=field):
                cfg.replace(**{field: value})
        with self.assertRaises(ValueError):
            cfg.replace(faec=np.full(9, 1.5))
        with self.assertRaises(ValueError):
            cfg.replace(topt=cfg.tmax + 1)       # tmin < topt < tmax
        cfg.replace(m=0.41, pH=5.9)            # valid values pass

    def test_named_pools_and_fluxes_in_results(self):
        r = core.KeylinkModel(make_config(10)).run()
        np.testing.assert_array_equal(r.pool('SOM'), r.pools[:, kn.SOM])
        np.testing.assert_array_equal(r.pool(kn.Pool.BACT), r.pools[:, 0])
        np.testing.assert_array_equal(r.flux('SOMeaten'), r.pools[:, kn.SOM_EATEN])
        np.testing.assert_array_equal(r.flux(kn.Flux.BACT_RESP), r.pools[:, 13])
        d = r.as_dict()
        self.assertEqual(len(d), 13 + 8)
        np.testing.assert_array_equal(d['LITeatenEng'], r.pools[:, 20])
        self.assertEqual([int(p) for p in kn.Pool], list(range(13)))
        self.assertEqual(kn.POOL_KEYS[kn.ENG], 'eng')


if __name__ == '__main__':
    unittest.main()
