"""
Tests for the soil water modifier on substrate availability (keylink_functions.calcAvail).

Reference: Flores et al. (2021) PeerJ 9:e10707, Eqs. 14-16:
    mH2Otot = 4 (V/Pvol)(1 - V/Pvol)          for V/Pvol < 0.5, else 1        (Eq. 14)
    mH2O    = PvolA / (PvolA + PvolW) mH2Otot  for the partially filled class  (Eq. 15)
    mH2O    = PvolW / (PvolA + PvolW) mH2Otot  for the pore class one above    (Eq. 16)
with PvolW the water filled and PvolA the aerated pore volume of the class. All other
classes (completely dry or completely flooded) have mH2O = 0.

Run with:  python run_tests.py   (all tests) from the model directory.
"""
import unittest
import numpy as np

from tests.helpers import mf

# pore volumes (l/m3) of the five classes: inaccessible, bacterial, micro, meso, macro
PV = np.array([45., 37., 37., 200., 6.])
SOMINI = 11470.


def expected_mwater(PV, PW):
    """Straightforward transcription of Eqs. 14-16. The water level is the largest matrix
    class holding water (macropores only count when the mesopores are full); that class
    gets Eq. 15, the class above it Eq. 16. A class that is exactly full therefore passes
    all activity to the class above."""
    V, Pvol = PW.sum(), PV.sum()
    mtot = 4*(V/Pvol)*(1-V/Pvol) if V/Pvol < 0.5 else 1.0
    mw = np.zeros(5)
    filled = [i for i in range(5) if PW[i] > 0 and (i < 4 or PW[3] >= PV[3])]
    if filled:
        k = max(filled)
        PvolW, PvolA = PW[k], PV[k]-PW[k]
        mw[k] = PvolA/(PvolA+PvolW)*mtot
        if k < 4:
            mw[k+1] = PvolW/(PvolA+PvolW)*mtot
    return mw


def availSOMbact_from_mwater(PV, mw):
    """Bacteria can access all classes except the inaccessible one (volume weighted)."""
    return np.sum(PV[1:]*mw[1:])/PV[1:].sum()


class TestWaterModifier(unittest.TestCase):

    def test_partially_filled_class_eq15_and_16(self):
        """Micropores filled to a fraction f: the micropore modifier is (1-f)*mH2Otot (Eq. 15)
        and the mesopores, one class above, get f*mH2Otot (Eq. 16)."""
        for f in (0.25, 0.5, 0.75, 0.9):
            PW = np.array([45., 37., f*37., 0., 0.])
            V_over_P = PW.sum()/PV.sum()
            mtot = 4*V_over_P*(1-V_over_P) if V_over_P < 0.5 else 1.0
            expected = (PV[2]*(1-f)*mtot + PV[3]*f*mtot)/PV[1:].sum()
            got = mf.calcAvail(PV, PW.copy(), SOMINI)[0]
            self.assertAlmostEqual(got, expected, places=12, msg="filling fraction %.2f" % f)
            # detritivores (meso + macro only) see exactly the Eq. 16 share
            self.assertAlmostEqual(mf.calcAvail(PV, PW.copy(), SOMINI)[3], PV[3]*f*mtot/PV[1:].sum(), places=12)

    def test_half_filled_gives_half_mtot(self):
        """The clearest check of the volume used: at 50 % filling Eq. 15 gives exactly
        0.5*mH2Otot for the filling class (the old code gave 2/3*mH2Otot)."""
        PW = np.array([45., 37., 18.5, 0., 0.])
        V_over_P = PW.sum()/PV.sum()
        mtot = 4*V_over_P*(1-V_over_P)
        mw = expected_mwater(PV, PW)
        self.assertAlmostEqual(mw[2], 0.5*mtot, places=12)
        self.assertAlmostEqual(mw[3], 0.5*mtot, places=12)
        self.assertAlmostEqual(mf.calcAvail(PV, PW.copy(), SOMINI)[0], availSOMbact_from_mwater(PV, mw), places=12)
        self.assertNotAlmostEqual(mw[2], PV[2]/(PV[2]+PW[2])*mtot, places=3)

    def test_modifier_is_continuous_across_class_boundaries(self):
        """As the water content rises through a class boundary the availability must not jump:
        a full micropore class (f=1) hands all activity to the mesopores, which is the same
        state as mesopores 'just starting to fill'."""
        prev = None
        for V in np.linspace(83., 250., 300):   # from 'micropores just filling' to 'mesopores 80 % full'
            PW = np.array([45., 37., min(37., V-82.), max(0., V-119.), 0.])
            a = mf.calcAvail(PV, PW.copy(), SOMINI)[0]
            if prev is not None:
                self.assertLess(abs(a-prev), 0.02, "jump at V=%.1f" % V)
            prev = a
        # full micropores, empty mesopores (field capacity): activity is in the mesopores, not zero
        PW = np.array([45., 37., 37., 0., 0.])
        a = mf.calcAvail(PV, PW.copy(), SOMINI)
        r = PW.sum()/PV.sum(); mtot = 4*r*(1-r) if r < 0.5 else 1.0
        self.assertAlmostEqual(a[0], PV[3]/PV[1:].sum()*mtot, places=12)
        self.assertGreater(a[3], 0)                                        # detritivores active too

    def test_macropore_water_only_counts_when_mesopores_full(self):
        PW = np.array([45., 37., 20., 0., 3.])        # transient macropore water, matrix level in micropores
        mw = expected_mwater(PV, PW)
        self.assertGreater(mw[2], 0); self.assertGreater(mw[3], 0); self.assertEqual(mw[4], 0)
        np.testing.assert_allclose(mf.calcAvail(PV, PW.copy(), SOMINI)[0], availSOMbact_from_mwater(PV, mw))
        PW = np.array([45., 37., 37., 200., 3.])      # saturated matrix: macropores are the filling class
        mw = expected_mwater(PV, PW)
        self.assertEqual(mw[3], 0); self.assertAlmostEqual(mw[4], (6.-3.)/6.)
        np.testing.assert_allclose(mf.calcAvail(PV, PW.copy(), SOMINI)[0], availSOMbact_from_mwater(PV, mw))

    def test_class_above_eq16(self):
        """Eq. 16: when mesopores are partially filled, macropores get
        PvolW/(PvolA+PvolW)*mH2Otot; the two modifiers sum to mH2Otot."""
        PW = np.array([45., 37., 37., 150., 0.])   # wet: mH2Otot = 1
        mw = expected_mwater(PV, PW)
        self.assertAlmostEqual(mw[3], 0.25)
        self.assertAlmostEqual(mw[4], 0.75)
        self.assertAlmostEqual(mw[3]+mw[4], 1.0)
        got = mf.calcAvail(PV, PW.copy(), SOMINI)
        self.assertAlmostEqual(got[0], availSOMbact_from_mwater(PV, mw), places=12)

    def test_full_availability_vector_against_reference(self):
        """All availability entries must equal a direct transcription of
        Eqs. 14-16 for several water states (dry, medium, wet)."""
        masks = {  # pore classes each consumer can use (indices into PV[1:] or PV[2:] / PV[3:])
            0: ([1, 1, 1, 1], 1), 1: ([0, 1, 1, 1], 1), 2: ([1, 1, 1, 1], 1), 3: ([0, 0, 1, 1], 1),
            4: ([0, 1, 1, 1], 1), 5: ([1, 1, 1], 2), 6: ([0, 1, 1], 2), 7: ([0, 1, 1], 2),
            8: ([0, 1, 1], 2), 9: ([1, 1], 3)}
        for PW in (np.array([45., 10., 0., 0., 0.]),
                   np.array([45., 37., 30., 0., 0.]),
                   np.array([45., 37., 37., 0., 0.]),
                   np.array([45., 37., 37., 120., 0.]),
                   np.array([45., 37., 37., 200., 2.])):
            got = mf.calcAvail(PV, PW.copy(), SOMINI)
            mw = expected_mwater(PV, PW)
            for k, (mask, start) in masks.items():
                exp = np.sum(np.array(mask)*PV[start:]*mw[start:])/PV[start:].sum()
                self.assertAlmostEqual(got[k], exp, places=12, msg="entry %d, PW=%s" % (k, PW))
            self.assertEqual(got[10], 1)                          # engineers always available
            self.assertAlmostEqual(got[11], PV[0]/PV.sum()*SOMINI)  # SOM in inaccessible pores


if __name__ == '__main__':
    unittest.main()
