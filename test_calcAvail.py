"""
Tests for the soil water modifier on substrate availability (keylink_functions.calcAvail).

Reference: Flores et al. (2021) PeerJ 9:e10707, Eqs. 14-16:
    mH2Otot = 4 (V/Pvol)(1 - V/Pvol)          for V/Pvol < 0.5, else 1        (Eq. 14)
    mH2O    = PvolA / (PvolA + PvolW) mH2Otot  for the partially filled class  (Eq. 15)
    mH2O    = PvolW / (PvolA + PvolW) mH2Otot  for the pore class one above    (Eq. 16)
with PvolW the water filled and PvolA the aerated pore volume of the class.

Run with:  python -m unittest test_calcAvail   (or pytest)
"""
import unittest
import numpy as np

import keylink_functions as mf

# pore volumes (l/m3) of the five classes: inaccessible, bacterial, micro, meso, macro
PV = np.array([45., 37., 37., 200., 6.])
SOMINI = 11470.


def expected_mwater(PV, PW):
    """Straightforward transcription of Eqs. 14-16 (for classes 0-3 and, via
    Eq. 16, the macropore class when mesopores are partially filled)."""
    V, Pvol = PW.sum(), PV.sum()
    mtot = 4*(V/Pvol)*(1-V/Pvol) if V/Pvol < 0.5 else 1.0
    mw = np.zeros(5)
    for i in range(4):
        PvolW = PW[i]
        PvolA = PV[i]-PW[i]
        if 0 < PvolW < PV[i]:
            mw[i] = PvolA/(PvolA+PvolW)*mtot
            if i == 3:
                mw[i+1] = PvolW/(PvolA+PvolW)*mtot
    return mw


def availSOMbact_from_mwater(PV, mw):
    """Bacteria can access all classes except the inaccessible one (volume weighted)."""
    return np.sum(PV[1:]*mw[1:])/PV[1:].sum()


class TestWaterModifier(unittest.TestCase):

    def test_partially_filled_class_eq15(self):
        """Eq. 15: with micropores filled to a fraction f the modifier of that
        class must be (1-f)*mH2Otot, i.e. PvolA/(PvolA+PvolW)*mH2Otot."""
        for f in (0.25, 0.5, 0.75, 0.9):
            PW = np.array([45., 37., f*37., 0., 0.])
            V_over_P = PW.sum()/PV.sum()
            mtot = 4*V_over_P*(1-V_over_P) if V_over_P < 0.5 else 1.0
            # only the micropore class contributes to availSOMbact here
            expected = PV[2]*(1-f)*mtot/PV[1:].sum()
            got = mf.calcAvail(PV, PW.copy(), SOMINI)[0]
            self.assertAlmostEqual(got, expected, places=12,
                                   msg="filling fraction %.2f" % f)

    def test_half_filled_gives_half_mtot(self):
        """The clearest check of the volume used: at 50 %% filling Eq. 15 gives
        exactly 0.5*mH2Otot (the old code gave 2/3*mH2Otot)."""
        PW = np.array([45., 37., 18.5, 0., 0.])
        V_over_P = PW.sum()/PV.sum()
        mtot = 4*V_over_P*(1-V_over_P)
        availSOMbact = mf.calcAvail(PV, PW.copy(), SOMINI)[0]
        mwater_micro = availSOMbact*PV[1:].sum()/PV[2]
        self.assertAlmostEqual(mwater_micro, 0.5*mtot, places=12)
        self.assertNotAlmostEqual(mwater_micro, PV[2]/(PV[2]+PW[2])*mtot, places=3)

    def test_modifier_decreases_monotonically_to_zero(self):
        """As the class fills up, the modifier must go continuously to 0 (the
        old formula was bounded below by 0.5*mH2Otot and then jumped to 0)."""
        prev = np.inf
        for f in np.linspace(0.05, 0.999, 40):
            PW = np.array([45., 37., f*37., 0., 0.])
            a = mf.calcAvail(PV, PW.copy(), SOMINI)[0]
            self.assertLess(a, prev)
            prev = a
        self.assertLess(prev, 1e-3)  # nearly full -> nearly zero
        PW = np.array([45., 37., 37., 0., 0.])  # exactly full -> exactly zero
        self.assertEqual(mf.calcAvail(PV, PW.copy(), SOMINI)[0], 0.0)

    def test_class_above_eq16(self):
        """Eq. 16: when mesopores are partially filled, macropores get
        PvolW/(PvolA+PvolW)*mH2Otot; the two modifiers sum to mH2Otot."""
        PW = np.array([45., 37., 37., 150., 0.])   # wet: mH2Otot = 1
        avail = mf.calcAvail(PV, PW.copy(), SOMINI)
        mw = expected_mwater(PV, PW)
        self.assertAlmostEqual(mw[3], 0.25)
        self.assertAlmostEqual(mw[4], 0.75)
        self.assertAlmostEqual(mw[3]+mw[4], 1.0)
        # availSOMbact, availSOMfungi, availSOMeng, availSOMsap all see meso+macro
        for idx in range(4):
            self.assertAlmostEqual(avail[idx], availSOMbact_from_mwater(PV, mw) if idx in (0, 2)
                                   else np.sum(PV[1:]*mw[1:]*([0, 1, 1, 1] if idx == 1 else [0, 0, 1, 1]))/PV[1:].sum(),
                                   places=12)

    def test_full_availability_vector_against_reference(self):
        """All availability entries must equal a direct transcription of
        Eqs. 14-16 for several water states (dry, medium, wet)."""
        masks = {  # pore classes each consumer can use (indices into PV[1:] or PV[2:] / PV[3:])
            0: ([1, 1, 1, 1], 1), 1: ([0, 1, 1, 1], 1), 2: ([1, 1, 1, 1], 1), 3: ([0, 0, 1, 1], 1),
            4: ([0, 1, 1, 1], 1), 5: ([1, 1, 1], 2), 6: ([0, 1, 1], 2), 7: ([0, 1, 1], 2),
            8: ([0, 1, 1], 2), 9: ([1, 1], 3)}
        for PW in (np.array([45., 10., 0., 0., 0.]),
                   np.array([45., 37., 30., 0., 0.]),
                   np.array([45., 37., 37., 120., 0.])):
            got = mf.calcAvail(PV, PW.copy(), SOMINI)
            mw = expected_mwater(PV, PW)
            for k, (mask, start) in masks.items():
                exp = np.sum(np.array(mask)*PV[start:]*mw[start:])/PV[start:].sum()
                self.assertAlmostEqual(got[k], exp, places=12, msg="entry %d, PW=%s" % (k, PW))
            self.assertEqual(got[10], 1)                          # engineers always available
            self.assertAlmostEqual(got[11], PV[0]/PV.sum()*SOMINI)  # SOM in inaccessible pores


if __name__ == '__main__':
    unittest.main()
