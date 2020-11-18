import unittest
import os
import numpy as np
from CGMFtk import histories as fh

class TestGammaTimings(unittest.TestCase):

    # test the time coincidence window

    def test_time_coincidence_window(self):
        """
        test whether the timing window is correct (-1 when all timing is included)
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        #ht = fh.Histories('92235_18MeV_timing.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-th-events/histories.cgmf.parallel.0.reference',nevents=nread)
        ht = fh.Histories('../../utils/cgmf/tests/u235nf-time-events/histories.cgmf.parallel.0.reference',nevents=nread)
        self.assertTrue(h.getTimeCoincidenceWindow()>0)
        self.assertEqual(ht.getTimeCoincidenceWindow(),1.5e-7)

    def test_gamma_ages(self):
        """
        ensure the gamma times exist when they should be calculated
        """
        nread = 10
        #h = fh.Histories('92235_18MeV_timing.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-time-events/histories.cgmf.parallel.0.reference',nevents=nread)
        self.assertTrue(len(h.getGammaAges())>0)


if __name__ == '__main__':
    unittest.main()
