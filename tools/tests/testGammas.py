import unittest
import os
import numpy as np
from CGMFtk import histories as fh

class TestGammaProperties(unittest.TestCase):

    # nubarg calculations

    def test_nubarg_for_events(self):
        """
        test the nubargtot function calculation (average gammas per event)
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        nug = h.getNugtot()
        self.assertEqual(np.mean(nug),h.nubargtot())

    def test_nubarg_for_fragments(self):
        """
        test the nubarg function calculation (average gammas per fragment)
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        nug = h.getNug()
        self.assertEqual(np.mean(nug),h.nubarg())

    # we tested the mean list calculation in the neutron tests

    # P(nug) tests

    def test_gamma_multiplicity_distribution(self):
        """
        test the gamma multiplicity distribution, P(N_gamma)
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        nug,pnug = h.Pnug()
        nugvals = h.getNugtot()
        nugmax = np.max(nugvals)
        self.assertEqual(nug[-1],nugmax)
        for i in range(nugmax+1):
            prob = (len(nugvals[nugvals==i]))/float(nread)
            self.assertEqual(prob,pnug[i])

    # energy thresholds

    def test_gamma_spectrum_energy_threshold(self):
        """
        test the implementation of the gamma-ray energy thresholds
        specifically in the PFGS, but is constructed the same elsewhere
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        Eth = 0.5
        eng,pfgs = h.pfgs(Eth=Eth)
        eng0,pfgs0 = h.pfgs()
        emin = eng[pfgs>0]
        emin0 = eng0[pfgs0>0]
        if (Eth>0):
            self.assertTrue(emin0[0]<emin[0])
        else:
            self.assertEqual(emin0[0],emin[0])


if __name__ == '__main__':
    unittest.main()
