import unittest
import os
import numpy as np
from CGMFtk import histories as fh

class TestNeutronProperties(unittest.TestCase):

    # nubar calculations - for the event and for the fragments (doesn't include pre-fission neutrons)

    def test_nubar_for_events(self):
        """
        test the nubartot function calculation (average neutrons per event)
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        nuLF = h.getNuLF()
        nuHF = h.getNuHF()
        nuPF = h.getPreFissionNu()
        nutot = nuLF + nuHF + nuPF
        nutot = np.mean(nutot)
        self.assertEqual(h.nubartot(),nutot)

    def test_nubar_for_fragments(self):
        """
        test the nubar function calculation (average neutrons per fragment)
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        nuLF = h.getNuLF()
        nuHF = h.getNuHF()
        nutot = np.sum(nuLF)+np.sum(nuHF)
        nutot = nutot/(2.*len(nuLF))
        self.assertEqual(h.nubar(),nutot)

    # mean neutron energy calculations

    def test_mean_neutron_energy_in_lab_frame(self):
        """
        test the calculation of the average neutron energy in the lab frame
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        nE = h.getNeutronElab()
        nPFE = h.getPreFissionNeutronElab()
        nuPF = h.getPreFissionNu()
        l = []
        for x in nE:
            l+=x
        for i in range(len(nuPF)):
            if (nuPF[i]>0):
                l+=nPFE[i]
        self.assertEqual(np.mean(l),h.meanNeutronElab())

    # mean list test covers the rest of the neutron avearges -- they all call this function

    def test_mean_list_calculation(self):
        """
        test the calculation of the mean of a list
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        nE = h.getNeutronEcm()
        l = []
        for x in nE:
            l+=x
        self.assertEqual(np.mean(l),h.meanList(nE))

    # P(nu)

    def test_neutron_multiplicity_distribution(self):
        """
        test the P(nu) calculation
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        nu,pnu = h.Pnu()
        nuvals = h.getNutot()
        numax = np.max(nuvals)
        self.assertEqual(nu[-1],numax)
        for i in range(numax+1):
            prob = (len(nuvals[nuvals==i]))/float(nread)
            self.assertEqual(prob,pnu[i])

    # energy thresholds

    def test_neutron_spectrum_energy_threshold(self):
        """
        test the implementation of the neutron energy threshold
        specificially in the PFNS, but is constructed the same elsewhere
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        Eth = 0.5
        eng,pfns = h.pfns(Eth=Eth)
        eng0,pfns0 = h.pfns()
        emin = eng[pfns>0]
        emin0 = eng0[pfns0>0]
        if (Eth>0):
            self.assertTrue(emin0[0]<emin[0])
        else:
            self.assertEqual(emin0[0],emin[0])
        



if __name__ == '__main__':
    unittest.main()
