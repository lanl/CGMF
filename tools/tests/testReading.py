import unittest
import os
#import sys
#import site
#sys.path.append('/User/lovell/Library/Python/3.7/lib/python/site-packages/CGMFtk-1.0-py3.7.egg/CGMFtk/')
from CGMFtk import histories as fh
from CGMFtk import yields as yld

# will need to put a file to test in this folder as well #

class TestHistoryFileRead(unittest.TestCase):
    
    def test_nevents_read(self):
        """
        test that the correct number of events has been read
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        self.assertEqual(h.numberEvents,nread)

    def test_number_of_fragments(self):
        """
        test that the number of fragments is twice the number of events
        """
        nread = 10
        #h = fh.Histories('92235_18MeV.cgmf',nevents=nread)
        h = fh.Histories('../../utils/cgmf/tests/u235nf-18MeV-events/histories.cgmf.parallel.0.reference',nevents=nread)
        self.assertEqual(h.numberFragments,2*nread)

    


class TestYieldFileRead(unittest.TestCase):

    def test_nevents_read(self):
        """
        test that the correct number of events has been read
        """
        nread = 10
        #h = yld.Yields('92235_18MeV_yields.cgmf',nevents=nread)
        h = yld.Yields('../../utils/cgmf/tests/u235nf-10MeV-yields/yields.cgmf.serial.0.reference',nevents=nread)
        self.assertEqual(h.numberEvents,nread)

    def test_number_of_fragments(self):
        """
        test that the number of fragments is twice the number of events
        """
        nread = 10
        #h = yld.Yields('92235_18MeV_yields.cgmf',nevents=nread)
        h = yld.Yields('../../utils/cgmf/tests/u235nf-10MeV-yields/yields.cgmf.serial.0.reference',nevents=nread)
        self.assertEqual(h.numberFragments,2*nread)

        

if __name__ == '__main__':
    unittest.main()
