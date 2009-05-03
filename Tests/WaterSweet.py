#!/usr/bin/env python
import unittest,logging

from PyQuante.Molecule import Molecule
from PyQuante.PyQuante2 import SCF

class UnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.h2o = Molecule('h2o',[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))],
                            units="Angstrom")

    def testH2OHF(self):
        h2o_hf = SCF(self.h2o,method='HF')
        h2o_hf.iterate()
        self.assertAlmostEqual(h2o_hf.energy,-76.011755864850628,4)

def runsuite(verbose=True):
    # To use psyco, uncomment this line:
    #import psyco; psyco.full()
    if verbose: verbosity=2
    else: verbosity=1
    # If you want more output, uncomment this line:
    #logging.basicConfig(format="%(message)s",level=logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(UnitTests)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)
    # Running without verbosity is equivalent to replacing the above
    # two lines with the following:
    #unittest.main()
    return

def debugsuite():
    import cProfile,pstats
    cProfile.run('runsuite()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__':
    import sys
    if "-d" in sys.argv:
        debugsuite()
    else:
        runsuite()
    
