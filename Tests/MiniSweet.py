#!/usr/bin/env python
"""\
MiniSweet.py - Smaller version of UnitSweet that focuses on the
               routines I think are rate-determining.

<beavis>heh, heh, he said *unit*</beavis>
"""

import unittest,logging
from PyQuante.CI import CIS
from PyQuante.Molecule import Molecule
from PyQuante.MP import MP2
from PyQuante.OEP import oep_hf,oep_hf_an
from PyQuante.PyQuante2 import SCF,SubspaceSolver,DmatSolver

class UnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],
                           units='Angs')
        self.he = Molecule('He',atomlist = [(2,(0,0,0))])

    def testH2HF(self):
        h2_hf = SCF(self.h2,method='HF')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testH2BLYP(self):
        h2_blyp = SCF(self.h2,method="DFT",functional='BLYP')
        h2_blyp.iterate()
        self.assertAlmostEqual(h2_blyp.energy,-1.166286,4)

    def testHeHF(self):
        he_hf = SCF(self.he,method='HF')
        he_hf.iterate()
        self.assertAlmostEqual(he_hf.energy,-2.855260,3)

    def testH2LDA(self):
        h2_lda = SCF(self.h2,method='DFT',functional="SVWN")
        h2_lda.iterate()
        self.assertAlmostEqual(h2_lda.energy,-1.132799,4)


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
    import profile,pstats
    profile.run('runsuite()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__':
    #runsuite()
    debugsuite()
