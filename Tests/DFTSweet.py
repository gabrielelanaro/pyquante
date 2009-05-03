#!/usr/bin/env python
"""\
UnitSweet.py - Unit testing for Python.

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
        self.li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)
        self.li_p = Molecule('Li+',atomlist = [(3,(0,0,0))],charge=1)
        self.li_m = Molecule('Li-',atomlist = [(3,(0,0,0))],charge=-1)
        self.h2o = Molecule('h2o',[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))],
                            units="Angstrom")
        self.oh = Molecule('oh',[(8,(0,0,0)),(1,(1.,0,0))],
                            units="Angstrom")
        self.lih = Molecule('LiH',[(1,(0,0,1.5)),(3,(0,0,-1.5))],units='Bohr')

    def testH2BLYP(self):
        h2_blyp = SCF(self.h2,method="DFT",functional='BLYP')
        h2_blyp.iterate()
        self.assertAlmostEqual(h2_blyp.energy,-1.166286,4)

    def testH2LDA(self):
        h2_lda = SCF(self.h2,method='DFT',functional="SVWN")
        h2_lda.iterate()
        self.assertAlmostEqual(h2_lda.energy,-1.132799,4)

    def testLiLDA(self):
        li_lda = SCF(self.li,method='DFT',functional="SVWN")
        li_lda.iterate()
        self.assertAlmostEqual(li_lda.energy,-7.332050,4)

    def testLiUHF(self):
        li_uhf = SCF(self.li,method='UHF')
        li_uhf.iterate()
        self.assertAlmostEqual(li_uhf.energy,-7.431364,4)

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
    
