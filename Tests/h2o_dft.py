#!/usr/bin/env python
"H2O using Gaussians to test the DFT module"

import unittest, sciunittest

from PyQuante.dft import *
from PyQuante.Molecule import Molecule

# Getting -75.857635 with the old style grid patching
# Getting -75.905545 with the Becke projection operator.
# Getting -75.896366 with Becke projection and heteroatom corrections
# What does Jaguar give?

# Changed 2005-11-04 to reflect Becke projection + heteroatom corrections
# Changed 2005-11-08 to remove the spingrid option
energy = -75.896499
name = "H2O/DFT"

def main():
    r = 1./0.5291772                  # conversion factor modified
    h2o=Molecule('h2o',atomlist = [(8,(0,0,0)),(1,(r,0,0)),(1,(0,r,0))])
    en,orbe,orbs = dft(h2o)
    return en

class WaterDFTTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of H2O (DFT) close to -75.896499?"""
        E = main()
        self.assertInside(E, energy, 1e-4)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(WaterDFTTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())
