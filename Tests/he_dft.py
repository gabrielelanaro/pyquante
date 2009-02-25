#!/usr/bin/env python
"He using Gaussians to test the DFT module"

import unittest, sciunittest

from PyQuante.dft import *
from PyQuante.Molecule import Molecule

energy = -2.8266976
name = "He"

def main():
    en,orbe,orbs = dft(Molecule('He',atomlist=[(2,(0,0,0))]))
    return en

class HeDFTTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of He/DFT close to -2.82?"""
        E = main()
        self.assertInside(E, energy, 1e-6)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(HeDFTTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())
