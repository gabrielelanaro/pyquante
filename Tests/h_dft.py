#!/usr/bin/env python
import unittest, sciunittest

from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

name = "H atom/DFT"
energy = -0.4415033 # open shell average occupation

def main():
    atomlist = Molecule('H',atomlist = [(1,(0,0,0))],multiplicity=2)
    en,orbe,orbs = dft(atomlist)
    return en

class HDFTTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of H atom (DFT) close to -0.4415033?"""
        E = main()
        self.assertInside(E, energy, 1e-5)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(HDFTTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())
