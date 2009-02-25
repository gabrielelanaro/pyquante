#!/usr/bin/env python
import unittest, sciunittest

from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

name = "Li atom/DFT"
energy = -7.3321  # open shell average occupation

def main():
    atomlist = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)
    en,orbe,orbs = dft(atomlist,verbose=True)
    return en

class LiDFTTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of Li atom (DFT) close to -7.3321?"""
        E = main()
        self.assertInside(E, energy, 1e-4)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(LiDFTTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())                                   

