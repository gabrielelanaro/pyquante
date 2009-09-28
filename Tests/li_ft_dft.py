#!/usr/bin/env python
import unittest, sciunittest

from PyQuante.NumWrap import arange
from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

energy = -7.349422
name = "Li/FT/DFT"

def main():
    li = Molecule('Li',
                  atomlist = [(3,(0,0,0))],
                  units='Angs',
                  multiplicity=2)
    en,orbe,orbs = dft(li,ETemp=1e4)
    return en

class LiFTDFTTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of Li/FT/DFT close to -7.349422?"""
        E = main()
        self.assertInside(E, energy, 1e-4)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(LiFTDFTTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())                                   
