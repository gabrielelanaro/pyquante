#!/usr/bin/env python
import unittest, sciunittest

from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

energy = -1.132473 # For etemp=1e4
name = 'H2/FT/DFT'

def main():
    r = 0.7
    h2 = Molecule('h2',
                  atomlist = [(1,(0,0,r/2.)),
                              (1,(0,0,-r/2.))],
                  units='Angs')
    en,orbe,orbs = dft(h2,ETemp=1e4)
    return en

class H2FTDFTTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of H2/FT/DFT close to -1.132473?"""
        E = main()
        self.assertInside(E, energy, 1e-3)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(H2FTDFTTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())
