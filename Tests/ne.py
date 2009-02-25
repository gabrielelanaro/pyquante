#!/usr/bin/env python
"Neon using Gaussians"

import unittest, sciunittest

from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule

# GAMESS-UK HF Energy
# Energy -128.4744065199

energy = -128.474406 # Changed 2003-04-07 to reflect diis
name = "Ne"

def main():
    ne = Molecule('Ne',atomlist = [(10,(0,0,0))])
    en,orbe,orbs = rhf(ne)
    return en

class NeTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of Ne (using Gaussians) close to -128.474406?"""
        result = main()
        self.assertInside(result, energy, 1e-4)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(NeTest)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
