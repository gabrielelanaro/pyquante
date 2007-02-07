#!/usr/bin/env python
# Test for NO using UHF
# Jaguar ROHF: Rno = 1.125057 A = 2.12604 au; E=-129.24135247609
#        UHF:  Rno = 1.126917 A = 2.12955 au; E=-129.24788298992
# GAMESS-UK
# Energy -129.2478829288 (uhf at 2.12955 bohr)

import unittest, sciunittest

from PyQuante.hartree_fock import uhf
from PyQuante.Molecule import Molecule

name = "NO/UHF"
energy = -129.247769

def main():
    atomlist = Molecule('NO',atomlist = [(7,(0,0,0)),(8,(2.12955,0,0))],
                        multiplicity=2)
    en,orbe,orbs = uhf(atomlist)
    return en

class NOUHFTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of NO (UHF) close to -129.247769?"""
        result = main()
        self.assertInside(result, energy, 1e-4)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(NOUHFTest)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
