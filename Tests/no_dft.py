# Test for NO using spin-averaged DFT
# Jaguar ROHF: Rno = 1.125057 A = 2.12604 au; E=-129.24135247609
#        UHF:  Rno = 1.126917 A = 2.12955 au; E=-129.24788298992
# GAMESS-UK
# Energy -129.2478829288 (uhf at 2.12955 bohr)

import unittest, sciunittest

from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

name = "NO/DFT"
# Don't know whether this is right. This is what the energy
#  gave the first time I ran it.
energy = -128.897264

def main():
    atomlist = Molecule('NO',atomlist = [(7,(0,0,0)),(8,(2.12955,0,0))],
                        multiplicity=2)
    en,orbe,orbs = dft(atomlist)
    return en

class NODFTTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of NO (DFT) close to -128.897264?"""
        E = main()
        self.assertInside(E, energy, 1e-4)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(NODFTTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())                                   

