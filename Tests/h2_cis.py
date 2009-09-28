#!/usr/bin/env python
"H2 first excited state energy"

import unittest, sciunittest

from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule
from PyQuante.CI import CIS

energy = -0.645236
name = "H2_CIS"

def main():
    atoms = Molecule('h2',[(1,(1.,0,0)),(1,(-1.,0,0))])
    bfs = getbasis(atoms)
    S,h,Ints = getints(bfs,atoms)
    en,orbe,orbs = rhf(atoms,integrals=(S,h,Ints))
    occs = [1.]+[0.]*9

    Ecis = CIS(Ints,orbs,orbe,1,9,en)
    return Ecis[0]

class H2CISTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of H2 first excited state (CIS) close to -0.645236?"""
        E = main()
        self.assertInside(E, energy, 1e-4)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(H2CISTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())


    

    
