#!/usr/bin/env python
"H2 correlation energy via MP2"

import unittest, sciunittest

from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule
from PyQuante.MP import MP2

# GAMESS-UK SCF ENERGY -1.0882670560
# GAMESS-UK MP2 ENERGY -1.11591790

energy = -1.1153708         #Original value
energy = -1.11591861694     #After increasing accuracy of HF convergence
name = "H2_MP2"

def main():
    atoms = Molecule('h2',[(1,(1.,0,0)),(1,(-1.,0,0))])
    bfs = getbasis(atoms)
    nel = atoms.get_nel()
    nbf = len(bfs)
    nocc = nel/2
    S,h,Ints = getints(bfs,atoms)
    en,orbe,orbs = rhf(atoms,integrals=(S,h,Ints),)

    emp2 = MP2(Ints,orbs,orbe,nocc,nbf-nocc)  
    return en+emp2

class H2MP2Test(sciunittest.TestCase):
    def runTest(self):
        """Energy of H2 (MP2) close to -1.111591861694?"""
        E = main()
        self.assertInside(E, energy, 1e-6)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(H2MP2Test)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite()) 
