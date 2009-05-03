#!/usr/bin/env python
"H2O using Gaussians"

import unittest, sciunittest

from PyQuante.MINDO3 import scf
from PyQuante.Molecule import Molecule

energy = -48.825159
name = "H2O MINDO/3"

def main():
    atomlist = Molecule('h2o',
                        atomlist = [(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))],
                        units = 'Angstrom')
    en = scf(atomlist)
    return en


def profmain():
    import cProfile,pstats
    cProfile.run('main()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

class WaterMindoTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of H2O (Mindo) close to -48.825159?"""
        E = main()
        self.assertInside(E, energy, 1e-7)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(WaterMindoTest)

if __name__ == '__main__':
    import unittest
    unittest.TextTestRunner(verbosity=2).run(suite())
