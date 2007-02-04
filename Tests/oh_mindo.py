#!/usr/bin/env python
"OH using Mindo"

import unittest, sciunittest

from PyQuante.MINDO3 import scf
from PyQuante.Molecule import Molecule

energy = 18.127533
name = "OH MINDO/3"

def main():
    atomlist = Molecule('oh',atomlist = [(8,(0,0,0)),(1,(1.,0,0))])
    en = scf(atomlist)
    return en

def profmain():
    import profile,pstats
    profile.run('main()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

class OHTest(sciunittest.TestCase):
    def runTest(self):
        """Energy of OH (using Mindo) close to 18.128?"""
        result = main()
        self.assertInside(result, energy, 1e-6)

def suite():
    return unittest.TestLoader().loadTestsFromTestCase(OHTest)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
