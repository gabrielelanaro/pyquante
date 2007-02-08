#!/usr/bin/env python
"""\
Rick's stab at unit testing

<beavis>heh, heh, he said *unit*</beavis>
"""

import unittest
from PyQuante.CI import CIS
from PyQuante.Ints import getbasis,getints
from PyQuante.Molecule import Molecule
from PyQuante.MP import MP2
from PyQuante.SCF import SCF,USCF

tol = 1e-5

class H2Test(unittest.TestCase):
    def setUp(self):
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],
                           units='Angs')

    def testHF(self):
        "HF energy of H2 close to -1.130501?"
        results = SCF(self.h2)
        self.assertAlmostEqual(results[0],-1.130501,4)

    def testSVWN(self):
        "SVWN energy of H2 close to -1.132710?"
        results = SCF(self.h2,functional='SVWN')
        self.assertAlmostEqual(results[0],-1.132710,4)

    def testBLYP(self):
        "BLYP energy of H2 close to -1.166221?"
        # This is slooooow.
        results = SCF(self.h2,functional='BLYP')
        self.assertAlmostEqual(results[0],-1.166221,4)

    def testHFFT(self):
        "HF energy of H2 at temp=1e4 is close to -1.130502?"
        results = SCF(self.h2,etemp=1e4)
        self.assertAlmostEqual(results[0],-1.130502,4)

    def testSVWNFT(self):
        "SVW energy of H2 at temp=1e4 is close to -1.132473?"
        results = SCF(self.h2,functional='SVWN',etemp=1e4)
        self.assertAlmostEqual(results[0],-1.132473,4)

    def testMP2(self):
        "MP2 energy of H2 is close to -1.156769"
        bfs = getbasis(self.h2)
        S,h,Ints = getints(bfs,self.h2)
        nbf = len(bfs)
        nclosed,nopen = self.h2.get_closedopen()
        en,orbe,orbs = SCF(self.h2,integrals=(S,h,Ints))
        emp2 = MP2(Ints,orbs,orbe,nclosed,nbf-nclosed)
        self.assertAlmostEqual(en+emp2,-1.156769,4)        

    def testCIS(self):
        "CIS energy of H2 is close to -0.559048"
        bfs = getbasis(self.h2)
        S,h,Ints = getints(bfs,self.h2)
        en,orbe,orbs = SCF(self.h2,integrals=(S,h,Ints))
        occs = [1.]+[0.]*9
        Ecis = CIS(Ints,orbs,orbe,occs,en)
        self.assertAlmostEqual(Ecis[0],-0.559048,4)
    
class LiTest(unittest.TestCase):
    def setUp(self):
        self.li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)
        self.li_p = Molecule('Li+',atomlist = [(3,(0,0,0))],charge=1)
        self.li_m = Molecule('Li-',atomlist = [(3,(0,0,0))],charge=-1)

    def testUHF(self):
        "UHF energy of Li is close to -7.431364"
        results = USCF(self.li)
        self.assertAlmostEqual(results[0],-7.431364,4)

    def testHFCation(self):
        "HF energy of Li+ is close to -7.235536"
        results = SCF(self.li_p)
        self.assertAlmostEqual(results[0],-7.235536,4)

    def testHFAnion(self):
        "HF energy of Li- is close to -7.407030"
        results = SCF(self.li_m)
        self.assertAlmostEqual(results[0],-7.407030,4)

    def testSVWN(self):
        "SVWN energy of Li is close to -7.332050"
        results = SCF(self.li,functional='SVWN')
        self.assertAlmostEqual(results[0],-7.332050,4)

    def testSVWNFT(self):
        "SVWN energy of Li at temp=1e4 is close to -7.349422"
        results = SCF(self.li,functional='SVWN',etemp=1e4)
        self.assertAlmostEqual(results[0],-7.349422,4)
        
if __name__ == '__main__':
    fullsuite = unittest.TestSuite()
    fullsuite.addTest(unittest.TestLoader().loadTestsFromTestCase(H2Test))
    fullsuite.addTest(unittest.TestLoader().loadTestsFromTestCase(LiTest))
    unittest.TextTestRunner(verbosity=2).run(fullsuite)

