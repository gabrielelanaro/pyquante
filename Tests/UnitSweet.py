#!/usr/bin/env python
"""\
UnitSweet.py - Unit testing for Python.

<beavis>heh, heh, he said *unit*</beavis>
"""

import unittest
from PyQuante.CI import CIS
from PyQuante.Molecule import Molecule
from PyQuante.MP import MP2
from PyQuante.HFSolver import HFSolver
from PyQuante.UHFSolver import UHFSolver
from PyQuante.DFTSolver import DFTSolver
from PyQuante.OEP import oep_hf,oep_hf_an
from PyQuante.MINDO3 import scf

class HFUnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],units='Angs')
        self.he = Molecule('He',atomlist = [(2,(0,0,0))])
        self.li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)
        self.li_p = Molecule('Li+',atomlist = [(3,(0,0,0))],charge=1)
        self.li_m = Molecule('Li-',atomlist = [(3,(0,0,0))],charge=-1)

    def testH2HF(self):
        h2_hf = HFSolver(self.h2)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testH2HFFT(self):
        h2_hf = HFSolver(self.h2)
        h2_hf.iterate(etemp=1e4)
        self.assertAlmostEqual(h2_hf.energy,-1.130502,4)

    def testHeHF(self):
        he_hf = HFSolver(self.he)
        he_hf.iterate()
        self.assertAlmostEqual(he_hf.energy,-2.855260,4)

    def testLiHF(self):
        li_hf = HFSolver(self.li)
        li_hf.iterate()
        return

    def testLipHF(self):
        li_p_hf = HFSolver(self.li_p)
        li_p_hf.iterate()
        self.assertAlmostEqual(li_p_hf.energy,-7.235536,4)

    def testLimHF(self):
        li_m_hf = HFSolver(self.li_m)
        li_m_hf.iterate()
        self.assertAlmostEqual(li_m_hf.energy,-7.407030,4)

class DFTUnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],units='Angs')
        self.he = Molecule('He',atomlist = [(2,(0,0,0))])
        self.li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)

    def testLiLDA(self):
        li_lda = DFTSolver(self.li)
        li_lda.iterate()
        self.assertAlmostEqual(li_lda.energy,-7.332050,4)

    def testLiLDAFT(self):
        li_lda = DFTSolver(self.li)
        li_lda.iterate(etemp=1e4)
        self.assertAlmostEqual(li_lda.energy,-7.349422,4)

    def testH2LDA(self):
        h2_lda = DFTSolver(self.h2)
        h2_lda.iterate()
        self.assertAlmostEqual(h2_lda.energy,-1.132710,4)

    def testH2LDAFT(self):
        h2_lda = DFTSolver(self.h2)
        h2_lda.iterate(etemp=1e4)
        self.assertAlmostEqual(h2_lda.energy,-1.132473,4)

    def testH2BLYP(self):
        h2_blyp = DFTSolver(self.h2,functional='BLYP')
        h2_blyp.iterate()
        self.assertAlmostEqual(h2_blyp.energy,-1.166221,4)

class UHFUnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)

    def testLiUHF(self):
        li_uhf = UHFSolver(self.li)
        li_uhf.iterate()
        self.assertAlmostEqual(li_uhf.energy,-7.431364,4)

    def testLiUHFFT(self):
        li_uhf = UHFSolver(self.li)
        li_uhf.iterate(etemp=1e4)
        return

class OtherUnitTests(unittest.TestCase):
    def setUp(self):
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],
                           units='Angs')
        self.lih = Molecule('LiH',[(1,(0,0,1.5)),(3,(0,0,-1.5))],units = 'Bohr')
        self.h2o = Molecule('h2o',[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))])
        self.ohm = Molecule('oh',[(8,(0,0,0)),(1,(1.,0,0))])

    def testH2OMINDO(self):
        en = scf(self.h2o)
        self.assertAlmostEqual(en,-48.825159,4)

    def testOHMMINDO(self):
        en = scf(self.ohm)
        self.assertAlmostEqual(en,18.127533,4)

    def testMP2(self):
        h2_hf = HFSolver(self.h2)
        h2_hf.iterate()
        nclosed,nopen = self.h2.get_closedopen()
        nbf = len(h2_hf.bfs)
        emp2 = MP2(h2_hf.Ints,h2_hf.orbs,h2_hf.orbe,nclosed,nbf-nclosed)
        self.assertAlmostEqual(h2_hf.energy+emp2,-1.156769,4)        

    def testCIS(self):
        h2_hf = HFSolver(self.h2)
        h2_hf.iterate()
        occs = [1.]+[0.]*9
        Ecis = CIS(h2_hf.Ints,h2_hf.orbs,h2_hf.orbe,occs,h2_hf.energy)
        self.assertAlmostEqual(Ecis[0],-0.559048,4)

    def testLiH_OEP_AN(self):
        do_oep_an = True
        lih_hf = HFSolver(self.lih)
        lih_hf.iterate()
        ints = lih_hf.S,lih_hf.h,lih_hf.Ints
        E_exx,orbe_exx,orbs_exx = oep_hf_an(self.lih,lih_hf.orbs,
                                            bfs=lih_hf.bfs,
                                            integrals=ints)
        self.assertAlmostEqual(E_exx,-7.981044,4)

    # Don't think it necessary to test both OEP and OEP_AN...
    #def testLiH_OEP(self):
    #    do_oep_an = True
    #    lih_hf = HFSolver(self.lih)
    #    lih_hf.iterate()
    #    ints = lih_hf.S,lih_hf.h,lih_hf.Ints
    #    E_exx,orbe_exx,orbs_exx = oep_hf(self.lih,lih_hf.orbs,
    #                                     bfs=lih_hf.bfs,
    #                                     integrals=ints)
    #    self.assertAlmostEqual(E_exx,-7.981044,4)
    
def suite():
    import ricks_unit
    return unittest.TestLoader().loadTestsFromModule(ricks_unit)
        
if __name__ == '__main__':
    fullsuite = suite()
    unittest.TextTestRunner(verbosity=2).run(fullsuite)
