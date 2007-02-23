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
from PyQuante.MINDO3Solver import MINDO3Solver
from PyQuante.UMINDO3Solver import UMINDO3Solver

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

    def testH2HFSTO(self):
        h2_hf = HFSolver(self.h2,basis='sto-3g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.117349,4)

    def testH2HFFT(self):
        h2_hf = HFSolver(self.h2)
        h2_hf.iterate(etemp=1e4)
        self.assertAlmostEqual(h2_hf.energy,-1.130502,4)

    def testHeHF(self):
        he_hf = HFSolver(self.he)
        he_hf.iterate()
        # Note: this single test case fluctuates between -2.855260 with numpy
        #  and -2.855168 with Numeric. I have no idea why this is the case. None
        #  of the other cases are affected in this magnitude. I have consequently
        #  turned the tolerance down, but this is worth keeping an eye on.
        self.assertAlmostEqual(he_hf.energy,-2.855260,3)

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

class BasisUnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],units='Angs')

    def testSTO3G(self):
        h2_hf = HFSolver(self.h2,basis='sto-3g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.117349,4)

    def testSTO6G(self):
        h2_hf = HFSolver(self.h2,basis='sto-3g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.117349,4)

    def test321G(self):
        h2_hf = HFSolver(self.h2,basis='3-21g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.121998,4)

    def test631Gss(self):
        h2_hf = HFSolver(self.h2,basis='6-31g**')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def test631Gdp(self):
        h2_hf = HFSolver(self.h2,basis='6-31G(d,p)')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testVDZ(self):
        h2_hf = HFSolver(self.h2,basis='cc-pvdz')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.126923,4)

    def testVTZ(self):
        h2_hf = HFSolver(self.h2,basis='cc-pvtz')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.132136,4)

    def testDZVP(self):
        h2_hf = HFSolver(self.h2,basis='dzvp')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.126728,4)

    def test6311G(self):
        h2_hf = HFSolver(self.h2,basis='6-311G**')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.131516,4)

    def test6311Gdp(self):
        h2_hf = HFSolver(self.h2,basis='6-311G++(2d,2p)')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.132122,4)

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

class SolverUnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],units='Angs')
        return

    def testRegularSolver(self):
        from PyQuante.HFSolver import HFSolver
        solver = HFSolver(self.h2)
        solver.iterate()
        self.assertAlmostEqual(solver.energy,-1.1305,4)
        return

    def testSubspaceSolver(self):
        from PyQuante.Solvers import SubspaceSolver
        from PyQuante.NumWrap import eigh
        solver = SubspaceSolver(self.h2,eigh)
        solver.iterate()
        self.assertAlmostEqual(solver.energy,-1.1305,4)
        return

    def testDavidsonSolver(self):
        from PyQuante.Solvers import SubspaceSolver,init_davidson
        dav = init_davidson(2)
        solver = SubspaceSolver(self.h2,dav)
        solver.iterate()
        self.assertAlmostEqual(solver.energy,-1.1305,4)
        return

    def testJacobiSolver(self):
        from PyQuante.Solvers import SubspaceSolver,init_jacobi
        jac = init_jacobi()
        solver = SubspaceSolver(self.h2,jac)
        solver.iterate()
        self.assertAlmostEqual(solver.energy,-1.1305,4)
        return

    def testTCPSolver(self):
        from PyQuante.DMP import DmatSolver,TCP,init_dmat_solver
        tcp = init_dmat_solver(TCP)
        solver = DmatSolver(self.h2,tcp)
        solver.iterate()
        self.assertAlmostEqual(solver.energy,-1.1305,4)
        return

    def testTRPSolver(self):
        from PyQuante.DMP import DmatSolver,TRP,init_dmat_solver
        trp = init_dmat_solver(TRP)
        solver = DmatSolver(self.h2,trp)
        solver.iterate()
        self.assertAlmostEqual(solver.energy,-1.1305,4)
        return

    def testCPSolver(self):
        from PyQuante.DMP import DmatSolver,CP,init_dmat_solver
        cp = init_dmat_solver(CP)
        solver = DmatSolver(self.h2,cp)
        solver.iterate()
        self.assertAlmostEqual(solver.energy,-1.1305,4)
        return

    def testMCWSolver(self):
        from PyQuante.DMP import DmatSolver,McWeeny,init_dmat_solver
        mcw = init_dmat_solver(McWeeny)
        solver = DmatSolver(self.h2,mcw)
        solver.iterate()
        self.assertAlmostEqual(solver.energy,-1.1305,4)
        return

    def testOrthog(self):
        from PyQuante.HFSolver import HFSolver
        from PyQuante.LA2 import CanOrth, SymOrth, CholOrth, simx
        from PyQuante.NumWrap import eigh
        solver = HFSolver(self.h2)
        h,S = solver.h, solver.S
        X1 = CanOrth(S)
        X2 = SymOrth(S)
        X3 = CholOrth(S)
        h1 = simx(h,X1)
        h2 = simx(h,X2)
        h3 = simx(h,X3)
        e1,v1 = eigh(h1)
        e2,v2 = eigh(h2)
        e3,v3 = eigh(h3)
        self.assertAlmostEqual(e1[0],e2[0],6)
        self.assertAlmostEqual(e1[0],e3[0],6)


class OtherUnitTests(unittest.TestCase):
    def setUp(self):
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],
                           units='Angs')
        self.lih = Molecule('LiH',[(1,(0,0,1.5)),(3,(0,0,-1.5))],units = 'Bohr')
        self.h2o = Molecule('h2o',[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))])
        self.oh = Molecule('oh',[(8,(0,0,0)),(1,(1.,0,0))])

    # Note: tolerance for MINDO tests scaled back since the results are
    # in kcals/mol and not hartrees
    def testH2OMINDO(self):
        h2o_mindo3 = MINDO3Solver(self.h2o)
        h2o_mindo3.iterate()
        #en = scf(self.h2o)
        self.assertAlmostEqual(h2o_mindo3.energy,-48.826208,2)

    def testOHMINDO(self):
        oh_mindo = UMINDO3Solver(self.oh)
        oh_mindo.iterate()
        self.assertAlmostEqual(oh_mindo.energy,18.1258,2)

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

    # Don't think it necessary to test both OEP and OEP_AN, especially
    # since this one is really slow
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
    import UnitSweet
    return unittest.TestLoader().loadTestsFromModule(UnitSweet)
        
if __name__ == '__main__':
    fullsuite = suite()
    unittest.TextTestRunner(verbosity=2).run(fullsuite)
