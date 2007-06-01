#!/usr/bin/env python
"""\
UnitSweet.py - Unit testing for Python.

<beavis>heh, heh, he said *unit*</beavis>
"""

import unittest,logging
from PyQuante.CI import CIS
from PyQuante.Molecule import Molecule
from PyQuante.MP import MP2
from PyQuante.OEP import oep_hf,oep_hf_an
from PyQuante.PyQuante2 import SCF,SubspaceSolver,DmatSolver

class UnitTests(unittest.TestCase):
    def setUp(self):
        from PyQuante.Molecule import Molecule
        self.h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],
                           units='Angs')
        self.he = Molecule('He',atomlist = [(2,(0,0,0))])
        self.li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)
        self.li_p = Molecule('Li+',atomlist = [(3,(0,0,0))],charge=1)
        self.li_m = Molecule('Li-',atomlist = [(3,(0,0,0))],charge=-1)
        self.h2o = Molecule('h2o',[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))],
                            units="Angstrom")
        self.oh = Molecule('oh',[(8,(0,0,0)),(1,(1.,0,0))],
                            units="Angstrom")
        self.lih = Molecule('LiH',[(1,(0,0,1.5)),(3,(0,0,-1.5))],units='Bohr')

    def testH2HF(self):
        h2_hf = SCF(self.h2,method='HF')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testH2HFFT(self):
        h2_hf = SCF(self.h2,method='HF',etemp=1e4)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130502,4)

    def testH2BLYP(self):
        h2_blyp = SCF(self.h2,method="DFT",functional='BLYP')
        h2_blyp.iterate()
        self.assertAlmostEqual(h2_blyp.energy,-1.166286,4)

    def testHeHF(self):
        he_hf = SCF(self.he,method='HF')
        he_hf.iterate()
        self.assertAlmostEqual(he_hf.energy,-2.855260,3)

    def testLiHF(self):
        li_hf = SCF(self.li,method='HF')
        li_hf.iterate()
        return

    def testLipHF(self):
        li_p_hf = SCF(self.li_p,method='HF')
        li_p_hf.iterate()
        self.assertAlmostEqual(li_p_hf.energy,-7.235536,4)

    def testLimHF(self):
        li_m_hf = SCF(self.li_m,method='HF')
        li_m_hf.iterate()
        self.assertAlmostEqual(li_m_hf.energy,-7.407030,4)

    def testH2LDA(self):
        h2_lda = SCF(self.h2,method='DFT',functional="SVWN")
        h2_lda.iterate()
        self.assertAlmostEqual(h2_lda.energy,-1.132799,4)

    def testH2LDAFT(self):
        h2_lda = SCF(self.h2,method='DFT',functional="SVWN",etemp=1e4)
        h2_lda.iterate()
        self.assertAlmostEqual(h2_lda.energy,-1.132558,4)

    def testLiLDA(self):
        li_lda = SCF(self.li,method='DFT',functional="SVWN")
        li_lda.iterate()
        self.assertAlmostEqual(li_lda.energy,-7.332050,4)

    def testLiLDAFT(self):
        li_lda = SCF(self.li,method='DFT',functional="SVWN",etemp=1e4)
        li_lda.iterate()
        self.assertAlmostEqual(li_lda.energy,-7.349422,4)

    def testLiUHF(self):
        li_uhf = SCF(self.li,method='UHF')
        li_uhf.iterate()
        self.assertAlmostEqual(li_uhf.energy,-7.431364,4)

    def testLiUHFFT(self):
        li_uhf = SCF(self.li,method="UHF",etemp=1e4)
        li_uhf.iterate()
        # No test, since I don't really know what the energy should be:
        # finite temperature HF is kind of a hack. But this at least
        # tests that the program runs
        return

    ########## Solver tests ##########

    def testSubspaceSolver(self):
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=SubspaceSolver)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)
    
    def testDavidsonSolver(self):
        from PyQuante.Solvers import davidson
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=SubspaceSolver,
                    solver=davidson,pass_nroots=True)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testJacobiSolver(self):
        from PyQuante.Solvers import jacobi
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=SubspaceSolver,
                    solver=jacobi)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)
    
    def testTCPSolver(self):
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=DmatSolver)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testTRPSolver(self):
        from PyQuante.DMP import TRP
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=DmatSolver,
                    solver=TRP)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testCPSolver(self):
        from PyQuante.DMP import CP
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=DmatSolver,
                    solver=CP)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)
    
    def testMCWSolver(self):
        from PyQuante.DMP import McWeeny
        h2_hf = SCF(self.h2,method='HF',SolverConstructor=DmatSolver,
                    solver=McWeeny)
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    ########## Basis set tests ##########

    def testSTO3G(self):
        h2_hf = SCF(self.h2,method='HF',basis='sto-3g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.117349,4)

    def testSTO6G(self):
        h2_hf = SCF(self.h2, method="HF",basis='sto-3g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.117349,4)

    def test321G(self):
        h2_hf = SCF(self.h2, method="HF",basis='3-21g')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.121998,4)

    def test631Gss(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-31g**')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def test631Gppss(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-31g++**')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130553,4)

    def test631Gdp(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-31G(d,p)')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.130501,4)

    def testVDZ(self):
        h2_hf = SCF(self.h2, method="HF",basis='cc-pvdz')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.126923,4)

    def testVTZ(self):
        h2_hf = SCF(self.h2, method="HF",basis='cc-pvtz')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.132136,4)

    def testDZVP(self):
        h2_hf = SCF(self.h2, method="HF",basis='dzvp')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.126728,4)

    def test6311G(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-311G**')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.131516,4)

    def test6311Gdp(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-311G++(2d,2p)')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.132122,4)

    def test6311G3d3p(self):
        h2_hf = SCF(self.h2, method="HF",basis='6-311G++(3d,3p)')
        h2_hf.iterate()
        self.assertAlmostEqual(h2_hf.energy,-1.132166,4)

    ########## MINDO3 Tests ##########

    def testH2OMINDO(self):
        h2o_mindo3 = SCF(self.h2o,method="MINDO3")
        h2o_mindo3.iterate()
        self.assertAlmostEqual(h2o_mindo3.energy,-48.826208,2)

    def testOHMINDO(self):
        oh_mindo = SCF(self.oh,method="UMINDO3")
        oh_mindo.iterate()
        self.assertAlmostEqual(oh_mindo.energy,18.1258,2)

    ########## Misc Tests ##########

    def testOrthog(self):
        from PyQuante.LA2 import CanOrth, SymOrth, CholOrth, simx
        from PyQuante.NumWrap import eigh
        solver = SCF(self.h2,method="HF")
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

    def testMP2(self):
        h2_hf = SCF(self.h2,method="HF")
        h2_hf.iterate()
        nclosed,nopen = self.h2.get_closedopen()
        nbf = len(h2_hf.basis_set.get())
        emp2 = MP2(h2_hf.ERI,h2_hf.solver.orbs,h2_hf.solver.orbe,nclosed,nbf-nclosed)
        self.assertAlmostEqual(h2_hf.energy+emp2,-1.156769,4)        

    def testCIS(self):
        h2_hf = SCF(self.h2,method="HF")
        h2_hf.iterate()
        occs = [1.]+[0.]*9
        Ecis = CIS(h2_hf.ERI,h2_hf.solver.orbs,h2_hf.solver.orbe,occs,h2_hf.energy)
        self.assertAlmostEqual(Ecis[0],-0.559115,4)

    def testLiH_OEP_AN(self):
        do_oep_an = True
        lih_hf = SCF(self.lih,method="HF")
        lih_hf.iterate()
        ints = lih_hf.S,lih_hf.h,lih_hf.ERI
        E_exx,orbe_exx,orbs_exx = oep_hf_an(self.lih,lih_hf.solver.orbs,
                                            bfs=lih_hf.basis_set.get(),
                                            integrals=ints)
        self.assertAlmostEqual(E_exx,-7.981044,4)

if __name__ == '__main__':
    #suite = unittest.TestLoader().loadTestsFromTestCase(UnitTests)
    #unittest.TextTestRunner(verbosity=2).run(suite)
    #logging.basicConfig(format="%(message)s",level=logging.DEBUG)
    unittest.main()
