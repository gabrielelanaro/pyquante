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

# Import test molecules
from PyQuante.TestMolecules import h2,he,li,li_p,li_m,h2o,oh,lih

class UnitTests(unittest.TestCase):
    def testH2BLYP(self):
        h2_blyp = SCF(h2,method="DFT",functional='BLYP')
        h2_blyp.iterate()
        self.assertAlmostEqual(h2_blyp.energy,-1.167767,4)

    def testH2LDA(self):
        h2_lda = SCF(h2,method='DFT',functional="SVWN")
        h2_lda.iterate()
        self.assertAlmostEqual(h2_lda.energy,-1.135061,4)

    def testLiLDA(self):
        li_lda = SCF(li,method='DFT',functional="SVWN")
        li_lda.iterate()
        self.assertAlmostEqual(li_lda.energy,-7.332050,4)

    def testLiUHF(self):
        li_uhf = SCF(li,method='UHF')
        li_uhf.iterate()
        self.assertAlmostEqual(li_uhf.energy,-7.431364,4)

    def testLiROHF(self):
        li_uhf = SCF(li,method='ROHF')
        li_uhf.iterate()
        self.assertAlmostEqual(li_uhf.energy,-7.431369,4)

    def testLiUHFFT(self):
        li_uhf = SCF(li,method="UHF",etemp=1e4)
        li_uhf.iterate()
        # No test, since I don't really know what the energy should be:
        # finite temperature HF is kind of a hack. But this at least
        # tests that the program runs
        return

    ########## Solver tests ##########

    def testSubspaceSolver(self):
        solv = SCF(h2,method='HF',SolverConstructor=SubspaceSolver)
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)
    
    def testDavidsonSolver(self):
        from PyQuante.Solvers import davidson
        solv = SCF(h2,method='HF',SolverConstructor=SubspaceSolver,
                    solver=davidson,pass_nroots=True)
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)

    def testJacobiSolver(self):
        from PyQuante.Solvers import jacobi
        solv = SCF(h2,method='HF',SolverConstructor=SubspaceSolver,
                    solver=jacobi)
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)
    
    def testTCPSolver(self):
        solv = SCF(h2,method='HF',SolverConstructor=DmatSolver)
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)

    def testNOTCPSolver(self):
        from PyQuante.DMP import NOTCP
        solv = SCF(h2,method='HF',SolverConstructor=DmatSolver,
                    solver=NOTCP)
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)

    def testTRPSolver(self):
        from PyQuante.DMP import TRP
        solv = SCF(h2,method='HF',SolverConstructor=DmatSolver,
                    solver=TRP)
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)

    def testCPSolver(self):
        from PyQuante.DMP import CP
        solv = SCF(h2,method='HF',SolverConstructor=DmatSolver,
                    solver=CP)
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)
    
    def testMCWSolver(self):
        from PyQuante.DMP import McWeeny
        solv = SCF(h2,method='HF',SolverConstructor=DmatSolver,
                    solver=McWeeny)
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)

    ########## Basis set tests ##########

    def testSTO3G(self):
        solv = SCF(h2,method='HF',basis='sto-3g')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.117099,4)

    def testSTO6G(self):
        solv = SCF(h2, method="HF",basis='sto-3g')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.117099,4)

    def test321G(self):
        solv = SCF(h2, method="HF",basis='3-21g')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.122956,4)

    def test631Gss(self):
        solv = SCF(h2, method="HF",basis='6-31g**')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131334,4)

    def test631Gppss(self):
        solv = SCF(h2, method="HF",basis='6-31g++**')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131403,4)

    def test631Gdp(self):
        solv = SCF(h2, method="HF",basis='6-31G(d,p)')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.131344,4)

    def testVDZ(self):
        solv = SCF(h2, method="HF",basis='cc-pvdz')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.128571,4)

    def testVTZ(self):
        solv = SCF(h2, method="HF",basis='cc-pvtz')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.133009,4)

    def testDZVP(self):
        solv = SCF(h2, method="HF",basis='dzvp')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.127306,4)

    def test6311G(self):
        solv = SCF(h2, method="HF",basis='6-311G**')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.132485,4)

    def test6311Gdp(self):
        solv = SCF(h2, method="HF",basis='6-311G++(2d,2p)')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.133063,4)

    def test6311G3d3p(self):
        solv = SCF(h2, method="HF",basis='6-311G++(3d,3p)')
        solv.iterate()
        self.assertAlmostEqual(solv.energy,-1.133023,4)

    ########## MINDO3 Tests ##########

    def testH2OMINDO(self):
        h2o_mindo3 = SCF(h2o,method="MINDO3")
        h2o_mindo3.iterate()
        self.assertAlmostEqual(h2o_mindo3.energy,-53.5176,2)

    def testOHMINDO(self):
        oh_mindo = SCF(oh,method="UMINDO3")
        oh_mindo.iterate()
        self.assertAlmostEqual(oh_mindo.energy,16.49043,2)

    ########## Misc Tests ##########

    def testOrthog(self):
        from PyQuante.LA2 import CanOrth, SymOrth, CholOrth, simx
        from PyQuante.NumWrap import eigh
        solver = SCF(h2,method="HF")
        h,S = solver.h, solver.S
        X1 = CanOrth(S)
        X2 = SymOrth(S)
        X3 = CholOrth(S)
        h1 = simx(h,X1)
        ha = simx(h,X2)
        h3 = simx(h,X3)
        e1,v1 = eigh(h1)
        e2,v2 = eigh(ha)
        e3,v3 = eigh(h3)
        self.assertAlmostEqual(e1[0],e2[0],6)
        self.assertAlmostEqual(e1[0],e3[0],6)

    def testMP2(self):
        solv = SCF(h2,method="HF")
        solv.iterate()
        nclosed,nopen = h2.get_closedopen()
        nbf = len(solv.basis_set.get())
        emp2 = MP2(solv.ERI,solv.solver.orbs,solv.solver.orbe,nclosed,nbf-nclosed)
        self.assertAlmostEqual(solv.energy+emp2,-1.157660,4)        

    def testCIS(self):
        solv = SCF(h2,method="HF")
        solv.iterate()
        nclosed,nopen = h2.get_closedopen()
        nbf = len(solv.basis_set.get())
        nocc = nclosed+nopen
        nvirt = nbf-nocc
        Ecis = CIS(solv.ERI,solv.solver.orbs,solv.solver.orbe,nocc,
                   nvirt,solv.energy)
        self.assertAlmostEqual(Ecis[0],-0.573134,3)

    def testLiH_OEP_AN(self):
        do_oep_an = True
        lih_hf = SCF(lih,method="HF")
        lih_hf.iterate()
        ints = lih_hf.S,lih_hf.h,lih_hf.ERI
        E_exx,orbe_exx,orbs_exx = oep_hf_an(lih,lih_hf.solver.orbs,
                                            bfs=lih_hf.basis_set.get(),
                                            integrals=ints)
        self.assertAlmostEqual(E_exx,-7.981282,4)

def runsuite(verbose=True):
    # To use psyco, uncomment this line:
    #import psyco; psyco.full()
    if verbose: verbosity=2
    else: verbosity=1
    # If you want more output, uncomment this line:
    #logging.basicConfig(format="%(message)s",level=logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(UnitTests)
    unittest.TextTestRunner(verbosity=verbosity).run(suite)
    # Running without verbosity is equivalent to replacing the above
    # two lines with the following:
    #unittest.main()
    return

def debugsuite():
    import cProfile,pstats
    cProfile.run('runsuite()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__':
    import sys
    if "-d" in sys.argv:
        debugsuite()
    else:
        runsuite()
    
