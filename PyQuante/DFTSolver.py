import unittest
from PyQuante.HFSolver import HFSolver
class DFTSolver(HFSolver):
    """\
    solver = DFTSolver(molecule,**options)
    
    Create a solver that can perform a DFT calculation on *molecule*.

    General options

    Option        Value   Description
    --------      -----   -----------
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    basis         None    The name of a basis set, e.g. '6-31g**',
                          'sto-3g','cc-pVTZ'
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not none, the guess orbitals
    functional    None    Use HF for the calculation (default)
                  SVWN    Use the SVWN (LDA) DFT functional
                  S0      Use the Slater Xalpha DFT functional
                  BLYP    Use the BLYP GGA DFT functional
                  PBE     Use the PBE DFT functional
    grid_nrad     32      Number of radial shells per atom
    grid_fineness 1       Radial shell fineness. 0->coarse, 1->medium, 2->fine

    Options passed into solver.iterate(**options):

    Options:      Value   Description
    --------      -----   -----------
    etol          1e-5    Energy convergence criteria
    max_iter      20      Maximum SCF iterations
    do_averaging  True    Use DIIS for accelerated convergence (default)
                  False   No convergence acceleration
    etemp         False   Use etemp value for finite temperature DFT (default)
                  float   Use (float) for the electron temperature
    """
    def setup(self,**opts):
        from PyQuante.DFunctionals import need_gradients
        self.nel = self.molecule.get_nel()
        self.nclosed,self.nopen = self.molecule.get_closedopen()
        self.Enuke = self.molecule.get_enuke()
        self.energy = 0
        self.etemp=0 # This may be overwritten in iterate
        self.setup_basis(**opts)
        self.setup_integrals(**opts)
        self.setup_guess(**opts)
        self.functional = opts.get('functional','SVWN')
        self.method = self.functional
        opts['do_grad_dens'] = need_gradients[self.functional]
        self.setup_grid(**opts)
        return

    def setup_grid(self,**opts):
        from PyQuante.MolecularGrid import MolecularGrid
        self.grid_nrad = opts.get('grid_nrad',32)
        self.grid_fineness = opts.get('grid_fineness',1)
        self.gr = MolecularGrid(self.molecule,self.grid_nrad,
                                self.grid_fineness,**opts) 
        self.gr.set_bf_amps(self.bfs)
        self.bfgrid = self.gr.allbfs() # bfs over all grid points
        return

    def update_density(self):
        HFSolver.update_density(self)
        self.gr.setdens(self.D)

    def update_XC(self):
        from PyQuante.dft import getXC
        self.Exc,self.XC = getXC(self.gr,self.nel,self.bfgrid,
                                 functional=self.functional)
        return
    
    def update_fock(self):
        self.update_J()
        self.update_XC()
        self.F = self.h+2*self.J+self.XC
        return

class UnitTests(unittest.TestCase):
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
        
def test_suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(UnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)

def test():
    from PyQuante import Molecule
    al = Molecule('Al',atomlist = [(13,(0,0,0))],multiplicity=2)
    al_lda = DFTSolver(al)
    al_lda.iterate()
    
if __name__ == '__main__': test()
