import unittest
from PyQuante.AbstractSolver import AbstractSolver
from PyQuante import logging

class HFSolver(AbstractSolver):
    """\
    solver = HFSolver(molecule,**options)
    
    Create a solver that can perform a HF calculation on *molecule*.

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

    Options passed into solver.iterate(**options):

    Options:      Value   Description
    --------      -----   -----------
    etol          1e-4    Energy convergence criteria
    max_iter      20      Maximum SCF iterations
    do_averaging  True    Use DIIS for accelerated convergence (default)
                  False   No convergence acceleration
    etemp         False   Use etemp value for finite temperature DFT (default)
                  float   Use (float) for the electron temperature
    """
    def setup_iterations(self,**opts):
        AbstractSolver.setup_iterations(self,**opts)
        self.do_averaging = opts.get('do_averaging',True)
        if self.do_averaging: self.setup_averaging()
        return

    def setup_averaging(self,**opts):
        from PyQuante.Convergence import DIIS
        self.avg = DIIS(self.S)
        return        

    def print_setup_info(self):
        AbstractSolver.print_setup_info(self)
        logging.debug("Nbf = %d" % len(self.bfs))
        return

    def print_iteration_info(self):
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
                  (self.iter,self.energy,self.Eone,self.Ej,self.Exc,self.Enuke))
        return

    def print_pre_iteration_info(self):
        if self.etemp:
            logging.info("Electron temperature = %.0f" % self.etemp)
        if self.do_averaging:
            logging.debug("Using DIIS averaging")
        AbstractSolver.print_pre_iteration_info(self)
        return
    
    def setup(self,**opts):
        AbstractSolver.setup(self,**opts)
        self.method = "HF"
        self.etemp=0 # This may be overwritten in iterate
        self.setup_basis(**opts)
        self.setup_integrals(**opts)
        self.setup_guess(**opts)
        self.update_density()
        return

    def setup_basis(self,**opts):
        from PyQuante.Ints import getbasis
        from PyQuante.Basis.Tools import get_basis_data
        self.bfs = opts.get('bfs',None)
        if not self.bfs:
            basis_data = opts.get('basis_data',None)
            if not basis_data:
                basis = opts.get('basis',None)
                if basis:
                    basis_data = get_basis_data(basis)
            self.bfs = getbasis(self.molecule,basis_data)
        return

    def setup_integrals(self,**opts):
        from PyQuante.Ints import getints
        integrals = opts.get('integrals',None)
        if integrals:
            self.S,self.h,self.Ints = integrals
        else:
            self.S,self.h,self.Ints = getints(self.bfs,self.molecule)
        return

    def setup_guess(self,**opts):
        from PyQuante.LA2 import geigh
        self.orbs = opts.get('orbs',None)
        if self.orbs is None:
            self.orbe,self.orbs = geigh(self.h,self.S)
        return

    def update_density(self):
        from PyQuante.LA2 import mkdens_spinavg
        from PyQuante.fermi_dirac import mkdens_fermi
        if self.etemp:
            self.D,self.entropy = mkdens_fermi(self.nel,self.orbe,self.orbs,self.etemp)
        else:
            self.D = mkdens_spinavg(self.orbs,self.nclosed,self.nopen)
            self.entropy=0
        return

    def update_J(self):
        from PyQuante.LA2 import trace2
        from PyQuante.Ints import getJ
        self.J = getJ(self.Ints,self.D)
        self.Ej = 2*trace2(self.D,self.J)
        return

    def update_K(self):
        from PyQuante.Ints import getK
        from PyQuante.LA2 import trace2
        self.K = getK(self.Ints,self.D)
        self.Exc = -trace2(self.D,self.K)
        return
    
    def update_fock(self):
        self.update_J()
        self.update_K()
        self.F = self.h + 2*self.J-self.K
        if self.do_averaging:
            self.F = self.avg.getF(self.F,self.D)
        return

    def solve_fock(self):
        from PyQuante.LA2 import geigh
        self.orbe,self.orbs = geigh(self.F,self.S)
        self.update_density()
        return

    def calculate_energy(self):
        from PyQuante.LA2 import trace2
        self.Eone = 2*trace2(self.D,self.h)
        self.energy = self.Eone + self.Ej + self.Exc + self.Enuke + self.entropy
        return        

    def is_converged(self):
        if not self.energy_history:
            self.energy_history.append(self.energy)
            return False
        is_converged_value = abs(self.energy-self.energy_history[-1]) < self.etol
        self.energy_history.append(self.energy)
        return is_converged_value

class UnitTests(unittest.TestCase):
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

def test():
    suite = unittest.TestLoader().loadTestsFromTestCase(UnitTests)
    unittest.TextTestRunner(verbosity=2).run(suite)
    
if __name__ == '__main__': test()
