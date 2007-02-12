from PyQuante import logging
from PyQuante.AbstractSolver import AbstractSolver

class MINDO3Solver(AbstractSolver):
    def setup(self):
        from PyQuante.MINDO3 import initialize, get_nbf, get_reference_energy,\
             get_F0, get_nel,get_open_closed,get_enuke
        self.charge = self.molecule.charge
        self.multiplicity = self.molecule.multiplicity
        # This is an ugly-ish hack to deal with the brain-dead
        #  way that MINDO3.get_open_closed works
        if self.multiplicity == 1: self.multiplicity=None
        # Ultimately I should subclass Atom for MINDO3Atom
        self.molecule = initialize(self.molecule)
        self.nel = get_nel(self.molecule,self.charge)
        self.nclosed,self.nopen = get_open_closed(self.nel,self.multiplicity)
        self.Enuke = get_enuke(self.molecule)
        self.energy = 0
        self.method = "MINDO3"
        self.nbf = get_nbf(self.molecule)
        self.eref = get_reference_energy(self.molecule)
        self.F0 = get_F0(self.molecule)
        self.get_guess()
        return

    def get_guess(self):
        from PyQuante.MINDO3 import get_guess_D
        self.D = get_guess_D(self.molecule)
        return

    def update_fock(self):
        from PyQuante.MINDO3 import get_F1, get_F2
        self.F1 = get_F1(self.molecule,self.D)
        self.F2 = get_F2(self.molecule,self.D)
        self.F = self.F0+self.F1+self.F2
        
    def solve_fock(self):
        from PyQuante.NumWrap import eigh
        from PyQuante.LA2 import mkdens
        self.orbe,self.orbs = eigh(self.F)
        self.D = 2*mkdens(self.orbs,0,self.nclosed)
        
    def calculate_energy(self):
        from PyQuante.LA2 import TraceProperty
        from PyQuante.MINDO3 import ev2kcal
        self.Eel = 0.5*TraceProperty(self.D,self.F0+self.F)
        self.Etot = self.Eel+self.Enuke
        self.energy = self.Etot*ev2kcal+self.eref
        
    def print_iteration_info(self):
        logging.debug("%d %10.4f %10.4f %10.4f" %
                  (self.iter,self.energy,self.Eel,self.Enuke))
        return
        
    def print_pre_iteration_info(self):
        logging.debug("Beginning SCF Optimization")
        logging.debug("It   Etot        Eel        Enuke")
        logging.debug("--   -------     ------     ------")
        return

    def is_converged(self):
        if not self.energy_history:
            self.energy_history.append(self.energy)
            return False
        is_converged_value = abs(self.energy-self.energy_history[-1]) < self.etol
        self.energy_history.append(self.energy)
        return is_converged_value

def test():
    from PyQuante import Molecule
    logging.basicConfig(level=logging.DEBUG,
                        format="%(message)s")
    h2o = Molecule('h2o',[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))])

    h2o_mindo = MINDO3Solver(h2o)
    h2o_mindo.iterate()
    h2o_mindo.met_target(-48.825159)

    from PyQuante.MINDO3 import scf
    scf(h2o,verbose=True)
    

if __name__ == '__main__': test()
