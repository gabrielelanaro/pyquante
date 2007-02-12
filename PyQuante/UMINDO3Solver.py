from PyQuante import logging
from MINDO3Solver import MINDO3Solver

class UMINDO3Solver(MINDO3Solver):
    def setup(self):
        MINDO3Solver.setup(self)
        self.nalpha = self.nclosed+self.nopen
        self.nbeta = self.nclosed
        return
    
    def get_guess(self):
        from PyQuante.MINDO3 import get_guess_D
        self.Da = self.Db = 0.5*get_guess_D(self.molecule)
        return

    def update_fock(self):
        from PyQuante.MINDO3 import get_F1_open, get_F2_open
        F1a = get_F1_open(self.molecule,self.Da,self.Db)
        F1b = get_F1_open(self.molecule,self.Db,self.Da)
        F2a = get_F2_open(self.molecule,self.Da,self.Db)
        F2b = get_F2_open(self.molecule,self.Db,self.Da)
        self.Fa = self.F0+F1a+F2a
        self.Fb = self.F0+F1b+F2b
        return

    def solve_fock(self):
        from PyQuante.NumWrap import eigh
        from PyQuante.LA2 import mkdens
        self.orbea,self.orbsa = eigh(self.Fa)
        self.orbeb,self.orbsb = eigh(self.Fb)
        self.Da = mkdens(self.orbsa,0,self.nalpha)
        self.Db = mkdens(self.orbsb,0,self.nbeta)

    def calculate_energy(self):
        from PyQuante.LA2 import TraceProperty
        from PyQuante.MINDO3 import ev2kcal
        self.Eel = 0.5*TraceProperty(self.Da,self.F0+self.Fa)+\
                   0.5*TraceProperty(self.Db,self.F0+self.Fb)
        self.Etot = self.Eel+self.Enuke
        self.energy = self.Etot*ev2kcal+self.eref

def test():
    from PyQuante import Molecule
    logging.basicConfig(level=logging.DEBUG,
                        format="%(message)s")
    oh = Molecule('oh',[(8,(0,0,0)),(1,(1.,0,0))],multiplicity=2)

    oh_mindo = UMINDO3Solver(oh)
    oh_mindo.iterate()
    oh_mindo.met_target(18.127533)

    from PyQuante.MINDO3 import scf
    scf(oh,verbose=True)

if __name__ == '__main__': test()
