from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

name = "Li atom/DFT"
energy = -7.3321  # open shell average occupation

def main():
    atomlist = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)
    en,orbe,orbs = dft(atomlist,verbose=True)
    return en

if __name__ == '__main__': main()
