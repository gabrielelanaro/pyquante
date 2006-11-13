from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

name = "H atom/DFT"
energy = -0.4415033 # open shell average occupation

def main():
    atomlist = Molecule('H',atomlist = [(1,(0,0,0))],multiplicity=2)
    en,orbe,orbs = dft(atomlist)
    return en

if __name__ == '__main__': main()
