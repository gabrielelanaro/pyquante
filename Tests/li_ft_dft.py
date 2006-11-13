#!/usr/bin/env python

from PyQuante.NumWrap import arange
from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

energy = -7.349422
name = "Li/FT/DFT"

def main():
    li = Molecule('Li',
                  atomlist = [(3,(0,0,0))],
                  units='Angs',
                  multiplicity=2)
    en,orbe,orbs = dft(li,ETemp=1e4)
    return en

if __name__ == '__main__': main()
