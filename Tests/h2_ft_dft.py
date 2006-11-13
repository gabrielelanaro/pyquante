#!/usr/bin/env python

from PyQuante.NumWrap import arange
from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

energy = -1.132473 # For etemp=1e4
name = 'H2/FT/DFT'

def main():
    r = 0.7
    h2 = Molecule('h2',
                  atomlist = [(1,(0,0,r/2.)),
                              (1,(0,0,-r/2.))],
                  units='Angs')
    en,orbe,orbs = dft(h2,ETemp=1e4)
    return en

if __name__ == '__main__': main()
    
                                   
    
