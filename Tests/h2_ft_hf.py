#!/usr/bin/env python
from PyQuante.NumWrap import arange
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule

energy = -1.130284
name = 'H2/FT/HF'

def main():
    h2 = Molecule('h2',
                  atomlist = [(1,(0.,0.,0.7)),(1,(0.,0.,-0.7))],
                  units = 'Bohr')
    en,orbe,orbs = rhf(h2,ETemp=1e4)
    print "Energy ",en,abs(en-energy)
    print "Spectrum ",orbe
    return en

if __name__ == '__main__': main()
    
                                   
    
