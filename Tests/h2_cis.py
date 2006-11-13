#!/usr/bin/env python
"H2 first excited state energy"

from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule
from PyQuante.CI import CIS

energy = -0.645236
name = "H2_CIS"

def main():
    atoms = Molecule('h2',[(1,(1.,0,0)),(1,(-1.,0,0))])
    bfs = getbasis(atoms)
    S,h,Ints = getints(bfs,atoms)
    en,orbe,orbs = rhf(atoms,integrals=(S,h,Ints))
    occs = [1.]+[0.]*9

    Ecis = CIS(Ints,orbs,orbe,occs,en)
    return Ecis[0]

if __name__ == '__main__': main()


    

    
