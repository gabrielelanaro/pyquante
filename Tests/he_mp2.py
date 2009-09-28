#!/usr/bin/env python
from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule
from PyQuante.MP import MP2

# jag results:
# scf: -2.855 160
# dmp2 -0.025 477
# emp2 -2.880 637

energy = -2.880637
name = "He_MP2"

def main():
    atoms = Molecule('He',
                     [(2,( .0000000000, .0000000000, .0000000000))],
                     units='Angstroms')
    bfs = getbasis(atoms)
    nel = atoms.get_nel()
    nbf = len(bfs)
    nocc = nel/2
    S,h,Ints = getints(bfs,atoms)
    en,orbe,orbs = rhf(atoms,integrals=(S,h,Ints))

    emp2 = MP2(Ints,orbs,orbe,nocc,nbf-nocc)
    return en+emp2

def isnear(a,b,tol=1e-4): return abs(a-b)<tol

if __name__ == '__main__':
    print isnear(main(),energy)


    

    
