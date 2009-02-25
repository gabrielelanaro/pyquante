#!/usr/bin/env python
"""\
 Methane correlation energy via MP2. This file is not included in the
 test suite because it takes a long time to perform the integral transform.
"""

from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule
from PyQuante.MP import MP2

energy = -40.362844          # Original number in file
energy = -40.36985496        # (GAMESS_UK)
mp2_correction = -0.16815019 # (GAMESS_UK)
name = "CH4_MP2"

def main():
    atoms = Molecule('ch4',
                     [(6,( .0000000000, .0000000000, .0000000000)),
                      (1,( .0000000000, .0000000000,1.0836058890)),
                      (1,(1.0216334297, .0000000000,-.3612019630)),
                      (1,(-.5108167148, .8847605034,-.3612019630)),
                      (1,(-.5108167148,-.8847605034,-.3612019630))],
                     units='Angstroms')
    bfs = getbasis(atoms)
    nel = atoms.get_nel()
    nbf = len(bfs)
    nocc = nel/2
    S,h,Ints = getints(bfs,atoms)
    en,orbe,orbs = rhf(atoms,integrals=(S,h,Ints))
    print "SCF completed, E = ",en 

    emp2 = MP2(Ints,orbs,orbe,nocc,nbf-nocc)
    print "MP2 correction = ",emp2 
    print "Final energy = ",en+emp2 
    return en+emp2

if __name__ == '__main__': main() 


    

    
