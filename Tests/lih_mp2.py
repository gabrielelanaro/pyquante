#!/usr/bin/env python
from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule
from PyQuante.MP import MP2

# jaguar scf energy: -7.981340
# jaguar mp2 corr    -0.020364
# jaguar mp2 energy: -8.001704

# Gamess-UK energy 
# Energy -8.00219125

energy = -8.001704    # from Jaguar
energy = -8.0021911   # Result obtained 2003-12-18
name = "LiH_MP2"

def main():
    LiH = Molecule('lih',
                     [(3,( .0000000000, .0000000000, .0000000000)),
                      (1,( .0000000000, .0000000000,1.629912))],
                     units='Angstroms')
    bfs = getbasis(LiH)
    nbf = len(bfs)
    nocc,nopen = LiH.get_closedopen()
    assert nopen==0
    S,h,Ints = getints(bfs,LiH)
    en,orbe,orbs = rhf(LiH,integrals=(S,h,Ints))
    print "SCF completed, E = ",en 

    emp2 = MP2(Ints,orbs,orbe,nocc,nbf-nocc)
    print "MP2 correction = ",emp2 
    print "Final energy = ",en+emp2 
    return en+emp2

if __name__ == '__main__': main() 


    

    
