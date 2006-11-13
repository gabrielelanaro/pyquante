#!/usr/bin/env python
"Neon using Gaussians"

from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule

# GAMESS-UK HF Energy
# Energy -128.4744065199

energy = -128.474406 # Changed 2003-04-07 to reflect diis
name = "Ne"

def main():
    ne = Molecule('Ne',atomlist = [(10,(0,0,0))])
    en,orbe,orbs = rhf(ne)
    return en

if __name__ == '__main__': main()    
    
