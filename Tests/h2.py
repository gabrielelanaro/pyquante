#!/usr/bin/env python
"H2 using Gaussians"

# Jaguar results for 6-31G** basis and 1 bohr distance:
# Energy   -1.082 099
# Enuke     1.000
# Eone     -2.826 426
# Etwo      0.744 333

from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule

energy = -1.082098
name = "H2"

def main():
    h2 = Molecule('h2',atomlist=[(1,(0,0,0)),(1,(1.0,0,0))])
    en,orbe,orbs = rhf(h2)
    return en

def profmain():
    import profile,pstats
    profile.run('main()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__': main()


    

    
