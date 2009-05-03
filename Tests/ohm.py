#!/usr/bin/env python
"H2O using Gaussians"

# Jaguar results for 6-31G** basis and 0.957939 angstrom distance:
# Energy -75.332657

import time
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule
from PyQuante.Constants import ang2bohr

energy = -75.332606
name = "OH-"

def main():
    oh=Molecule('OH-',
                atomlist = [(8,(0,0,0)),(1,(0.957939*ang2bohr,0,0))],
                charge=-1)
    en,orbe,orbs = rhf(oh)
    return en


def profmain():
    import cProfile,pstats
    cProfile.run('main()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__':
    t0 = time.time()
    en = main()
    #profmain()
    t1 = time.time()
    print en,energy,t1-t0


    
