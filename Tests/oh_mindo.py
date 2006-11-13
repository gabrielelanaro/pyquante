#!/usr/bin/env python
"H2O using Gaussians"

import time
from PyQuante.MINDO3 import scf
from PyQuante.Molecule import Molecule

energy = 18.127533
name = "OH MINDO/3"

def main():
    atomlist = Molecule('oh',atomlist = [(8,(0,0,0)),(1,(1.,0,0))])
    en = scf(atomlist)
    return en


def profmain():
    import profile,pstats
    profile.run('main()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__': main()


    
