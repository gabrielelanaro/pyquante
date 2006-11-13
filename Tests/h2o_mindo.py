#!/usr/bin/env python
"H2O using Gaussians"

import time
from PyQuante.MINDO3 import scf
from PyQuante.Molecule import Molecule

energy = -48.825159
name = "H2O MINDO/3"

def main():
    atomlist = Molecule('h2o',
                        atomlist = [(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))])
    en = scf(atomlist)
    return en


def profmain():
    import profile,pstats
    profile.run('main()','prof')
    prof = pstats.Stats('prof')
    prof.strip_dirs().sort_stats('time').print_stats(15)

if __name__ == '__main__': main()


    
