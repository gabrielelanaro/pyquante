#!/usr/bin/env python
"""\
Reproduce the data from Table 1 in Johnson, Gill and Pople's Paper:
Performance of a Family of Density Functional Methods.
"""

from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

def main(**opts):
    functional = opts.get('functional','BLYP')
    h = Molecule('H',atomlist=[(1,(0,0,0))],
                 multiplicity=2)
    he = Molecule('He',atomlist=[(2,(0,0,0))])
    li = Molecule('Li',atomlist=[(3,(0,0,0))],
                  multiplicity=2)
    be = Molecule('Be',atomlist=[(4,(0,0,0))])
    b = Molecule('B',atomlist=[(5,(0,0,0))],
                 multiplicity=2)
    c = Molecule('C',atomlist=[(6,(0,0,0))],
                 multiplicity=3)
    n = Molecule('N',atomlist=[(7,(0,0,0))],
                 multiplicity=4)
    o = Molecule('O',atomlist=[(8,(0,0,0))],
                 multiplicity=3)
    f = Molecule('F',atomlist=[(9,(0,0,0))],
                 multiplicity=2)
    ne = Molecule('Ne',atomlist=[(10,(0,0,0))])
    #for atoms in [h,he,li,be,b,c,n,o,f,ne]:
    # Just do singlets now, since higher multiplicities are off
    #  b/c of Average Open Shell issues
    print "Running Pople atom tests using %s functional" % functional
    for atoms in [he,be,ne]:
        en,orbe,orbs = dft(atoms,functional=functional)
        print atoms.name, en
    return

if __name__ == '__main__': main()

    

    
