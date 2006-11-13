#!/usr/bin/env python

from PyQuante.basis_631ss import basis
from PyQuante.hartree_fock import *

# Print out matlab files to display how the J and K matrices of H2 fall
#  off with distance.

def h2_jk():
    for r in range(1,10):
        atomlist = [(1,(0,0,r/2.)),(1,(0,0,-r/2.))]
        bfs = getbasis(atomlist,basis)
        S,h = get1ints(bfs,atomlist)
        Ints = get2ints(bfs)
        en,orbe,orbs = scf(atomlist,S,h,Ints)
        D = mkdens(orbs,0,1)
        J = getJ(Ints,D)
        K = getK(Ints,D)
        write_matlab('j_%d.dat' % r,J)
        write_matlab('k_%d.dat' % r,K)
    return

def write_matlab(filename,A):
    file = open(filename,'w')
    n,m = A.shape
    file.write('%d %d\n' % (n,m))
    for i in range(n):
        for j in range(m):
            file.write('%16.10f ' % A[i,j])
        file.write('\n')
    file.close()
    return

if __name__ == '__main__': h2_jk()
