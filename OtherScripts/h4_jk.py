#!/usr/bin/env python

from PyQuante.basis_631ss import basis
from PyQuante.hartree_fock import *

# Print out matlab files to display how the J and K matrices of 2 H2
#  molecules fall off with distance.

def getK0(Ints,D,bfs,centermax=2):
    "Form the few-center K approximation"
    nbf = D.shape[0]
    K = zeros((nbf,nbf),Float)
    for i in range(nbf):
        icenter = bfs[i].origin()
        for j in range(i+1):
            jcenter = bfs[j].origin()
            for k in range(nbf):
                kcenter = bfs[k].origin()
                for l in range(nbf):
                    lcenter = bfs[l].origin()
                    centers = [icenter]
                    for xcenter in [jcenter,kcenter,lcenter]:
                        if xcenter not in centers:
                            centers.append(xcenter)
                    ncenter = len(centers)
                    if ncenter <= centermax:
                        index_k1 = ijkl2intindex(i,k,j,l)
                        index_k2 = ijkl2intindex(i,l,k,j)
                        K[i,j] = K[i,j] + \
                                 D[k,l]*0.5*(Ints[index_k1]+Ints[index_k2])
                        K[j,i] = K[i,j]
    return K

def h4_jk():
    for r in range(2,20,2):
        atomlist = [(1,(-r/2.0,0,0.7)),(1,(-r/2.0,0,-0.7)),
                    (1,(r/2.0,0,0.7)),(1,(r/2.0,0,-0.7))]
        bfs = getbasis(atomlist,basis)
        S,h = get1ints(bfs,atomlist)
        Ints = get2ints(bfs)
        en,orbe,orbs = scf(atomlist,S,h,Ints)
        enuke = get_enuke_from_atomlist(atomlist)
        D = mkdens(orbs,0,2)
        J = getJ(Ints,D)
        K = getK(Ints,D)
        K0 = getK0(Ints,D,bfs,2)
        en0 = get_energy(h,h+2*J-K0,D,enuke)
        print "H2-H2 separation= %10.4f, Ens = %10.4f, %10.4f" %\
              (r*0.52918,en,en0)
        write_mat('j%2d.dat' % r,J)
        write_mat('k%2d.dat' % r,K)
        write_mat('k0%2d.dat' % r,K0)
    return

def write_mat(filename,A):
    file = open(filename,'w')
    n,m = A.shape
    file.write('%d %d\n' % (n,m))
    for i in range(n):
        for j in range(m):
            file.write('%16.10f ' % A[i,j])
        file.write('\n')
    file.close()
    return

if __name__ == '__main__': h4_jk()
