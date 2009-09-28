#!/usr/bin/env python
"""\
 Test the 2e integral code.

"""

from PyQuante.pyints import contr_coulomb as pycc
from PyQuante.rys import contr_coulomb as ryscc
from PyQuante.hgp import contr_coulomb as hgpcc
from PyQuante.cints import contr_coulomb as ccc
from PyQuante.crys import contr_coulomb as cryscc
from PyQuante.chgp import contr_coulomb as chgpcc
from PyQuante.CGBF import CGBF
from PyQuante.cints import ijkl2intindex
from PyQuante.Ints import getbasis
from PyQuante.Molecule import Molecule

from time import time

def get2ints(bfs,coul_func):
    """Store integrals in a long array in the form (ij|kl) (chemists
    notation. We only need i>=j, k>=l, and ij <= kl"""
    from array import array
    nbf = len(bfs)
    totlen = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
    Ints = array('d',[0]*totlen)
    for i in range(nbf):
        for j in range(i+1):
            ij = i*(i+1)/2+j
            for k in range(nbf):
                for l in range(k+1):
                    kl = k*(k+1)/2+l
                    if ij >= kl:
                        ijkl = ijkl2intindex(i,j,k,l)
                        Ints[ijkl] = coulomb(bfs[i],bfs[j],bfs[k],bfs[l],
                                             coul_func)
    return Ints

def coulomb(a,b,c,d,coul_func):
    "Coulomb interaction between 4 contracted Gaussians"
    Jij = coul_func(a.exps(),a.coefs(),a.pnorms(),a.origin(),a.powers(),
                    b.exps(),b.coefs(),b.pnorms(),b.origin(),b.powers(),
                    c.exps(),c.coefs(),c.pnorms(),c.origin(),c.powers(),
                    d.exps(),d.coefs(),d.pnorms(),d.origin(),d.powers())
    return a.norm()*b.norm()*c.norm()*d.norm()*Jij

def maxdiff(a,b):
    md = -1e10
    for i in range(len(a)):
        md = max(md,a[i]-b[i])
    return md

def test():
    #from PyQuante.Basis.sto3g import basis_data
    from PyQuante.Basis.p631ss import basis_data
    r = 1/0.52918
    atoms=Molecule('h2o',atomlist = [(8,(0,0,0)),(1,(r,0,0)),(1,(0,r,0))])

    inttol = 1e-6 # Tolerance to which integrals must be equal
    
    bfs = getbasis(atoms,basis_data)
    print "Int times: "
    t0 = time()

    int0 = get2ints(bfs,chgpcc)
    t1 = time()
    print "CHGP Ints:    ",t1-t0

    int1 = get2ints(bfs,cryscc)
    t2 = time()
    print "CRys Ints:   ",t2-t1
    assert maxdiff(int0,int1)<inttol
    

    int1 = get2ints(bfs,ccc)
    t3 = time()
    print "CINTS Ints:   ",t3-t2
    assert maxdiff(int0,int1)<inttol

    ints1 = get2ints(bfs,hgpcc)
    t4 = time()
    print "HGP Ints:     ",t4-t3
    assert maxdiff(int0,int1)<inttol

    int1 = get2ints(bfs,ryscc)
    t5 = time()
    print "Rys Ints:  ",t5-t4
    assert maxdiff(int0,int1)<inttol

    int1 = get2ints(bfs,pycc)
    t6 = time()
    print "Py Ints:  ",t6-t5
    assert maxdiff(int0,int1)<inttol


if __name__ == '__main__': test()

    
# Sample times (9/28/09, Macbook Air:)
# Int times: 
# CHGP Ints:     3.02386283875
# CRys Ints:    2.28243303299
# CINTS Ints:    6.17023396492
# HGP Ints:      250.576164007
# Rys Ints:   204.740512133
# Py Ints:   295.842331886
