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
from PyQuante.Ints import getbasis,ijkl2intindex
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

def test_timings():
    from PyQuante.basis_sto3g import basis_data
    r = 1/0.52918
    atoms=Molecule('h2o',atomlist = [(8,(0,0,0)),(1,(r,0,0)),(1,(0,r,0))])
    bfs = getbasis(atoms,basis_data)
    #print "Int times: "
    #t0 = time()
    #get2ints(bfs,chgpcc)
    #t1 = time()
    #print "CHGP Ints:    ",t1-t0
    #get2ints(bfs,cryscc)
    #t2 = time()
    #print "CRys Ints:   ",t2-t1
    #get2ints(bfs,ccc)
    t3 = time()
    #print "CINTS Ints:   ",t3-t2
    get2ints(bfs,hgpcc)
    t4 = time()
    print "HGP Ints:     ",t4-t3
    #get2ints(bfs,ryscc)
    #t5 = time()
    #print "Rys Ints:  ",t5-t4
    #get2ints(bfs,pycc)
    #t6 = time()
    #print "Py Ints:  ",t6-t5

def test_vrr():
    from PyQuante.hgp import vrr,contr_vrr
    from chgp import vrr as cvrr

    xyza = (0,0,0)
    xyzb = (0,0,0)
    xyzc = (0,0,0)
    xyzd = (0,0,0)

    lmna = (2,0,0)
    lmnc = (0,0,2)

    alphaa = 1.
    alphab = 2.
    alphac = 3.
    alphad = 1.

    norma = normb = normc = normd = 1.

    M = 0

    vrr1 = vrr(xyza,norma,lmna,alphaa,
               xyzb,normb,alphab,
               xyzc,normc,lmnc,alphac,
               xyzd,normd,alphad,M)
    vrr2 = cvrr(xyza,norma,lmna,alphaa,
               xyzb,normb,alphab,
               xyzc,normc,lmnc,alphac,
               xyzd,normd,alphad,M)

    print vrr1, vrr2

def test():
    r = 1/0.52918
    atoms=Molecule('name',atomlist = [(8,(0,0,0))])
    bfs = getbasis(atoms)
    j1 = coulomb(bfs[2],bfs[2],bfs[11],bfs[11],chgpcc)
    j2 = coulomb(bfs[11],bfs[11],bfs[2],bfs[2],chgpcc)
    print j1,j2
    

def test_all():
    r = 1/0.52918
    atoms=Molecule('name',atomlist = [(8,(0,0,0))])
    bfs = getbasis(atoms)

    nbf = len(bfs)
    for i in range(nbf):
        for j in range(i+1):
            for k in range(nbf):
                for l in range(k+1):
                    ij = i*(i+1)/2+j
                    kl = k*(k+1)/2+l
                    if ij < kl:
                        j1 = coulomb(bfs[i],bfs[j],bfs[k],bfs[l],chgpcc)
                        j2 = coulomb(bfs[k],bfs[l],bfs[i],bfs[j],chgpcc)
                        if abs(j1-j2)>1e-6:
                            print i,j,k,l
                            print bfs[i].powers(),bfs[j].powers(),\
                                  bfs[k].powers(),bfs[k].powers()
    

if __name__ == '__main__': test()

    
