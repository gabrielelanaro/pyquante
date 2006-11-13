#!/usr/bin/env python
"""\
 CI.py: Configuration Interaction Routines.

 Currently, the only one implemented is CI-S.

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

import os,sys
from PyQuante.cints import ijkl2intindex
from NumWrap import zeros,dot,matrixmultiply
from NumWrap import Heigenvectors
from Ints import getbasis, get2ints

def SingleExcitations(occs,virts):
    singles = []
    for occ in occs:
        for virt in virts:
            singles.append((occ,virt))
    return singles

def DoubleExcitations(occs,virts):
    doubles = []
    for occa in occs:
        for occb in occs:
            for virta in virts:
                for virtb in virts:
                    doubles.append((occa,occb,virta,virtb))
    return doubles

def CIS(Ints,orbs,orbe,occs,ehf):
    CIMatrix = CISMatrix(Ints,orbs,ehf,orbe,occs)
    Ecis,Vectors = Heigenvectors(CIMatrix)
    return Ecis

def get_occ_unocc(occs):
    ntot = len(occs)
    nocc = 0
    for fi in occs:
        if fi>0.01: nocc += 1
    return nocc,ntot-nocc

def CISMatrix(Ints,orbs,Ehf,orbe,occs):
    "Naive implementation: Int xfrm + slow formation"
    # The best reference for this stuff is Chap 4 of Szabo/Ostlund

    nocc, nvirt = get_occ_unocc(occs)
    singles = SingleExcitations(range(nocc),range(nocc,nocc+nvirt))
    nex = len(singles)

    MOInts = TransformInts(Ints,orbs)

    # see Szabo/Ostlund Table 4.1
    CIMatrix = zeros((nex,nex),'d')
    for ar in range(nex):
        a,r = singles[ar]
        for bs in range(nex):
            b,s = singles[bs]
            rabs = ijkl2intindex(r,a,b,s)
            rsba = ijkl2intindex(r,s,b,a)
            CIMatrix[ar,bs] = 2*MOInts[rabs] - MOInts[rsba]
            if r==s and a==b: CIMatrix[ar,bs] += Ehf+orbe[r]-orbe[a]
            CIMatrix[bs,ar] = CIMatrix[ar,bs]
    return CIMatrix

def CISDMatrix(Ints,orbs,Ehf,orbe,occs):
    nocc, nvirt = get_occ_unocc(occs)
    singles = SingleExcitations(range(nocc),range(nocc,nocc+nvirt))
    #doubles = DoubleExcitations(range(nocc),range(nocc,nocc+nvirt))
    doubles = []
    nsin = len(singles)
    ndoub = len(doubles)
    nex = nsin+ndoub
    
    MOInts = TransformInts(Ints,orbs)

    # see Szabo/Ostlund Table 4.1
    CIMatrix = zeros((nex,nex),'d')
    for ar in range(nsin):
        a,r = singles[ar]
        for bs in range(nsin):
            b,s = singles[bs]
            rabs = ijkl2intindex(r,a,b,s)
            rsba = ijkl2intindex(r,s,b,a)
            CIMatrix[ar,bs] = 2*MOInts[rabs] - MOInts[rsba]
            if r==s and a==b: CIMatrix[ar,bs] += Ehf+orbe[r]-orbe[a]
            CIMatrix[bs,ar] = CIMatrix[ar,bs]
    return CIMatrix
    
def TransformIntsSlo(Ints,orbs):
    """This is really, really slow (N^8). TransformInts should be used
    instead. This routine only exists for purposes of comparison."""

    nbf,nmo = orbs.shape # probably wrong if not square
    totlen = nmo*(nmo+1)*(nmo*nmo+nmo+2)/8
    MOInts = zeros(totlen,'d')
    
    for i in range(nmo):
        for j in range(i+1):
            ij = i*(i+1)/2+j
            for k in range(nmo):
                for l in range(k+1):
                    kl = k*(k+1)/2+l
                    if ij >= kl:
                        ijkl = ijkl2intindex(i,j,k,l)
                        val = 0
                        for mu in range(nbf):
                            for nu in range(nbf):
                                for sigma in range(nbf):
                                    for eta in range(nbf):
                                        mnse = ijkl2intindex(mu,nu,sigma,eta)
                                        val += Ints[mnse]*\
                                               orbs[i,mu]*orbs[j,nu]*\
                                               orbs[k,sigma]*orbs[l,eta]
                        MOInts[ijkl] = val
    return MOInts

def TransformInts(Ints,orbs):
    """O(N^5) 4-index transformation of the two-electron integrals. Not as
    efficient as it could be, since it inflates to the full rectangular
    matrices rather than keeping them compressed. But at least it gets the
    correct result."""

    from time import time

    t0 = time()

    nbf,nmo = orbs.shape
    totlen = nmo*(nmo+1)*(nmo*nmo+nmo+2)/8

    temp = zeros((nbf,nbf,nbf,nmo),'d')
    tempvec = zeros(nbf,'d')
    temp2 = zeros((nbf,nbf,nmo,nmo),'d')

    mos = range(nmo) # preform so we don't form inside loops
    bfs = range(nbf)

    # Start with (mu,nu|sigma,eta)
    # Unpack aoints and transform eta -> l
    for mu in bfs:
        for nu in bfs:
            for sigma in bfs:
                for l in mos:
                    for eta in bfs:
                        tempvec[eta] = Ints[ijkl2intindex(mu,nu,sigma,eta)]
                    temp[mu,nu,sigma,l] = dot(orbs[l,:],tempvec)

    # Transform sigma -> k
    for mu in bfs:
        for nu in bfs:
            for l in mos:
                for k in mos:
                    temp2[mu,nu,k,l] = dot(orbs[k,:],temp[mu,nu,:,l])

    # Transform nu -> j
    for mu in bfs:
        for k in mos:
            for l in mos:
                for j in mos:
                    temp[mu,j,k,l] = dot(orbs[j,:],temp2[mu,:,k,l])

    # Transform mu -> i and repack integrals:
    MOInts = zeros(totlen,'d')
    for i in mos:
        for j in range(i+1):
            ij = i*(i+1)/2+j
            for k in mos:
                for l in range(k+1):
                    kl = k*(k+1)/2+l
                    if ij >= kl:
                        ijkl = ijkl2intindex(i,j,k,l)
                        MOInts[ijkl] = dot(orbs[i,:],temp[:,j,k,l])

    del temp,temp2,tempvec #force garbage collection now
    return MOInts


def test():
    from Ints import getbasis,getints
    from hartree_fock import scf
    from IO import mtx2file
    from Molecule import Molecule

    atoms = Molecule('h2',[(1,(1.,0,0)),(1,(-1.,0,0))])
    bfs = getbasis(atoms)
    S,h,Ints = getints(bfs,atoms)
    en,orbe,orbs = scf(atoms,S,h,Ints,0,0.0001,10)
    print "SCF completed, E = ",en
    print " orbital energies ",orbe
    occs = [1.]+[0.]*9

    Ecis = CIS(Ints,orbs,orbe,occs,en)
    print Ecis

if __name__ == '__main__':
    test()

