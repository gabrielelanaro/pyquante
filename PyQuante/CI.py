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
from NumWrap import zeros,dot,matrixmultiply,eigh
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

def CIS(Ints,orbs,orbe,nocc,nvirt,ehf):
    CIMatrix = CISMatrix(Ints,orbs,ehf,orbe,nocc,nvirt)
    Ecis,Vectors = eigh(CIMatrix)
    return Ecis

def get_occ_unocc(occs):
    ntot = len(occs)
    nocc = 0
    for fi in occs:
        if fi>0.01: nocc += 1
    return nocc,ntot-nocc

def CISMatrix(Ints,orbs,Ehf,orbe,nocc,nvirt):
    "Naive implementation: Int xfrm + slow formation"
    # The best reference for this stuff is Chap 4 of Szabo/Ostlund

    # Generate the list, and the number of excitations:
    singles = SingleExcitations(range(nocc),range(nocc,nocc+nvirt))
    nex = len(singles)

    # Do the four-index transformation of the 2e ints. This is expensive!
    MOInts = TransformInts(Ints,orbs)

    # Build the CI matrix using the Slater Condon rules
    # see Szabo/Ostlund Table 4.1
    CIMatrix = zeros((nex,nex),'d')
    for ar in xrange(nex):
        a,r = singles[ar]
        for bs in xrange(nex):
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
    for ar in xrange(nsin):
        a,r = singles[ar]
        for bs in xrange(nsin):
            b,s = singles[bs]
            rabs = ijkl2intindex(r,a,b,s)
            rsba = ijkl2intindex(r,s,b,a)
            CIMatrix[ar,bs] = 2*MOInts[rabs] - MOInts[rsba]
            if r==s and a==b: CIMatrix[ar,bs] += Ehf+orbe[r]-orbe[a]
            CIMatrix[bs,ar] = CIMatrix[ar,bs]
    return CIMatrix
    
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
                    temp[mu,nu,sigma,l] = dot(orbs[:,l],tempvec)

    # Transform sigma -> k
    for mu in bfs:
        for nu in bfs:
            for l in mos:
                for k in mos:
                    temp2[mu,nu,k,l] = dot(orbs[:,k],temp[mu,nu,:,l])

    # Transform nu -> j
    for mu in bfs:
        for k in mos:
            for l in mos:
                for j in mos:
                    temp[mu,j,k,l] = dot(orbs[:,j],temp2[mu,:,k,l])

    # Transform mu -> i and repack integrals:
    MOInts = zeros(totlen,'d')
    for i in mos:
        for j in xrange(i+1):
            ij = i*(i+1)/2+j
            for k in mos:
                for l in xrange(k+1):
                    kl = k*(k+1)/2+l
                    if ij >= kl:
                        ijkl = ijkl2intindex(i,j,k,l)
                        MOInts[ijkl] = dot(orbs[:,i],temp[:,j,k,l])

    del temp,temp2,tempvec #force garbage collection now
    return MOInts


def test():
    from Ints import getbasis,getints
    from hartree_fock import rhf
    from Molecule import Molecule

    # Make a test molecule for the calculation
    h2 = Molecule('h2',[(1,(1.,0,0)),(1,(-1.,0,0))])

    # Get a basis set and compute the integrals.
    # normally the routine will do this automatically, but we
    # do it explicitly here so that we can pass the same set
    # of integrals into the CI code and thus not recompute them.
    bfs = getbasis(h2)
    S,h,Ints = getints(bfs,h2)

    # Compute the HF wave function for our molecule
    en,orbe,orbs = rhf(h2,
                       integrals=(S,h,Ints)
                       )
    print "SCF completed, E = ",en
    print " orbital energies ",orbe

    # Compute the occupied and unoccupied orbitals, used in the
    # CIS program to generate the excitations
    nclosed,nopen = h2.get_closedopen()
    nbf = len(bfs)
    nocc = nclosed+nopen
    nvirt = nbf-nocc

    # Call the CI program:
    Ecis = CIS(Ints,orbs,orbe,nocc,nvirt,en)
    print "Ecis = ",Ecis

if __name__ == '__main__':
    test()

