"""\
 MP.py: Moller Plesset Perturbation Theory

 MP2 has been tested for H2. Should work for all closed shell systems.

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from PyQuante.cints import ijkl2intindex
from NumWrap import zeros,dot

VERBOSE=0

def TransformIntsMP2(Ints,orbs,nclosed):
    """\
    O(N^5) 4-index transformation of the two-electron integrals.
    Only transform the ones needed for MP2, which reduces the
    scaling to O(nN^4), where n are the occs (<<N).
    """

    from time import time

    t0 = time()

    nbf,nmo = orbs.shape
    totlen = nmo*(nmo+1)*(nmo*nmo+nmo+2)/8

    occs = range(nclosed)
    mos = range(nmo)
    bfs = range(nbf)

    # Start with (mu,nu|sigma,eta)
    # Unpack aoints and transform sigma -> b
    temp = zeros((nbf,nbf,nclosed,nbf),'d')
    tempvec = zeros(nbf,'d')
    for mu in bfs:
        for nu in bfs:
            for eta in bfs:
                for b in occs:
                    for sigma in bfs:
                        tempvec[sigma] = Ints[ijkl2intindex(mu,nu,sigma,eta)]
                    temp[mu,nu,b,eta] = dot(orbs[:,b],tempvec)

    temp2 = zeros((nclosed,nbf,nclosed,nbf),'d')
    for nu in bfs:
        for eta in bfs:
            for b in occs:
                for a in occs:
                    temp2[a,nu,b,eta] = dot(orbs[:,a],temp[:,nu,b,eta])

    temp = zeros((nclosed,nbf,nclosed,nmo),'d')
    for a in occs:
        for nu in bfs:
            for b in occs:
                for j in mos:
                    temp[a,nu,b,j] = dot(orbs[:,j],temp2[a,nu,b,:])

    # Transform mu -> i and repack integrals:
    MOInts = zeros(totlen,'d')
    for a in occs:
        for j in mos:
            for b in occs:
                for i in mos:
                    aibj = ijkl2intindex(a,i,b,j)
                    MOInts[aibj] = dot(orbs[:,i],temp[a,:,b,j])

    #print "Integral transform time = ",time()-t0
    del temp,temp2,tempvec #force garbage collection now
    return MOInts

def MP2(aoints,orbs,orbe,nclosed,nvirt):
    #moints = TransformInts(aoints,orbs)
    moints = TransformIntsMP2(aoints,orbs,nclosed)
    occs = range(nclosed)
    unoccs = range(nclosed,nclosed+nvirt)
    nocc = len(occs)

    Epairs = zeros((nocc,nocc),'d')

    Emp2 = 0
    for a in occs:
        for b in occs:
            for r in unoccs:
                for s in unoccs:
                    arbs = moints[ijkl2intindex(a,r,b,s)]
                    asbr = moints[ijkl2intindex(a,s,b,r)]
                    Epairs[a,b] += arbs*(2*arbs-asbr)/\
                                   (orbe[a]+orbe[b]-orbe[r]-orbe[s])
    if VERBOSE:
        print "MP2 pair energies"
        for a in xrange(nocc):
            for b in xrange(a):
                print a,b,Epairs[a,b]+Epairs[b,a]
            print a,a,Epairs[a,a]
    return sum(sum(Epairs))

def EN2(aoints,orbs,orbe,nclosed,nvirt):
    moints = TransformIntsMP2(aoints,orbs,nclosed)
    occs = range(nclosed)
    unoccs = range(nclosed,nclosed+nvirt)
    nocc = len(occs)

    Epairs = zeros((nocc,nocc),'d')

    Emp2 = 0
    for a in occs:
        for b in occs:
            for r in unoccs:
                for s in unoccs:
                    arbs = moints[ijkl2intindex(a,r,b,s)]
                    asbr = moints[ijkl2intindex(a,s,b,r)]
                    Epairs[a,b] += arbs*(2*arbs-asbr)/\
                                   (orbe[a]+orbe[b]-orbe[r]-orbe[s])
    if VERBOSE:
        print "EN2 pair energies"
        for a in xrange(nocc):
            for b in xrange(a):
                print a,b,Epairs[a,b]+Epairs[b,a]
            print a,a,Epairs[a,a]
    return sum(sum(Epairs))

def test():
    # Jaguar gets -0.0276516 h for this:
    from Ints import getbasis,getints
    from hartree_fock import scf
    from IO import mtx2file
    from Molecule import Molecule

    atoms = Molecule('h2',[(1,(1.,0,0)),(1,(-1.,0,0))])
    bfs = getbasis(atoms)
    S,h,Ints = getints(bfs,atoms)
    en,orbe,orbs = scf(atoms,S,h,Ints,0,0.0001,10)
    print "SCF completed, E = ",en

    emp2 = MP2(Ints,orbs,orbe,1,9)
    print "MP2 correction = ",emp2
    print "Final energy = ",en+emp2
    return

if __name__ == '__main__': test()

                                                         
