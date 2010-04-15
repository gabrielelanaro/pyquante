"""\
 Tom Manz's Simplified Size-Consistent Configuration Interaction
 This may be equivalent to second-order Epstein-Nesbet pair
 correlation theory.

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Note: for this routine to work properly, we need to use the
#  cints (or perhaps the chgp) integrals from the various basis
#  function routines. This is not particularly elegantly done
#  right now -- one has to comment out the appropriate routines
#  in PGBF or CGBF.

from PyQuante.cints import ijkl2intindex
from PyQuante.LA2 import mkdens
from PyQuante.Ints import getbasis,getints, getJ, getK, getT
from PyQuante.hartree_fock import get_energy,uhf
from PyQuante.IO import mtx2file
from NumWrap import zeros,dot,sqrt,choose

from time import time

def TransformInts(Ints,orbs1,orbs2, nocc):

    nbf,nmo = orbs1.shape
    totlen = nmo*nmo*nmo*nmo

    occs = range(nocc)
    mos = range(nmo)
    bfs = range(nbf)

    # Start with (mu,nu|sigma,eta)
    # Unpack aoints and transform sigma -> b
    # Here sigma, b, nu, j are of first same spin group,
    #  others of second same spin group
    temp = zeros((nbf,nbf,nbf,nbf),'d')
    tempvec = zeros(nbf,'d')
    for mu in bfs:
        for nu in bfs:
            for eta in bfs:
                for b in bfs:
                    for sigma in bfs:
                        tempvec[sigma] = Ints[ijkl2intindex(mu,nu,sigma,eta)]
                    temp[mu,nu,b,eta] = dot(orbs1[b,:],tempvec)

    temp2 = zeros((nbf,nbf,nbf,nbf),'d')
    for nu in bfs:
        for eta in bfs:
            for b in bfs:
                for a in bfs:
                    temp2[a,nu,b,eta] = dot(orbs2[a,:],temp[:,nu,b,eta])

    temp = zeros((nbf,nbf,nbf,nbf),'d')
    for a in bfs:
        for nu in bfs:
            for b in bfs:
                for j in bfs:
                    temp[a,nu,b,j] = dot(orbs1[j,:],temp2[a,nu,b,:])

    # Transform mu -> i and repack integrals:
    MOInts = zeros(totlen,'d')
    for a in bfs:
        for j in bfs:
            for b in bfs:
                for i in bfs:
                    aibj = ijkl2intindex(a,i,b,j)
                    MOInts[aibj] = dot(orbs2[i,:],temp[a,:,b,j])

    del temp,temp2,tempvec #force garbage collection now
    return MOInts, nbf

def EN2(molecule,**opts):#
    "General wrapper for the simple CI method"
    nalpha,nbeta = molecule.get_alphabeta()
    bfs = getbasis(molecule)
    S,h,Ints = getints(bfs,molecule)
    energy,(orbea,orbeb),(orbsa,orbsb) = uhf(molecule,integrals=(S,h,Ints),
                                             bfs=bfs,**opts)
    EHF = energy
    print "The Hatree-Fock energy is ",EHF
    #compute the transformed molecular orbital integrals

    aamoints, nbf = TransformInts(Ints,orbsa,orbsa, nalpha)
    bbmoints, nbf = TransformInts(Ints,orbsb,orbsb, nbeta)
    abmoints, nbf = TransformInts(Ints,orbsa,orbsb, nalpha)
    
    #Initialize the fractional occupations:
    Yalpha = zeros((nbf),'d')
    Ybeta = zeros((nbf),'d')

    #set up the occupied and virtual orbitals
    aoccs = range(nalpha)
    boccs = range(nbeta)
    avirt = range(nalpha,nbf) #numbers of alpha virtual orbitals
    bvirt = range(nbeta,nbf) #numbers of beta virtual orbitals

    ########  Computation of the primary energy correction  ######### 
    #Set initial correction terms to zero
    Ec1 = 0.
    sum = 0.
    z = 1.
   
    #compute correction term for two alpha electrons

    for a in aoccs:
        for b in xrange(a):
            for r in avirt:
                for s in xrange(nalpha,r):
                    arbs = aamoints[ijkl2intindex(a,r,b,s)]
                    asbr = aamoints[ijkl2intindex(a,s,b,r)]
                    rraa = aamoints[ijkl2intindex(r,r,a,a)] - \
                           aamoints[ijkl2intindex(r,a,a,r)]
                    rrbb = aamoints[ijkl2intindex(r,r,b,b)] - \
                           aamoints[ijkl2intindex(r,b,b,r)]
                    ssaa = aamoints[ijkl2intindex(s,s,a,a)] - \
                           aamoints[ijkl2intindex(s,a,a,s)]
                    ssbb = aamoints[ijkl2intindex(s,s,b,b)] - \
                           aamoints[ijkl2intindex(s,b,b,s)]
                    rrss = aamoints[ijkl2intindex(r,r,s,s)] - \
                           aamoints[ijkl2intindex(r,s,s,r)]
                    aabb = aamoints[ijkl2intindex(a,a,b,b)] - \
                           aamoints[ijkl2intindex(a,b,b,a)]

                    eigendif = (orbea[r] + orbea[s] - orbea[a] - orbea[b])
                    delcorr = (-rraa - rrbb - ssaa - ssbb + rrss + aabb) 
                    delta = eigendif + delcorr*z

                    Eio = (arbs - asbr)

                    x = -Eio/delta
                    if abs(x) > 1:
                        print "Warning a large x value has been ",\
                              "discovered with x = ",x
                    x = choose(x < 1, (1,x))
                    x = choose(x > -1, (-1,x))                   
                    sum += x*x
                    Yalpha[a] -= x*x
                    Yalpha[b] -= x*x
                    Yalpha[r] += x*x
                    Yalpha[s] += x*x
                    Ec1 += x*Eio             


    #compute correction term for two beta electrons

    for a in boccs:
        for b in xrange(a):
            for r in bvirt:
                for s in xrange(nbeta,r):
                    arbs = bbmoints[ijkl2intindex(a,r,b,s)]
                    asbr = bbmoints[ijkl2intindex(a,s,b,r)]
                    rraa = bbmoints[ijkl2intindex(r,r,a,a)] - \
                           bbmoints[ijkl2intindex(r,a,a,r)]
                    rrbb = bbmoints[ijkl2intindex(r,r,b,b)] - \
                           bbmoints[ijkl2intindex(r,b,b,r)]
                    ssaa = bbmoints[ijkl2intindex(s,s,a,a)] - \
                           bbmoints[ijkl2intindex(s,a,a,s)]
                    ssbb = bbmoints[ijkl2intindex(s,s,b,b)] - \
                           bbmoints[ijkl2intindex(s,b,b,s)]
                    rrss = bbmoints[ijkl2intindex(r,r,s,s)] - \
                           bbmoints[ijkl2intindex(r,s,s,r)]
                    aabb = bbmoints[ijkl2intindex(a,a,b,b)] - \
                           bbmoints[ijkl2intindex(a,b,b,a)]


                    eigendif = (orbeb[r] + orbeb[s] - orbeb[a] - orbeb[b])
                    delcorr = (-rraa - rrbb - ssaa - ssbb + rrss + aabb) 
                    delta = eigendif + delcorr*z

                    Eio = (arbs - asbr)

                    x = -Eio/delta
                    if abs(x) > 1: print "Warning a large x value has ",\
                       "been discovered with x = ",x
                    x = choose(x < 1, (1,x))
                    x = choose(x > -1, (-1,x))                   
                    sum += x*x
                    Ybeta[a] -= x*x
                    Ybeta[b] -= x*x
                    Ybeta[r] += x*x
                    Ybeta[s] += x*x
                    Ec1 += x*Eio

    #compute correction term for one alpha and one beta electron

    for a in aoccs:
        for b in boccs:
            for r in avirt:
                for s in bvirt:
                    arbs = abmoints[ijkl2intindex(a,r,b,s)]
                    rraa = aamoints[ijkl2intindex(r,r,a,a)] - \
                           aamoints[ijkl2intindex(r,a,a,r)]
                    rrbb = abmoints[ijkl2intindex(r,r,b,b)]
                    aass = abmoints[ijkl2intindex(a,a,s,s)]
                    ssbb = bbmoints[ijkl2intindex(s,s,b,b)] - \
                           bbmoints[ijkl2intindex(s,b,b,s)]
                    rrss = abmoints[ijkl2intindex(r,r,s,s)]
                    aabb = abmoints[ijkl2intindex(a,a,b,b)]

                    eigendif = (orbea[r] + orbeb[s] - orbea[a] - orbeb[b])
                    delcorr = (-rraa - rrbb - aass - ssbb + rrss + aabb)
                    delta = eigendif + delcorr*z

                    Eio = arbs

                    x = -Eio/delta
                    if abs(x) > 1: print "Warning a large x value has ",\
                       "been discovered with x = ",x
                    x = choose(x < 1, (1,x))
                    x = choose(x > -1, (-1,x))                   
                    sum += x*x
                    Yalpha[a] -= x*x
                    Ybeta[b] -= x*x
                    Yalpha[r] += x*x
                    Ybeta[s] += x*x
                    Ec1 += x*Eio

    #compute the fractional occupations of the occupied orbitals
    for a in aoccs:
        Yalpha[a] = 1 + Yalpha[a]
    for b in boccs:
        Ybeta[b] = 1 + Ybeta[b]
    #for a in xrange(nbf):
        #print "For alpha = ",a,"the fractional occupation is ",Yalpha[a]
        #print "For beta  = ",a,"the fractional occupation is ",Ybeta[a]

    #print the energy and its corrections
    E = energy + Ec1
    print "The total sum of excitations is ",sum
    print "The primary correlation correction is ",Ec1
    print "The total energy is ", E
    return E



                                                         
