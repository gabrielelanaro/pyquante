#!/usr/bin/env python

"""\
 Code for restricted open-shell hartree fock programs in PyQuante.

 A good reference for the equations here is 'The Self-Consistent Field
 Equations for Generalized Valence Bond and Open-Shell Hartree-Fock
 Wave Functions', F. W. Bobrowicz and W. A. Goddard, III. in 'Methods
 of Electronic Structure Theory', H. F. Schaefer, III, ed., Plenum
 Publishing Company, 1977.

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
from PyQuante.Ints import getbasis,getints,getJ,getK,get2JmK
from PyQuante.LA2 import mkdens,GHeigenvectors,SimilarityTransform,\
     TraceProperty,SimilarityTransformT
from NumWrap import zeros,take,transpose,matrixmultiply
from NumWrap import Heigenvectors

def get_os_dens(orbs,f,noccsh):
    istart = iend = 0
    nsh = len(f)
    Ds = [ ]
    assert len(f) == len(noccsh)
    for ish in range(nsh):
        iend += noccsh[ish]
        Ds.append(f[ish]*mkdens(orbs,istart,iend))
        istart = iend
    return Ds

def get_os_hams(Ints,Ds):
    Hs = [get2JmK(Ints,Ds[0])]
    for D in Ds[1:]:
        Hs.append(getJ(Ints,D))
        Hs.append(getK(Ints,D))
    return Hs

def ocbse(orbs,h,Hs,f,a,b,noccsh):
    # Need to write this so that we don't need the orbs 3 times!
    nsh = len(noccsh)
    nocc = sum(noccsh)
    nbf = norb = h.shape[0]
    vstart,vend = (nocc,norb) # range limits for virtual orbs
    nvirt = vend-vstart
    orbs3 = zeros((nbf,nbf),'d')
    orbe = zeros(nbf,'d')
    for ish in range(nsh):
        # Form the range of orbitals for this shell
        istart = sum(noccsh[:ish])
        iend = sum(noccsh[:ish+1])
        # Form the Fock matrix for shell ish:
        F = f[ish]*h
        for jsh in range(nsh): F += a[ish,jsh]*Hs[jsh] + b[ish,jsh]*Hs[jsh]
        # form the orbital space of all of the orbs in ish plus the virts
        orbrange = range(istart,iend)+range(vstart,vend)
        T = take(orbs,orbrange,1)
        # Transform to MO space
        #F = SimilarityTransform(F,T) # SimilarityTransformT??
        FT = matrixmultiply(F,transpose(T))
        F = matrixmultiply(T,FT)
        orbe2,orbs2 = Heigenvectors(F)
        # Insert orbital energies into the right place
        orbe[istart:iend] = orbe2[:noccsh[ish]]
        orbe[vstart:vend] = orbe2[-nvirt:]
        # Map the orbs back to the occupied space
        for i in range(len(orbrange)):
            iorb = orbrange[i]
            coefs = orbs2[i]
            nbfi = len(coefs)
            for j in range(nbfi):
                orbs3[orbrange[i],:] += coefs[j]*orbs[orbrange[j],:]
        # Have to copy orbs3 -> orbs, here?
        orbs = orbs3
    return orbe,orbs

def rotion(orbs,h,Hs,nclosed,nopen):
    nham = (len(Hs)+1)/2
    if nham == 1: return orbs # No effect for closed shell systems
    rot = get_rot(h,Hs)
    erot = expmat(rot)

def get_rot(h,Hs,nclosed,nopen):
    rot = zeros((nclosed+nopen,nclosed+nopen))
    for j in range(nclosed,nclosed+nopen):
        for i in range(nclosed):
            Wij = -0.5*(h[i,j]+Hs[0][i,j])
            Wii = -0.5*(h[i,i]+Hs[0][i,i])
            Wjj = -0.5*(h[j,j]+Hs[0][j,j])
            for k in range(nham):
                Wij = Wij - 0.5*Hs[2*i+1][i,j]
                Wii = Wij - 0.5*Hs[2*i+1][i,i]
                Wjj = Wij - 0.5*Hs[2*i+1][j,j]
            jsh = j-nclosed
            Jij = Hs[2*jsh+1][i,i]
            Kij = Hs[2*jsh+2][i,i]
            gamma = Kij-0.5*(Kij+Jij)
            Xij = -Wij
            Bij = Wii-Wjj+gamma
            if Bij > 0:
                Rij = -Xij/Bij
            else:
                Rij = Xij/Bij
            rot[i,j] = rot[j,i] = Rij
    return rot

def get_noccsh(nclosed,nopen):
    # Get noccsh from open/closed
    noccsh = []
    if nclosed: noccsh.append(nclosed)
    if nopen:   noccsh.append(nopen)
    return noccsh

def get_f(nclosed,nopen):
    # Get f from noccsh
    f = []
    if nclosed: f.append(1.0)
    if nopen:   f.append(0.5)
    return f

def get_a(nsh,f):
    a = zeros((nsh,nsh),'d')
    for i in range(nsh):
        a[i,i] = f[i]
        for j in range(i):
            a[i,j] = 2.*f[i]*f[j]
            a[j,i] = a[i,j]
    return a

def get_b(nsh,f):
    b = zeros((nsh,nsh),'d')
    for i in range(nsh):
        b[i,i] = 0
        for j in range(i):
            if f[i] == 0.5 and f[j] == 0.5:
                b[i,j] = 0.5
            else:
                b[i,j] = -f[i]*f[j]
            b[j,i] = b[i,j]
    return b

def rohf(atoms,noccsh=None,f=None,a=None,b=None,**opts):
    """\
    rohf(atoms,noccsh=None,f=None,a=None,b=None,**opts):
        Restricted open shell HF driving routine

    atoms      A Molecule object containing the system of interest
    """
    ConvCriteria = opts.get('ConvCriteria',1e-4)
    MaxIter = opts.get('MaxIter',20)
    DoAveraging = opts.get('DoAveraging',False)
    verbose = opts.get('verbose',True)

    bfs = opts.get('bfs',None)
    if not bfs:
        basis_data = opts.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = opts.get('integrals', None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()

    orbs = opts.get('orbs',None)
    if not orbs: orbe,orbs = GHeigenvectors(h,S)

    nclosed,nopen = atoms.get_closedopen()
    nocc = nopen+nclosed
    if not noccsh: noccsh = get_noccsh(nclosed,nopen)
    nsh = len(noccsh)
    nbf = norb = len(bfs)
    if not f: f = get_f(nclosed,nopen)
    if not a: a = get_a(nsh,f)
    if not b: b = get_b(nsh,f)

    if verbose:
        print "ROHF calculation"
        print "nsh = ",nsh
        print "noccsh = ",noccsh
        print "f = ",f
        print "a_ij: "
        for i in range(nsh):
            for j in range(i+1):
                print a[i,j],
            print
        print "b_ij: "
        for i in range(nsh):
            for j in range(i+1):
                print b[i,j],
            print
    enuke = atoms.get_enuke()
    energy = eold = 0.
    for i in range(MaxIter):
        Ds = get_os_dens(orbs,f,noccsh)
        Hs = get_os_hams(Ints,Ds)
        #orbs = rotion(orbs,h,Hs,nclosed,nopen)
        orbe,orbs = ocbse(orbs,h,Hs,f,a,b,noccsh)
        # Compute the energy
        eone = 0
        for ish in range(nsh): eone += TraceProperty(Ds[ish],h)
        energy = enuke+eone+sum(orbe[:nocc])
        print "energy = ",energy
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    return energy,orbe,orbs

if __name__ == '__main__': 
    from PyQuante.Molecule import Molecule
    from PyQuante.hartree_fock import hf
    he = Molecule('He',[(2,(0,0,0))])
    li = Molecule('Li',[(3,(0,0,0))],multiplicity=2)
    be = Molecule('Be',[(4,(0,0,0))],multiplicity=3)
    mol = be
    print "RHF results (for comparison)"
    hfen,hforbe,hforbs = hf(mol,verbose=True)
    print "RHF Energy = ",hfen
    print "RHF spectrum: ",hforbe
    energy,orbe,orbs = rohf(mol)
    print "ROHF Energy = ",energy
    print "ROHF spectrum: ",orbe
    
