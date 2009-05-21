"""\
 Code for restricted open-shell hartree fock programs in PyQuante.

 A good reference for the equations here is 'The Self-Consistent Field
 Equations for Generalized Valence Bond and Open-Shell Hartree-Fock
 Wave Functions', F. W. Bobrowicz and W. A. Goddard, III. in 'Methods
 of Electronic Structure Theory', H. F. Schaefer, III, ed., Plenum
 Publishing Company, 1977.

 This program is part of the PyQuante quantum chemistry program suite.

 Status: Closed shell cases work with the open shell code. 

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
from PyQuante.Ints import getbasis,getints,getJ,getK,get2JmK
from PyQuante.LA2 import mkdens,geigh,trace2,simx
from PyQuante.NumWrap import zeros,take,transpose,matrixmultiply,eigh,dot
from PyQuante.NumWrap import identity
from math import sqrt

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
    # GVB2P5 did this a little more efficiently; they stored
    # 2J-K for the core, then J,K for each open shell. Didn't
    # seem worth it here, so I'm jst storing J,K separately
    Hs = []
    for D in Ds:
        Hs.append(getJ(Ints,D))
        Hs.append(getK(Ints,D))
    return Hs

def get_orbs_in_shell(ish,noccsh,norb):
    # Construct the list of orbitals that must be
    # considered in the active space for the ith shell
    nocc = sum(noccsh)
    vstart,vend = nocc,norb
    istart = sum(noccsh[:ish])
    iend = sum(noccsh[:ish+1])
    return range(istart,iend)+range(vstart,vend)

def get_open_shell_fock(ish,nsh,f,a,b,h,Hs,**kwargs):
    nof = kwargs.get('nof',False)
    # Form the Fock matrix for shell ish:
    if nof:
        F = h
    else:
        F = f[ish]*h
    for jsh in range(nsh):
        if nof:
            F += a[ish,jsh]*Hs[2*jsh]/f[ish]+b[ish,jsh]*Hs[2*jsh+1]/f[ish]
        else:
            F += a[ish,jsh]*Hs[2*jsh]+b[ish,jsh]*Hs[2*jsh+1]
    return F

def update_orbe(orbs_in_shell,orbe,mo_orbe):
    for i,iorb in enumerate(orbs_in_shell):
        orbe[iorb] = mo_orbe[i]
    return

def update_orbs(orbs_in_shell,orbs,mo_orbs):
    nbf,nmo = orbs.shape
    # Map the orbs back to the occupied space
    for i,iorb in enumerate(orbs_in_shell):
        vec = zeros(nbf,'d')
        for j,cj in enumerate(mo_orbs[i,:]):
            vec += cj*orbs[:,orbs_in_shell[j]]
        orbs[:,iorb] = vec
    return 
        
def ocbse(orbs,h,Hs,f,a,b,noccsh):
    # Need to write this so that we don't need the orbs 3 times!
    nsh = len(noccsh)
    nbf = norb = h.shape[0]
    orbe = zeros(norb,'d')
    for ish in range(nsh):
        orbs_in_shell = get_orbs_in_shell(ish,noccsh,norb)
        F = get_open_shell_fock(ish,nsh,f,a,b,h,Hs)
        # form the orbital space of all of the orbs in ish plus the virts
        T = take(orbs,orbs_in_shell,1)
        # Transform to MO space
        Fmo = simx(F,T)
        mo_orbe,mo_orbs = eigh(Fmo)
        # Insert orbital energies into the right place
        update_orbe(orbs_in_shell,orbe,mo_orbe)
        update_orbs(orbs_in_shell,orbs,mo_orbs)
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

def get_fab(nclosed,nopen):
    f = []
    iopen_start = 0
    if nclosed:
        f.append(1.0)
        iopen_start = 1
    if nopen:
        f.append(0.5)
    nsh = len(f)
    a = zeros((nsh,nsh),'d')
    b = zeros((nsh,nsh),'d')

    for i in range(nsh):
        for j in range(nsh):
            a[i,j] = 2.*f[i]*f[j]
            b[i,j] = -f[i]*f[j]

    if nopen == 1:
        a[iopen_start,iopen_start] = 0
        b[iopen_start,iopen_start] = 0
    elif nopen > 1:
        b[iopen_start,iopen_start] = -0.5

    return f,a,b

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
    if not orbs: orbe,orbs = geigh(h,S)

    nclosed,nopen = atoms.get_closedopen()
    nocc = nopen+nclosed
    if not noccsh: noccsh = get_noccsh(nclosed,nopen)
    nsh = len(noccsh)
    nbf = norb = len(bfs)
    if not f:
        f,a,b = get_fab(nclosed,nopen)

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
        Smo = simx(S,orbs)
        print "Smo: \n",Smo[:4,:4]
        orbe,orbs = ocbse(orbs,h,Hs,f,a,b,noccsh)
        orthogonalize(orbs,S)
        Smo = simx(S,orbs)
        print "Smo: \n",Smo[:4,:4]
        # Compute the energy
        eone = 0
        for ish in range(nsh): eone += trace2(Ds[ish],h)
        energy = enuke+eone+sum(orbe[:nocc])
        print "energy = ",energy
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    return energy,orbe,orbs

def orthogonalize(orbs,S):
    nbf,norb = orbs.shape
    for i in range(norb):
        for j in range(i):
            Sij = dot(orbs[:,j],dot(S,orbs[:,i]))
            orbs[:,i] -= Sij*orbs[:,j]
        Sii = dot(orbs[:,i],dot(S,orbs[:,i]))
        orbs[:,i] /= sqrt(Sii)
    return 

if __name__ == '__main__': 
    from PyQuante.Molecule import Molecule
    from PyQuante.hartree_fock import hf
    h = Molecule('H',[(1,(0,0,0))],multiplicity=2)
    he = Molecule('He',[(2,(0,0,0))])
    li = Molecule('Li',[(3,(0,0,0))],multiplicity=2)
    be = Molecule('Be',[(4,(0,0,0))],multiplicity=3)
    mol = li
    print "HF results (for comparison)"
    hfen,hforbe,hforbs = hf(mol,verbose=True)
    print "RHF Energy = ",hfen
    print "RHF spectrum: ",hforbe
    energy,orbe,orbs = rohf(mol)
    print "ROHF Energy = ",energy
    print "ROHF spectrum: ",orbe
    
