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
from PyQuante.NumWrap import zeros,transpose,matrixmultiply,eigh,dot
from PyQuante.NumWrap import identity,take
from PyQuante.hartree_fock import get_energy
from math import sqrt

def get_os_dens(orbs,f,noccsh):
    istart = iend = 0
    nsh = len(f)
    Ds = [ ]
    assert len(f) == len(noccsh)
    for ish in xrange(nsh):
        iend += noccsh[ish]
        Ds.append(mkdens(orbs,istart,iend))
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

def get_os_fock(ish,nsh,f,a,b,h,Hs,**kwargs):
    nof = kwargs.get('nof',False)
    # Form the Fock matrix for shell ish:
    if nof:
        F = h
    else:
        F = f[ish]*h
    for jsh in xrange(nsh):
        if nof:
            F += a[ish,jsh]*Hs[2*jsh]/f[ish]+b[ish,jsh]*Hs[2*jsh+1]/f[ish]
        else:
            F += a[ish,jsh]*Hs[2*jsh]+b[ish,jsh]*Hs[2*jsh+1]
    return F

def update_orbe(orbs_in_shell,orbe,mo_orbe):
    for i,iorb in enumerate(orbs_in_shell):
        orbe[iorb] = mo_orbe[i]
    return

def update_orbs(orbs_in_shell,orbs,new_orbs):
    for i,iorb in enumerate(orbs_in_shell):
        orbs[:,iorb] = new_orbs[:,i]
    return 
        
def ocbse(orbs,h,Hs,f,a,b,noccsh):
    # Need to write this so that we don't need the orbs 3 times!
    nsh = len(noccsh)
    nbf = norb = h.shape[0]
    orbe = zeros(norb,'d')
    for ish in xrange(nsh):
        orbs_in_shell = get_orbs_in_shell(ish,noccsh,norb)
        F = get_os_fock(ish,nsh,f,a,b,h,Hs)
        # form the orbital space of all of the orbs in ish plus the virts
        T = orbs.take(orbs_in_shell,1)
        #print "take worked? ",(T==get_orbs(orbs,orbs_in_shell)).all()
        # Transform to MO space
        Fmo = ao2mo(F,T)
        mo_orbe,mo_orbs = eigh(Fmo)
        T = matrixmultiply(T,mo_orbs)
        # Insert orbital energies into the right place
        update_orbe(orbs_in_shell,orbe,mo_orbe)
        update_orbs(orbs_in_shell,orbs,T)
    return orbe,orbs

def get_orbs(orbs,orbs_in_shell):
    "This should do the same thing as take(orbs,orbs_in_shell,1)"
    A = zeros((orbs.shape[0],len(orbs_in_shell)),'d')
    for i,iorb in enumerate(orbs_in_shell):
        A[:,i] = orbs[:,iorb]
    return A

def rotion(orbs,h,Hs,f,a,b,noccsh):
    nsh = len(noccsh)
    nocc = sum(noccsh)
    if nsh == 1: return orbs # No effect for closed shell systems
    rot = get_rot(h,Hs,f,a,b,noccsh)
    print "Rotation matrix:\n",rot
    erot = expmat(rot)
    print "Exp rotation matrix:\n",erot
    T = matrixmultiply(orbs[:,:nocc],erot)
    orbs[:,:nocc] = T
    return orbs

def expmats(A):
    # For testing agains scipy
    from scipy.linalg.matfuncs import expm
    return expm(A)

def expmat(A,**kwargs):
    nmax = kwargs.get('nmax',12)
    cut = kwargs.get('cut',1e-8)
    E = identity(A.shape[0],'d')
    D = E
    for i in xrange(1,nmax):
        D = matrixmultiply(D,A)/i
        E += D
        maxel = D.max()
        if abs(maxel) < cut:
            break
    else:
        print "Warning: expmat unconverged after %d iters: %g" % (nmax,maxel)
    return E

def get_sh(i,noccsh):
    nsh = len(noccsh)
    isum = 0
    for ish in xrange(nsh):
        isum += noccsh[ish]
        if i < isum:
            return i
    return None

def get_rot(h,Hs,f,a,b,noccsh):
    nocc = sum(noccsh)
    nsh = len(noccsh)
    rot = zeros((nocc,nocc),'d')
    for i in xrange(nocc):
        ish = get_sh(i,noccsh)
        for j in xrange(nocc):
            jsh = get_sh(j,noccsh)
            if jsh == ish: continue
            Wij = -0.5*(h[i,j]+Hs[0][i,j])
            Wii = -0.5*(h[i,i]+Hs[0][i,i])
            Wjj = -0.5*(h[j,j]+Hs[0][j,j])
            for k in xrange(nsh):
                Wij = Wij - 0.5*Hs[2*i+1][i,j]
                Wii = Wij - 0.5*Hs[2*i+1][i,i]
                Wjj = Wij - 0.5*Hs[2*i+1][j,j]
            Jij = Hs[2*jsh][i,i]
            Kij = Hs[2*jsh+1][i,i]
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

    for i in xrange(nsh):
        for j in xrange(nsh):
            a[i,j] = 2.*f[i]*f[j]
            b[i,j] = -f[i]*f[j]

    if nopen == 1:
        a[iopen_start,iopen_start] = 0
        b[iopen_start,iopen_start] = 0
    elif nopen > 1:
        b[iopen_start,iopen_start] = -0.5

    return f,a,b

def rohf_wag(atoms,noccsh=None,f=None,a=None,b=None,**kwargs):
    """\
    rohf(atoms,noccsh=None,f=None,a=None,b=None,**kwargs):
        Restricted open shell HF driving routine

    atoms      A Molecule object containing the system of interest
    """
    ConvCriteria = kwargs.get('ConvCriteria',1e-4)
    MaxIter = kwargs.get('MaxIter',25)
    DoAveraging = kwargs.get('DoAveraging',False)
    verbose = kwargs.get('verbose',True)

    bfs = kwargs.get('bfs',None)
    if not bfs:
        basis_data = kwargs.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = kwargs.get('integrals', None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()

    orbs = kwargs.get('orbs',None)
    if orbs is None:
        orbe,orbs = geigh(h,S)

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
        for i in xrange(nsh):
            for j in xrange(i+1):
                print a[i,j],
            print
        print "b_ij: "
        for i in xrange(nsh):
            for j in xrange(i+1):
                print b[i,j],
            print
    enuke = atoms.get_enuke()
    energy = eold = 0.
    for i in xrange(MaxIter):
        Ds = get_os_dens(orbs,f,noccsh)
        Hs = get_os_hams(Ints,Ds)
        orbs = rotion(orbs,h,Hs,f,a,b,noccsh)
        orbe,orbs = ocbse(orbs,h,Hs,f,a,b,noccsh)
        orthogonalize(orbs,S)
        # Compute the energy
        eone = sum(f[ish]*trace2(Ds[ish],h) for ish in xrange(nsh))
        energy = enuke+eone+sum(orbe[:nocc])
        print energy,eone
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    return energy,orbe,orbs

def rohf(atoms,**opts):
    """\
    rohf(atoms,**opts) - Restriced Open Shell Hartree Fock
    atoms       A Molecule object containing the molecule
    """

    ConvCriteria = opts.get('ConvCriteria',1e-5)
    MaxIter = opts.get('MaxIter',40)
    DoAveraging = opts.get('DoAveraging',True)
    averaging = opts.get('averaging',0.95)
    verbose = opts.get('verbose',True)

    bfs = opts.get('bfs',None)
    if not bfs:
        basis_data = opts.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)
    nbf = len(bfs)

    integrals = opts.get('integrals', None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()

    nalpha,nbeta = atoms.get_alphabeta()
    S,h,Ints = getints(bfs,atoms)

    orbs = opts.get('orbs',None)
    if orbs is None:
        orbe,orbs = geigh(h,S)
    norbs = nbf

    enuke = atoms.get_enuke()
    eold = 0.

    if verbose: print "ROHF calculation on %s" % atoms.name
    if verbose: print "Nbf = %d" % nbf
    if verbose: print "Nalpha = %d" % nalpha
    if verbose: print "Nbeta = %d" % nbeta
    if verbose: print "Averaging = %s" % DoAveraging
    print "Optimization of HF orbitals"

    for i in xrange(MaxIter):
        if verbose: print "SCF Iteration:",i,"Starting Energy:",eold
        Da = mkdens(orbs,0,nalpha)
        Db = mkdens(orbs,0,nbeta)
        if DoAveraging:
            if i: 
                Da = averaging*Da + (1-averaging)*Da0
                Db = averaging*Db + (1-averaging)*Db0
            Da0 = Da
            Db0 = Db

        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)
        Ka = getK(Ints,Da)
        Kb = getK(Ints,Db)

        Fa = h+Ja+Jb-Ka
        Fb = h+Ja+Jb-Kb
        energya = get_energy(h,Fa,Da)
        energyb = get_energy(h,Fb,Db)
        eone = (trace2(Da,h) + trace2(Db,h))/2
        etwo = (trace2(Da,Fa) + trace2(Db,Fb))/2
        energy = (energya+energyb)/2 + enuke
        print i,energy,eone,etwo,enuke
        if abs(energy-eold) < ConvCriteria: break
        eold = energy

        Fa = ao2mo(Fa,orbs)
        Fb = ao2mo(Fb,orbs)

        # Building the approximate Fock matrices in the MO basis
        F = 0.5*(Fa+Fb)
        K = Fb-Fa

        # The Fock matrix now looks like
        #      F-K    |  F + K/2  |    F
        #   ---------------------------------
        #    F + K/2  |     F     |  F - K/2
        #   ---------------------------------
        #       F     |  F - K/2  |  F + K

        # Make explicit slice objects to simplify this
        do = slice(0,nbeta)
        so = slice(nbeta,nalpha)
        uo = slice(nalpha,norbs)
        F[do,do] -= K[do,do]
        F[uo,uo] += K[uo,uo]
        F[do,so] += 0.5*K[do,so]
        F[so,do] += 0.5*K[so,do]
        F[so,uo] -= 0.5*K[so,uo]
        F[uo,so] -= 0.5*K[uo,so]

        orbe,mo_orbs = eigh(F)
        orbs = matrixmultiply(orbs,mo_orbs)
        
    if verbose:
        print "Final ROHF energy for system %s is %f" % (atoms.name,energy)
    return energy,orbe,orbs

def ao2mo(M,C): return simx(M,C)
def mo2ao(M,C,S):
    SC = matrixmultiply(S,C)
    return simx(M,SC,'t')

def symmetrize(A): return (A+A.T)/2

def printmat(mat,name='mat',**kwargs):
    istart = kwargs.get('istart',0)
    istop = kwargs.get('istop',4)
    jstart = kwargs.get('jstart',0)
    jstop = kwargs.get('jstop',4)
    suppress = kwargs.get('suppress',True)
    print name,'\n',mat[istart:istop,jstart:jstop]
    return

def orthogonalize(orbs,S):
    nbf,norb = orbs.shape
    Smax = 0
    for i in xrange(norb):
        for j in xrange(i):
            Sij = dot(orbs[:,j],dot(S,orbs[:,i]))
            Smax = max(Smax,abs(Sij))
            orbs[:,i] -= Sij*orbs[:,j]
        Sii = dot(orbs[:,i],dot(S,orbs[:,i]))
        orbs[:,i] /= sqrt(Sii)
    print "Max orthogonalized element = ",Smax
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
    energy,orbe,orbs = rohf(mol)
    print "ROHF Energy = ",energy
    print "ROHF spectrum: \n",orbe
    #energy,orbe,orbs = rohf_wag(mol,orbs=orbs)
    energy,orbe,orbs = rohf_wag(mol)
    print "ROHF Energy = ",energy
    print "ROHF spectrum: \n",orbe
    
