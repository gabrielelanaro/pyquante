"""\
 Ints.py Basic routines for integrals in the PyQuante framework

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from CGBF import CGBF,coulomb
#from contracted_gto import coulomb
from NumWrap import zeros,dot,reshape
from PyQuante.cints import ijkl2intindex as intindex
from PyQuante.Basis.Tools import get_basis_data
import settings
import logging

logger = logging.getLogger("pyquante")

sym2powerlist = {
    'S' : [(0,0,0)],
    'P' : [(1,0,0),(0,1,0),(0,0,1)],
    'D' : [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)],
    'F' : [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),
           (0,3,0),(0,2,1),(0,1,2), (0,0,3)]
    }

sorted = True
jints = {}
kints = {}

def getbasis(atoms,basis_data=None,**opts):
    """\
    bfs = getbasis(atoms,basis_data=None)
    
    Given a Molecule object and a basis library, form a basis set
    constructed as a list of CGBF basis functions objects.
    """
    from PyQuante.Basis.basis import BasisSet
    return BasisSet(atoms, basis_data, **opts)
    # Option to omit f basis functions from imported basis sets
    omit_f = opts.get('omit_f',False)
    if not basis_data:
        from PyQuante.Basis.p631ss import basis_data
    elif type(basis_data) == type(''):
        # Assume this is a name of a basis set, e.g. '6-31g**'
        #  and import dynamically
        basis_data = get_basis_data(basis_data)
    bfs = []
    for atom in atoms:
        bs = basis_data[atom.atno]
        for sym,prims in bs:
            if omit_f and sym == "F": continue
            for power in sym2powerlist[sym]:
                bf = CGBF(atom.pos(),power,atom.atid)
                for expnt,coef in prims:
                    bf.add_primitive(expnt,coef)
                bf.normalize()
                bfs.append(bf)
    return bfs


def getints(bfs,atoms):
    logger.info("Calculating Integrals...")
    S,h = get1ints(bfs,atoms)
    Ints = get2ints(bfs)
    logger.info("Integrals Calculated.")
    return S,h,Ints

def get1ints(bfs,atoms):
    "Form the overlap S and h=t+vN one-electron Hamiltonian matrices"
    nbf = len(bfs)
    S = zeros((nbf,nbf),'d')
    h = zeros((nbf,nbf),'d')

    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(nbf):
            bfj = bfs[j]
            S[i,j] = bfi.overlap(bfj)
            h[i,j] = bfi.kinetic(bfj)
            for atom in atoms:
                h[i,j] = h[i,j] + atom.atno*bfi.nuclear(bfj,atom.pos())
    return S,h

def getT(bfs):
    "Form the kinetic energy matrix"
    nbf = len(bfs)
    T = zeros((nbf,nbf),'d')

    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(nbf):
            bfj = bfs[j]
            T[i,j] = bfi.kinetic(bfj)
    return T

def getS(bfs):
    "Form the overlap matrix"
    nbf = len(bfs)
    S = zeros((nbf,nbf),'d')

    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(nbf):
            bfj = bfs[j]
            S[i,j] = bfi.overlap(bfj)
    return S

def getV(bfs,atoms):
    "Form the nuclear attraction matrix V"
    nbf = len(bfs)
    V = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        bfi = bfs[i]
        for j in xrange(nbf):
            bfj = bfs[j]
            for atom in atoms:
                V[i,j] = V[i,j] + atom.atno*bfi.nuclear(bfj,atom.pos())
    return V

if settings.libint_enabled == True:
    # Libint Integrals
    import numpy as np
    import clibint
    
    def get2ints(basis):
        lenbasis = len(basis.bfs)
        
        Ints = np.zeros((lenbasis**4),dtype=np.float64)

        for i,a in enumerate(basis.shells):
            for j,b in enumerate(basis.shells[:i+1]):
                for k,c in enumerate(basis.shells):
                    for l,d in enumerate(basis.shells[:k+1]):
                        if (i+j)>=(k+l):
                            clibint.shell_compute_eri(a,b,c,d,Ints)
        if sorted:
            sortints(lenbasis,Ints)
        return Ints
else:
    # PyQuante Integrals
    def get2ints(bfs):
        """Store integrals in a long array in the form (ij|kl) (chemists
        notation. We only need i>=j, k>=l, and ij <= kl"""

        from array import array
        nbf = len(bfs)
        totlen = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
        Ints = array('d',[0]*totlen)

        for i in xrange(nbf):
            for j in xrange(i+1):
                ij = i*(i+1)/2+j
                for k in xrange(nbf):
                    for l in xrange(k+1):
                        kl = k*(k+1)/2+l
                        if ij <= kl:
                            Ints[intindex(i,j,k,l)] = coulomb(bfs[i],bfs[j],
                                                              bfs[k],bfs[l])

        if sorted:
            sortints(nbf,Ints)
        return Ints

def sortints(nbf,Ints):
    for i in xrange(nbf):
        for j in xrange(i+1):
            jints[i,j] = fetch_jints(Ints,i,j,nbf)
            kints[i,j] = fetch_kints(Ints,i,j,nbf)
    return

def fetch_jints(Ints,i,j,nbf):
    temp = zeros(nbf*nbf,'d')
    kl = 0
    for k in xrange(nbf):
        for l in xrange(nbf):
            index = intindex(i,j,k,l)
            temp[kl] = Ints[index]
            kl += 1
    return temp

def fetch_kints(Ints,i,j,nbf):
    temp = zeros(nbf*nbf,'d')
    kl = 0
    for k in xrange(nbf):
        for l in xrange(nbf):
            temp[kl] = Ints[intindex(i,k,j,l)]
            kl += 1
    return temp

def getJ(Ints,D):
    "Form the Coulomb operator corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    J = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        for j in xrange(i+1):
            if sorted:
                temp = jints[i,j]
            else:
                temp = fetch_jints(Ints,i,j,nbf)
            J[i,j] = dot(temp,D1d)
            J[j,i] = J[i,j]
    return J

def getK(Ints,D):
    "Form the exchange operator corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    K = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        for j in xrange(i+1):
            if sorted:
                temp = kints[i,j]
            else:
                temp = fetch_kints(Ints,i,j,nbf)
            K[i,j] = dot(temp,D1d)
            K[j,i] = K[i,j]
    return K

def get2JmK(Ints,D):
    "Form the 2J-K integrals corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    G = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        for j in xrange(i+1):
            if sorted:
                temp = 2*jints[i,j]-kints[i,j]
            else:
                temp = 2*fetch_jints(Ints,i,j,nbf)-fetch_kints(Ints,i,j,nbf)
            G[i,j] = dot(temp,D1d)
            G[j,i] = G[i,j]
    return G

