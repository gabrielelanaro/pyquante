"""\
 fermi_dirac.py: Utilities for finite temperature Fermi-Dirac occupations.

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
import sys
from NumWrap import matrixmultiply,transpose
from math import exp,log
from Constants import Kboltz
from LA2 import mkdens
import logging
logger = logging.getLogger("pyquante")

def mkdens_fermi(nel,orbe,orbs,e_temp):
    """
    mkdens_fermi(nel,orbe,orbs,e_temp)

    Create a density matrix from the orbitals, Orbs, and the Fermi-Dirac
    occupations, Occs, derived from the orbital energies, Orbe, given the
    electron temperature, e_temp.

    D = Orbs*Occs(Orbe)Orbs^T

    Arguments:
    nel     Number of electrons in the system
    orbe    The orbital energies
    orbs    The orbitals
    e_temp  The electron temperature
    """
    efermi = get_efermi(nel,orbe,e_temp)
    occs = get_fermi_occs(efermi,orbe,e_temp)
    D = mkdens_occs(orbs,occs)
    entropy = get_entropy(occs,e_temp)
    return D,entropy

def mkdens_occs(c,occs,**opts):
    "Density matrix from a set of occupations (e.g. from FD expression)."
    tol = opts.get('tol',1e-5)
    verbose = opts.get('verbose',False)
    # Determine how many orbs have occupations greater than 0
    norb = 0
    for fi in occs:
        if fi < tol: break
        norb += 1
    if verbose:
        print "mkdens_occs: %d occupied orbitals found" % norb
    # Determine how many doubly occupied orbitals we have
    nclosed = 0
    for i in xrange(norb):
        if abs(1.-occs[i]) > tol: break
        nclosed += 1
    if verbose:
        print "mkdens_occs: %d closed-shell orbitals found" % nclosed
    D = mkdens(c,0,nclosed)
    for i in xrange(nclosed,norb):
        D = D + occs[i]*matrixmultiply(c[:,i:i+1],transpose(c[:,i:i+1]))
    return D
    
def get_fermi_occ(efermi,en,temp):
    kT = Kboltz*temp
    x = (en-efermi)/kT
    if x < -50.: return 1.
    elif x > 50.: return 0
    return 1/(1+exp(x))

def get_entropy(occs,temp):
    kT = Kboltz*temp
    entropy = 0
    for fi in occs:
        if abs(fi) < 1e-10: break # stop summing when occs get small
        if fi > 1e-10:
            entropy += kT*fi*log(fi)
        if (1-fi) > 1e-10:
            entropy += kT*(1.-fi)*log(1.-fi)
    return entropy
    

def get_fermi_occs(efermi,orbe,temp):
    occs = []
    for en in orbe:
        occs.append(get_fermi_occ(efermi,en,temp))
    return occs

def get_t0_occs(nel,nbf):
    occs = [0]*nbf
    nc,no = divmod(nel,2)
    for i in xrange(nc): occs[i] = 1.
    for i in xrange(nc,nc+no): occs[i] = 0.5
    return occs

def get_efermi(nel,orbe,temp,**opts):
    "Bisection method to get Fermi energy from Fermi-Dirac dist"
    tol = opts.get('tol',1e-9)
    verbose = opts.get('verbose',True)

    elow,ehigh = orbe[0]-100.,orbe[-1]
    nlow = 2*sum(get_fermi_occs(elow,orbe,temp))
    nhigh = 2*sum(get_fermi_occs(ehigh,orbe,temp))

    if nlow > nel:
        logger.error("elow incorrect %f -> %f " % (elow,nlow))
        raise Exception("elow incorrect %f -> %f " % (elow,nlow))
    if nhigh < nel:
        logger.error("ehigh incorrect %f -> %f " % (ehigh,nhigh))
        raise Exception("ehigh incorrect %f -> %f " % (ehigh,nhigh))

    for i in xrange(100):
        efermi = (elow+ehigh)/2
        n = 2*sum(get_fermi_occs(efermi,orbe,temp))
        if abs(n-nel) < tol:
            break
        elif n < nel:
            elow = efermi
        else:
            ehigh = efermi
    else:
        print "get_fd_occs: Too many iterations"
    return efermi

