#!/usr/bin/env python
"""\
 hartree_fock.py Basic routines for HF programs in the PyQuante framework

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

import string,sys,time

from fermi_dirac import get_efermi, get_fermi_occs,mkdens_occs,get_entropy
from LA2 import GHeigenvectors,mkdens,TraceProperty
from Ints import get2JmK,getbasis,getints,getJ,getK
from Convergence import DIIS
from PyQuante import logging

from math import sqrt,pow
from PyQuante.cints import dist

def get_fock(D,Ints,h):
    "Return the HF Fock matrix"
    return h + get2JmK(Ints,D)

def get_energy(h,F,D,enuke=0.,**opts):
    "Form the total energy of the closed-shell wave function."
    eone = TraceProperty(D,h)
    etwo = TraceProperty(D,F)
    return eone+etwo+enuke

def get_guess(h,S):
    "Form an initial guess from the one-electron Hamiltonian"
    evals,evecs = GHeigenvectors(h,S)
    return evecs

def get_nel(atoms,charge=0):
    print "Warning, hartree_fock.get_nel deprecated"
    print "Use the Molecular instance function"
    return atoms.get_nel(charge)

def get_enuke(atoms):
    print "Warning, hartree_fock.get_enuke deprecated"
    print "Use the Molecular instance function"
    return atoms.get_enuke()

def hf(atoms,**opts):
    """\
    Experimental inteface to rhf/uhf routines. See those
    routines for more information.
    """
    # Ultimately this will check spin_type to determine whether
    # to call uhf or rohf for open shell
    nclosed,nopen = atoms.get_closedopen()
    if nopen: return uhf(atoms,**opts)
    return rhf(atoms,**opts)

def rhf(atoms,**opts):
    """\
    rhf(atoms,**opts) - Closed-shell HF driving routine
    
    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    ConvCriteria  1e-4    Convergence Criteria
    MaxIter       20      Maximum SCF iterations
    DoAveraging   True    Use DIIS for accelerated convergence (default)
                  False   No convergence acceleration
    ETemp         False   Use ETemp value for finite temperature DFT (default)
                  float   Use (float) for the electron temperature
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not none, the guess orbitals
    """
    ConvCriteria = opts.get('ConvCriteria',1e-4)
    MaxIter = opts.get('MaxIter',20)
    DoAveraging = opts.get('DoAveraging',False)
    ETemp = opts.get('ETemp',False)

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
    if orbs is None: orbe,orbs = GHeigenvectors(h,S)

    nclosed,nopen = atoms.get_closedopen()
    nocc = nclosed
    assert(nopen == 0), "SCF currently only works for closed-shell systems"
    enuke = atoms.get_enuke()
    eold = 0.

    logging.info("RHF calculation on %s" % atoms.name)
    logging.info("Nbf = %d" % len(bfs))
    logging.info("Nclosed = %d" % nclosed)

    if DoAveraging:
        logging.info("Using DIIS averaging")
        avg = DIIS(S)
    logging.debug("Optimization of HF orbitals")
    for i in range(MaxIter):
        if ETemp:
            efermi = get_efermi(nel,orbe,ETemp)
            occs = get_fermi_occs(efermi,orbe,ETemp)
            D = mkdens_occs(orbs,occs)
            entropy = get_entropy(occs,ETemp)
        else:
            D = mkdens(orbs,0,nocc)
        G = get2JmK(Ints,D)
        F = h+G
        if DoAveraging: F = avg.getF(F,D)
        orbe,orbs = GHeigenvectors(F,S)
        energy = get_energy(h,F,D,enuke)
        if ETemp:
            energy += entropy
        logging.debug("%d %f" % (i,energy))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final HF energy for system %s is %f" % (atoms.name,energy))
    return energy,orbe,orbs

def uhf(atoms,**opts):
    """\
    uhf(atoms,**opts) - Unrestriced Open Shell Hartree Fock
    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    ConvCriteria  1e-4    Convergence Criteria
    MaxIter       20      Maximum SCF iterations
    DoAveraging   True    Use DIIS averaging for convergence acceleration
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not None, the guess orbitals
    """
    ConvCriteria = opts.get('ConvCriteria',1e-5)
    MaxIter = opts.get('MaxIter',20)
    DoAveraging = opts.get('DoAveraging',True)
    averaging = opts.get('averaging',0.5)
    ETemp = opts.get('ETemp',False)

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

    nalpha,nbeta = atoms.get_alphabeta() #pass in opts for multiplicity
    S,h,Ints = getints(bfs,atoms)

    orbs = opts.get('orbs',None)
    if not orbs:
        orbe,orbs = GHeigenvectors(h,S)
        orbea = orbeb = orbe
    orbsa = orbsb = orbs
    
    enuke = atoms.get_enuke()
    eold = 0.

    logging.info("UHF calculation on %s" % atoms.name)
    logging.info("Nbf = %d" % len(bfs))
    logging.info("Nalpha = %d" % nalpha)
    logging.info("Nbeta = %d" % nbeta)
    logging.info("Averaging = %s" % DoAveraging)
    logging.debug("Optimization of HF orbitals")
    for i in range(MaxIter):
        if ETemp:
            # We have to multiply nalpha and nbeta by 2
            #  to get the Fermi energies to come out correct:
            efermia = get_efermi(2.0*nalpha,orbea,ETemp)
            occsa = get_fermi_occs(efermia,orbea,ETemp)
            #print "occsa = ",occsa
            Da = mkdens_occs(orbsa,occsa)
            efermib = get_efermi(2.0*nbeta,orbeb,ETemp)
            occsb = get_fermi_occs(efermib,orbeb,ETemp)
            #print "occsb = ",occsb
            Db = mkdens_occs(orbsb,occsb)
            entropy = 0.5*(get_entropy(occsa,ETemp)+get_entropy(occsb,ETemp))
        else:
            Da = mkdens(orbsa,0,nalpha)
            Db = mkdens(orbsb,0,nbeta)
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
        orbea,orbsa = GHeigenvectors(Fa,S)
        orbeb,orbsb = GHeigenvectors(Fb,S)
        energya = get_energy(h,Fa,Da)
        energyb = get_energy(h,Fb,Db)
        energy = (energya+energyb)/2+enuke
        if ETemp: energy += entropy
        logging.debug("%d %f" % (i,energy))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final UHF energy for system %s is %f" % (atoms.name,energy))
    return energy,(orbea,orbeb),(orbsa,orbsb)

if __name__ == '__main__':
    from Molecule import Molecule
    h2 = Molecule('H2',atomlist=[(1,(1.,0,0)),(1,(-1.,0,0))])
    en = rhf(h2)
    print en
