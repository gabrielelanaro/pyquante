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
from LA2 import geigh,mkdens,trace2
from Ints import get2JmK,getbasis,getints,getJ,getK
from Convergence import DIIS
import logging

logger = logging.getLogger("pyquante") # Hack!!!

from math import sqrt,pow,fabs
from PyQuante.cints import dist

from NumWrap import array2string, dot, transpose


def get_fock(D,Ints,h):
    "Return the HF Fock matrix"
    return h + get2JmK(Ints,D)

def get_energy(h,F,D,enuke=0.,**opts):
    "Form the total energy of the closed-shell wave function."
    eone = trace2(D,h)
    etwo = trace2(D,F)
    #print "eone,etwo,enuke",(2*eone),(etwo-eone),enuke
    return eone+etwo+enuke

def get_guess(h,S):
    "Form an initial guess from the one-electron Hamiltonian"
    evals,evecs = geigh(h,S)
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
    
    logger.info("RHF calculation on %s" % atoms.name)

    bfs = opts.get('bfs',None)
    if not bfs:
        basis_data = opts.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    nclosed,nopen = atoms.get_closedopen()
    nocc = nclosed
    assert(nopen == 0), "SCF currently only works for closed-shell systems"
    
    logger.info("Nbf = %d" % len(bfs))
    logger.info("Nclosed = %d" % nclosed)

    integrals = opts.get('integrals', None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()

    orbs = opts.get('orbs',None)
    if orbs is None: orbe,orbs = geigh(h,S)

    enuke = atoms.get_enuke()
    eold = 0.

    
    if DoAveraging:
        logger.info("Using DIIS averaging")
        avg = DIIS(S)
    logging.debug("Optimization of HF orbitals")
    for i in xrange(MaxIter):
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
        orbe,orbs = geigh(F,S)
        energy = get_energy(h,F,D,enuke)
        if ETemp:
            energy += entropy
        logging.debug("%d %f" % (i,energy))
        if abs(energy-eold) < ConvCriteria: break
        logger.info("Iteration: %d    Energy: %f    EnergyVar: %f"%(i,energy,abs(energy-eold)))
        eold = energy
    if i < MaxIter:
        logger.info("PyQuante converged in %d iterations" % i)
    else:
        logger.warning("PyQuante failed to converge after %d iterations"
                            % MaxIter)
    logger.info("Final HF energy for system %s is %f" % (atoms.name,energy))
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
    MaxIter = opts.get('MaxIter',40)
    DoAveraging = opts.get('DoAveraging',True)
    averaging = opts.get('averaging',0.5)
    ETemp = opts.get('ETemp',False)
    verbose = opts.get('verbose',False)

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
    if orbs!=None:
        #orbsa = orbsb = orbs
        orbsa = orbs[0]
        orbsb = orbs[1]
    else:
        orbe,orbs = geigh(h,S)
        orbea = orbeb = orbe
        orbsa = orbsb = orbs

    enuke = atoms.get_enuke()
    eold = 0.

    logger.info("UHF calculation on %s" % atoms.name)
    logger.info("Nbf = %d" % len(bfs))
    logger.info("Nalpha = %d" % nalpha)
    logger.info("Nbeta = %d" % nbeta)
    logger.info("Averaging = %s" % DoAveraging)
    logging.debug("Optimization of HF orbitals")
    for i in xrange(MaxIter):
        if verbose: print "SCF Iteration:",i,"Starting Energy:",eold
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
        orbea,orbsa = geigh(Fa,S)
        orbeb,orbsb = geigh(Fb,S)
        energya = get_energy(h,Fa,Da)
        energyb = get_energy(h,Fb,Db)
        energy = (energya+energyb)/2+enuke
        Dab = Da+Db
        Eone = trace2(Dab,h)
        Ej = 0.5*trace2(Dab,Ja+Jb)
        Ek = -0.5*(trace2(Da,Ka)+trace2(Db,Kb))
        if ETemp: energy += entropy
        logger.debug("%d %f %f %f %f" % (i,energy,Eone,Ej,Ek))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logger.info("Final UHF energy for system %s is %f" % (atoms.name,energy))
    return energy,(orbea,orbeb),(orbsa,orbsb)


#added by Hatem H Helal hhh23@cam.ac.uk
def mk_auger_dens(c, occ):
    "Forms a density matrix from a coef matrix c and occupations in occ"
    #count how many states we were given
    nstates = occ.shape[0]
    D = 0.0
    for i in xrange(nstates):
        D += occ[i]*dot( c[:,i:i+1], transpose(c[:,i:i+1]))
    #pad_out(D)
    return D

def uhf_fixed_occ(atoms,occa, occb,**opts):
    """\
    occa and occb represent the orbital occupation arrays for 
    the calculating spin orbitals with holes
    
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

    from biorthogonal import biorthogonalize,pad_out

    ConvCriteria = opts.get('ConvCriteria',1e-5)
    MaxIter = opts.get('MaxIter',40)
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

    orbsa = opts.get('orbsa',None)
    orbsb = opts.get('orbsb',None)
    if (orbsa!=None and orbsb!=None):
        orbsa = orbsa
        orbsb = orbsb
    else:
        orbe,orbs = geigh(h,S)
        orbea = orbeb = orbe
        orbsa = orbsb = orbs
    
    #print "A Trial Orbital Energies:\n", orbea

    print "A Trial Orbitals:\n"
    pad_out(orbsa)

    #print "B Trial Orbital Energies:\n",orbeb

    print "B Trial Orbitals:\n"
    pad_out(orbsb)
    
    enuke = atoms.get_enuke()
    eold = 0.

    for i in xrange(MaxIter):
        print "SCF Iteration:",i,"Starting Energy:",eold
        #save the starting orbitals
        oldorbs_a=orbsa
        oldorbs_b=orbsb

        Da = mk_auger_dens(orbsa,occa)
        Db = mk_auger_dens(orbsb,occb)
        #Da_std = mkdens(orbsa,0,nalpha)
        #Db_std = mkdens(orbsb,0,nbeta)
        #pad_out(Da - Da_std ) #use to test mk_aug_dens with ground state occupations
        #pad_out(Db - Db_std )
        
        
        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)
        Ka = getK(Ints,Da)
        Kb = getK(Ints,Db)
        Fa = h+Ja+Jb-Ka
        Fb = h+Ja+Jb-Kb

        orbea,orbsa = geigh(Fa,S)
        orbeb,orbsb = geigh(Fb,S)
        
        #save the new orbitals
        neworbs_a=orbsa
        neworbs_b=orbsb
        
        #now we biorthogonalize the new orbitals to the old ones
        #to setup occupation arrays for the next scf cycle
        orbsa = biorthogonalize(neworbs_a,oldorbs_a,S,nalpha,occa)
        orbsb = biorthogonalize(neworbs_b,oldorbs_b,S,nbeta,occb)
        
        energya = get_energy(h,Fa,Da)
        energyb = get_energy(h,Fb,Db)
        energy = (energya+energyb)/2+enuke
        Dab = Da+Db
        Eone = trace2(Dab,h)
        Ej = 0.5*trace2(Dab,Ja+Jb)
        Ek = -0.5*(trace2(Da,Ka)+trace2(Db,Kb))
        
        #print "%d %f %f %f %f" % (i,energy,Eone,Ej,Ek)
        
        logging.debug("%d %f %f %f %f" % (i,energy,Eone,Ej,Ek))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
        if i==(MaxIter-1):
            print "Warning: Reached maximum number of SCF cycles may want to rerun calculation with more SCF cycles"
    logger.info("Final UHF energy for system %s is %f" % (atoms.name,energy))
    return energy,(orbea,orbeb),(orbsa,orbsb)


if __name__ == '__main__':
    from Molecule import Molecule
    h2 = Molecule('H2',atomlist=[(1,(1.,0,0)),(1,(-1.,0,0))])
    en = rhf(h2)
    print en
