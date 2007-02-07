#!/usr/bin/env python
"""\
The SCF.py module is intened to be a unified interface to all of the
routines. Currently there are separate functions for restriced/average
spin SCF.py, and unrestricted/spin polarized USCF.py, but I hope to
unify these at some point in the near future.

The goal in unifying these routines is to standardize the input and
output, and hopefully to lower the barrier for people understanding and
being able to modify the code.

"""

# Todo:
# - Unify SCF/USCF
# - Include MINDO code?


from PyQuante import logging
from PyQuante.Convergence import DIIS
from PyQuante.DFunctionals import XC,need_gradients
from PyQuante.Ints import getbasis,getints,getJ,getK
from PyQuante.LA2 import GHeigenvectors,mkdens,mkdens_spinavg,TraceProperty
from PyQuante.MolecularGrid import MolecularGrid
from PyQuante.Molecule import Molecule
from PyQuante.dft import getXC
from PyQuante.fermi_dirac import mkdens_fermi

def SCF(atoms,**opts):
    """\
    energy, orbe, orbs = SCF(atoms, **opts)
    
    Compute the energy of the molecule in *atoms* using a unified
    interface to the HF and DFT routines in PyQuante.

    Options:      Value   Description
    --------      -----   -----------
    verbose       False   Output terse information to the logger (default)
                  True    Print out additional information 
    etol          1e-4    Energy convergence criteria
    max_iter      20      Maximum SCF iterations
    do_averaging  True    Use DIIS for accelerated convergence (default)
                  False   No convergence acceleration
    etemp         False   Use etemp value for finite temperature DFT (default)
                  float   Use (float) for the electron temperature
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not none, the guess orbitals
    functional    None    Use HF for the calculation (default)
                  SVWN    Use the SVWN (LDA) DFT functional
                  S0      Use the Slater Xalpha DFT functional
                  BLYP    Use the BLYP GGA DFT functional
                  PBE     Use the PBE DFT functional
    grid_nrad     32      Number of radial shells per atom
    grid_fineness 1       Radial shell fineness. 0->coarse, 1->medium, 2->fine
    etarget       None    Target energy for testing purposes
    targettol     1e-4    Tolerance for energy target sucess
    """
    # Get the various options from the kwargs flags
    verbose = opts.get('verbose',False) 
    etol = opts.get('etol',1e-4)
    max_iter = opts.get('max_iter',20)
    do_averaging = opts.get('do_averaging',True)
    etemp = opts.get('etemp',False)
    functional = opts.get('functional',None)
    etarget = opts.get('etarget',None)
    targettol = opts.get('targettol',1e-4)

    if verbose:
        oldloglevel = logging.root.level
        logging.root.setLevel(logging.DEBUG)

    if functional:
        do_dft = True
        opts['do_grad_dens'] = need_gradients[functional]
        method = functional
    else:
        method = "HF"
        do_dft = False

    bfs = opts.get('bfs',None)
    if not bfs:
        basis_data = opts.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = opts.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    Enuke = atoms.get_enuke()

    if do_dft:
        # default medium mesh
        grid_nrad = opts.get('grid_nrad',32)
        grid_fineness = opts.get('grid_fineness',1)

        gr = MolecularGrid(atoms,grid_nrad,grid_fineness,**opts) 
        gr.set_bf_amps(bfs)

        bfgrid = gr.allbfs() # bfs over all grid points

    # It would be nice to have a more intelligent treatment of the guess
    # so that I could pass in a density rather than a set of orbs.
    orbs = opts.get('orbs',None)
    if orbs is None: orbe,orbs = GHeigenvectors(h,S)

    nclosed,nopen = atoms.get_closedopen()

    logging.info("%s calculation on %s" % (method,atoms.name))

    logging.debug("Nbf = %d" % len(bfs))
    logging.debug("Nclosed = %d" % nclosed)
    logging.debug("Nopen = %d" % nopen)

    if nopen:
        if not functional:
            logging.warning("Warning! Running average occupation HF")
        else:
            logging.debug("Using average occupation spin")        

    if etemp:
        logging.info("Electron temperature set to %.0f" % etemp)
        if not functional:
            logging.warning("Warning! Running finite temperature HF")

    eold = 0.
    if do_averaging:
        logging.debug("Using DIIS averaging")
        avg=DIIS(S)

    logging.debug("Beginning SCF Optimization")
    logging.debug("It   Etot       Eone        Ej        Exc         Enuke")
    logging.debug("--   -------    -------     ------    -------     ------")
    for i in range(max_iter):
        if etemp:
            D,entropy = mkdens_fermi(nel,orbe,orbs,etemp)
        else:
            D = mkdens_spinavg(orbs,nclosed,nopen)

        if do_dft:
            gr.setdens(D)

        J = getJ(Ints,D)

        if do_dft:
            Exc,XC = getXC(gr,nel,bfgrid,**opts)
            F = h+2*J+XC
        else:
            K = getK(Ints,D)
            Exc = -TraceProperty(D,K)
            F = h + 2*J-K

        if do_averaging:
            F = avg.getF(F,D)

        orbe,orbs = GHeigenvectors(F,S)
            
        # Compute energy components
        Ej = 2*TraceProperty(D,J)
        Eone = 2*TraceProperty(D,h)
        energy = Eone + Ej + Exc + Enuke
        if etemp: energy += entropy
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
                  (i,energy,Eone,Ej,Exc,Enuke))
        if abs(energy-eold) < etol: break
        eold = energy
    logging.info("Final energy for system %s is %f"
                 % (atoms.name,energy))
    if etarget:
        worked = abs(energy-etarget) < targettol
        if worked:
            logging.info("Target energy met")
        else:
            logging.warning("Warning! energy should have been %f" % etarget)
    logging.info("\n")
    if verbose:
        # Set the logging back to what it was before
        logging.root.setLevel(oldloglevel)
    return energy,orbe,orbs

def USCF(atoms,**opts):
    """\
    energy, orbe, orbs = USCF(atoms, **opts)
    
    Compute the energy of the molecule in *atoms* using a unified
    interface to the HF and DFT routines in PyQuante.

    This routine uses an unrestricted/spin-polarized formalism

    Options:      Value   Description
    --------      -----   -----------
    verbose       False   Output terse information to the logger (default)
                  True    Print out additional information 
    etol          1e-4    Energy convergence criteria
    max_iter      20      Maximum SCF iterations
    do_averaging  True    Use DIIS for accelerated convergence (default)
                  False   No convergence acceleration
    etemp         False   Use etemp value for finite temperature DFT (default)
                  float   Use (float) for the electron temperature
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not none, the guess orbitals
    functional    None    Use HF for the calculation (default)
                  SVWN    Use the SVWN (LDA) DFT functional
                  S0      Use the Slater Xalpha DFT functional
                  BLYP    Use the BLYP GGA DFT functional
                  PBE     Use the PBE DFT functional
    grid_nrad     32      Number of radial shells per atom
    grid_fineness 1       Radial shell fineness. 0->coarse, 1->medium, 2->fine
    etarget       None    Target energy for testing purposes
    targettol     1e-4    Tolerance for energy target sucess
    """
    # Get the various options from the kwargs flags
    verbose = opts.get('verbose',False) 
    etol = opts.get('etol',1e-4)
    max_iter = opts.get('max_iter',20)
    do_averaging = opts.get('do_averaging',True)
    etemp = opts.get('etemp',False)
    functional = opts.get('functional',None)
    etarget = opts.get('etarget',None)
    targettol = opts.get('targettol',1e-4)

    if verbose:
        oldloglevel = logging.root.level
        logging.root.setLevel(logging.DEBUG)

    if functional:
        do_dft = True
        opts['do_grad_dens'] = need_gradients[functional]
        method = functional
    else:
        method = "HF"
        do_dft = False

    bfs = opts.get('bfs',None)
    if not bfs:
        basis_data = opts.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = opts.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    Enuke = atoms.get_enuke()

    if do_dft:
        # default medium mesh
        grid_nrad = opts.get('grid_nrad',32)
        grid_fineness = opts.get('grid_fineness',1)

        gr = MolecularGrid(atoms,grid_nrad,grid_fineness,**opts) 
        gr.set_bf_amps(bfs)

        bfgrid = gr.allbfs() # bfs over all grid points

    # It would be nice to have a more intelligent treatment of the guess
    # so that I could pass in a density rather than a set of orbs.
    orbs = opts.get('orbs',None)
    if orbs is None: orbe,orbs = GHeigenvectors(h,S)
    orbsa = orbsb = orbs
    orbea = orbeb = orbe

    nalpha,nbeta = atoms.get_alphabeta()

    logging.info("U%s calculation on %s" % (method,atoms.name))

    logging.debug("Nbf = %d" % len(bfs))
    logging.debug("Nalpha = %d" % nalpha)
    logging.debug("Nbeta = %d" % nbeta)

    if nalpha == nbeta:
        logging.warning("Warning! nalpha = nbeta: should you be using SCF instead?")

    if etemp:
        logging.info("Electron temperature set to %.0f" % etemp)
        if not functional:
            logging.warning("Warning! Running finite temperature HF")

    eold = 0.
    if do_averaging:
        logging.debug("Using DIIS averaging")
        avga=DIIS(S)
        avgb=DIIS(S)

    logging.debug("Beginning SCF Optimization")
    logging.debug("It   Etot       Eone        Ej        Exc         Enuke")
    logging.debug("--   -------    -------     ------    -------     ------")
    for i in range(max_iter):
        if etemp:
            Da,entropya = mkdens_fermi(2*nalpha,orbea,orbsa,etemp)
            Db,entropyb = mkdens_fermi(2*nbeta,orbeb,orbsb,etemp)
            entropy = (entropya+entropyb)/2.0
        else:
            Da = mkdens(orbsa,0,nalpha)
            Db = mkdens(orbsb,0,nbeta)

        #if do_averaging:
        #    if i:
        #        Da = 0.5*(Da+Da0)
        #        Db = 0.5*(Db+Db0)
        #    Da0 = Da
        #    Db0 = Db

        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)

        Ka = getK(Ints,Da)
        Kb = getK(Ints,Db)
        Fa = h + Ja+Jb-Ka
        Fb = h + Ja+Jb-Kb

        if do_averaging:
            Fa = avga.getF(Fa,Da)
            Fb = avgb.getF(Fb,Db)

        orbea,orbsa = GHeigenvectors(Fa,S)
        orbeb,orbsb = GHeigenvectors(Fb,S)
            
        # Compute energy components
        Dab = Da+Db
        Ej = TraceProperty(Dab,Ja+Jb)/2
        Eone = TraceProperty(Dab,h)
        Exc = -TraceProperty(Da,Ka)/2 - TraceProperty(Db,Kb)/2
        energy = Eone + Ej + Exc + Enuke
        if etemp: energy += entropy
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
                  (i,energy,Eone,Ej,Exc,Enuke))
        if abs(energy-eold) < etol: break
        eold = energy
    logging.info("Final energy for system %s is %f"
                 % (atoms.name,energy))
    if etarget:
        worked = abs(energy-etarget) < targettol
        if worked:
            logging.info("Target energy met")
        else:
            logging.warning("Warning! energy should have been %f" % etarget)

    logging.info("\n")
    if verbose:
        # Set the logging back to what it was before
        logging.root.setLevel(oldloglevel)
    return energy,(orbea,orbeb),(orbsa,orbsb)

def test():
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        #filename='SCF.log',filemode='w',
        )

    h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],units='Angs')
    he = Molecule('He',atomlist = [(2,(0,0,0))])
    li = Molecule('Li',atomlist = [(3,(0,0,0))],multiplicity=2)

    # RHF Tests:
    en,orbe,orbs = SCF(h2,etarget=-1.130501)
    en,orbe,orbs = SCF(he,etarget=-2.855160)
    en,orbe,orbs = SCF(li)

    # RDFT Tests:
    en,orbe,orbs = SCF(h2,functional='SVWN',etarget=-1.132710) 
    en,orbe,orbs = SCF(he,functional='SVWN',etarget=-2.826697)
    en,orbe,orbs = SCF(li,functional='SVWN',etarget=-7.332050)

    # UHF Tests
    # Jaguar gives an energy of -7.431365 for this (essentially what uhf() gives)
    # Energy should be -7.431364. Doesn't quite work...
    en,(orbea,orbeb),(orbsa,orbsb) = USCF(li,etarget=-7.431364)

    # FT Tests
    en,orbe,orbs = SCF(h2,etemp=1e4,etarget=-1.130502)
    en,orbe,orbs = SCF(h2,functional='SVWN',etemp=1e4,etarget=-1.132473)
    en,orbe,orbs = SCF(li,etemp=1e4)
    en,orbe,orbs = SCF(li,functional='SVWN',etemp=1e4,etarget=-7.349422)
    en,(orbea,orbeb),(orbsa,orbsb) = USCF(li,etemp=1e4)
    

if __name__ == '__main__': test()
