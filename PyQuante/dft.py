#!/usr/bin/env python
"""\
 dft.py Initial version of dft routines in PyQuante

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

#from math import *
from Ints import getbasis,getJ,getints
from MolecularGrid import MolecularGrid
from LA2 import GHeigenvectors,mkdens,mkdens_spinavg,TraceProperty
from fermi_dirac import get_efermi, get_fermi_occs,mkdens_occs, get_entropy
from NumWrap import zeros,dot,array,ravel,transpose
from DFunctionals import XC,need_gradients
from time import time
from Convergence import DIIS
from PyQuante.cints import dist
from PyQuante import logging

def getXC(gr,nel,bfgrid,**opts):
    "Form the exchange-correlation matrix"
    # Needs to be rewritten with a more intelligent
    # handling of the functionals
    verbose = opts.get('verbose',False)

    nbf = gr.nbf()
    Vxc = zeros((nbf,nbf),'d')
    dens = gr.dens() # densities at each grid point
    weight = gr.weights()  # weights at each grid point
    gamma = gr.gamma()

    exc,vxc = XC(dens,gamma,**opts)

    renorm_factor = nel/dot(weight,dens)
    weight *= renorm_factor # Renormalize to the proper # electrons

    Exc = dot(weight,exc)
    wv = weight*vxc  # Combine w*v in a vector for multiplication by bfs
    
    for a in range(nbf):
        wva = wv*bfgrid[:,a] 
        for b in range(nbf):
            Vxc[a,b] = dot(wva,bfgrid[:,b])
    return Exc,Vxc

def getXCnew(gr,nel,bfgrid,**opts):
    # going to need to do a
    from DFunctionals import XCNEW

    functional = opts.get('functional','SVWN')
    do_grad_dens = need_gradients[functional]

    dens = gr.dens()
    weight = gr.weights()
    gamma = gr.gamma()
    npts = len(dens)

    # This only happens for non-spin-polarized cases
    amdens = zeros((2,npts),'d')
    amgamma = zeros((3,npts),'d')
    amdens[0,:] = amdens[1,:] = 0.5*dens
    if gamma is not None:
    	amgamma[0,:] = amgamma[1,:] = amgamma[2,:] = 0.25*gamma

    fxc,dfxcdna,dfxcdnb,dfxcdgaa,dfxcdgab,dfxcdgbb = XCNEW(amdens,amgamma,
                                                           **opts)

    renorm_factor = nel/dot(weight,dens)
    weight *= renorm_factor # Renormalize to the proper # electrons

    Exc = dot(weight,fxc)
    wv = weight*dfxcdna  # Combine w*v in a vector for multiplication by bfs

    # Just do the Fxc_a summation here, since we're still looking only
    #  at the unpolarized special case. May have to multiply the result
    #  by two

    # First do the part that doesn't depend upon gamma
    nbf = gr.nbf()
    Fxc = zeros((nbf,nbf),'d')
    for a in range(nbf):
        wva = wv*bfgrid[:,a] 
        for b in range(nbf):
            Fxc[a,b] = dot(wva,bfgrid[:,b])

    # Now do the gamma-dependent part.
    # Fxc_a += dot(2 dfxcdgaa*graddensa + dfxcdgab*graddensb,grad(chia*chib))
    # We can do the dot product as
    #  sum_grid sum_xyz A[grid,xyz]*B[grid,xyz]
    # or a 2d trace product
    # Here A contains the dfxcdgaa stuff
    #      B contains the grad(chia*chib)
    # Yuk. This is ugly...
    if do_grad_dens:
        A = transpose(0.5*transpose(gr.grad())*(weight*(2*dfxcdgaa+dfxcdgab))) # should be (npts,3)
        for a in range(nbf):
            for b in range(nbf):
                B = gr.gradbfab(a,b)
		#test1 = A*B
		#print A.shape,B.shape,test1.shape
                Fxc[a,b] += sum(ravel(A*B))
    return Exc,Fxc # do we need 2*Fxc here?
    

def dft(atoms,**opts):
    """\
    dft(atoms,**opts) - DFT driving routine

    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    verbose       False   Output terse information to stdout (default)
                  True    Print out additional information 
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
    functional    SVWN    Use the SVWN (LDA) DFT functional (default)
                  S0      Use the Slater Xalpha DFT functional
                  BLYP    Use the BLYP GGA DFT functional
                  PBE     Use the PBE DFT functional
    grid_nrad     32      Number of radial shells per atom
    grid_fineness 1       Radial shell fineness. 0->coarse, 1->medium, 2->fine
    spin_type     A       Average occupation method for open shell (default)
                  R       Restricted open shell (not implemented yet)
                  U       Unrestricted open shell (aka spin-polarized dft)
                          Only A works now. Stay tuned.
    """
    verbose = opts.get('verbose',False) 
    ConvCriteria = opts.get('ConvCriteria',1e-4)
    MaxIter = opts.get('MaxIter',20)
    DoAveraging = opts.get('DoAveraging',True)
    ETemp = opts.get('ETemp',False)
    functional = opts.get('functional',None)
    opts['do_grad_dens'] = need_gradients[functional]

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
    enuke = atoms.get_enuke()

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

    logging.debug("DFT calculation on %s using functional %s"
                 % (atoms.name,functional))
    logging.debug("Nbf = %d" % len(bfs))
    logging.debug("Nclosed = %d" % nclosed)
    logging.debug("Nopen = %d" % nclosed)
    if nopen: logging.debug("Using spin-averaged dft for open shell calculation")

    eold = 0.
    if DoAveraging:
        logging.debug("Using DIIS averaging")
        avg=DIIS(S)

    # Converge the LDA density for the system:
    logging.debug("Optimization of DFT density")
    for i in range(MaxIter):
        if ETemp:
            efermi = get_efermi(nel,orbe,ETemp)
            occs = get_fermi_occs(efermi,orbe,ETemp)
            D = mkdens_occs(orbs,occs)
            entropy = get_entropy(occs,ETemp)
        else:
            D = mkdens_spinavg(orbs,nclosed,nopen)
    
        gr.setdens(D)

        J = getJ(Ints,D)

        Exc,XC = getXC(gr,nel,bfgrid,**opts)
            
        F = h+2*J+XC
        if DoAveraging: F = avg.getF(F,D)
        
        orbe,orbs = GHeigenvectors(F,S)
        
        Ej = 2*TraceProperty(D,J)
        Eone = 2*TraceProperty(D,h)
        energy = Eone + Ej + Exc + enuke
        if ETemp: energy += entropy
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
                  (i,energy,Eone,Ej,Exc,enuke))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final %s energy for system %s is %f"
                 % (functional,atoms.name,energy))
    return energy,orbe,orbs

def dftnew(atoms,**opts):
    """\
    dft(atoms,**opts) - DFT driving routine

    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    verbose       False   Output terse information to stdout (default)
                  True    Print out additional information 
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
    functional    SVWN    Use the SVWN (LDA) DFT functional (default)
                  S0      Use the Slater Xalpha DFT functional
                  BLYP    Use the BLYP GGA DFT functional
                  PBE     Use the PBE DFT functional
    grid_nrad     32      Number of radial shells per atom
    grid_fineness 1       Radial shell fineness. 0->coarse, 1->medium, 2->fine
    spin_type     A       Average occupation method for open shell (default)
                  R       Restricted open shell (not implemented yet)
                  U       Unrestricted open shell (aka spin-polarized dft)
                          Only A works now. Stay tuned.
    """
    verbose = opts.get('verbose',False) 
    ConvCriteria = opts.get('ConvCriteria',1e-4)
    MaxIter = opts.get('MaxIter',20)
    DoAveraging = opts.get('DoAveraging',True)
    ETemp = opts.get('ETemp',False)
    functional = opts.get('functional','SVWN')
    opts['do_grad_dens'] = need_gradients[functional]

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
    enuke = atoms.get_enuke()

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

    logging.debug("DFT calculation on %s using functional %s"
                 % (atoms.name,functional))
    logging.debug("Nbf = %d" % len(bfs))
    logging.debug("Nclosed = %d" % nclosed)
    logging.debug("Nopen = %d" % nclosed)
    if nopen: logging.debug("Using spin-averaged dft for open shell calculation")

    eold = 0.
    if DoAveraging:
        logging.debug("Using DIIS averaging")
        avg=DIIS(S)

    # Converge the LDA density for the system:
    logging.debug("Optimization of DFT density")
    for i in range(MaxIter):
        if ETemp:
            efermi = get_efermi(nel,orbe,ETemp)
            occs = get_fermi_occs(efermi,orbe,ETemp)
            D = mkdens_occs(orbs,occs)
            entropy = get_entropy(occs,ETemp)
        else:
            D = mkdens_spinavg(orbs,nclosed,nopen)
    
        gr.setdens(D)

        J = getJ(Ints,D)

        Exc,XC = getXCnew(gr,nel,bfgrid,**opts)
            
        F = h+2*J+XC
        if DoAveraging: F = avg.getF(F,D)
        
        orbe,orbs = GHeigenvectors(F,S)
        
        Ej = 2*TraceProperty(D,J)
        Eone = 2*TraceProperty(D,h)
        energy = Eone + Ej + Exc + enuke
        if ETemp: energy += entropy
        logging.debug("%d %14.8f %14.8f %14.8f %14.8f %14.8f" %
                  (i,energy,Eone,Ej,Exc,enuke))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final %s energy for system %s is %f"
                 % (functional,atoms.name,energy))
    return energy,orbe,orbs

def udft(atoms,**opts):
    """\
    udft(atoms,**opts) - Unrestricted spin DFT driving routine

    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    verbose       False   Output terse information to stdout (default)
                  True    Print out additional information 
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
    functional    SVWN    Use the SVWN (LDA) DFT functional (default)
                  S0      Use the Slater Xalpha DFT functional
                  BLYP    Use the BLYP GGA DFT functional
                  PBE     Use the PBE DFT functional
    grid_nrad     32      Number of radial shells per atom
    grid_fineness 1       Radial shell fineness. 0->coarse, 1->medium, 2->fine
    """
    verbose = opts.get('verbose',False) 
    ConvCriteria = opts.get('ConvCriteria',1e-4)
    MaxIter = opts.get('MaxIter',20)
    DoAveraging = opts.get('DoAveraging',True)
    ETemp = opts.get('ETemp',False)
    functional = opts.get('functional','SVWN')
    opts['do_grad_dens'] = need_gradients[functional]

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
    enuke = atoms.get_enuke()

    # default medium mesh
    grid_nrad = opts.get('grid_nrad',32)
    grid_fineness = opts.get('grid_fineness',1)

    gr = MolecularGrid(atoms,grid_nrad,grid_fineness,**opts) 
    gr.set_bf_amps(bfs)

    bfgrid = gr.allbfs() # bfs over all grid points

    # It would be nice to have a more intelligent treatment of the guess
    # so that I could pass in a density rather than a set of orbs.
    orbs = opts.get('orbs',None)
    if not orbs: orbe,orbs = GHeigenvectors(h,S)
    orbsa = orbsb = orbs

    nalpha,nbeta = atoms.get_alphabeta()

    logging.debug("UDFT calculation on %s using functional %s"
                 % (atoms.name,functional))
    logging.debug("Nbf = %d" % len(bfs))
    logging.debug("Nalpha = %d" % nalpha)
    logging.debug("Nbeta = %d" % nbeta)
        
    eold = 0.

    # Converge the LDA density for the system:
    logging.debug("Optimization of DFT density")
    for i in range(MaxIter):
        Da = mkdens(orbsa,0,nalpha)
        Db = mkdens(orbsb,0,nbeta)

        gr.setdens_sp(Da,Db)

        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)

        Exc,XCa,XCb = getXC_sp(gr,nel,bfgrid,**opts)
            
        Fa = h+Ja+Jb-Ka
        Fb = h+Ja+Jb-Kb
        
        orbea,orbsa = GHeigenvectors(Fa,S)
        orbeb,orbsb = GHeigenvectors(Fb,S)
        
        Eja = TraceProperty(D,Ja)
        Ejb = TraceProperty(D,Jb)
        Eone = 2*TraceProperty(D,h)
        energy = Eone + Eja + Ejb + Exc + enuke
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
                  (i,energy,Eone,Eja+Ejb,Exc,enuke))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final U%s energy for system %s is %f"
                 % (functional,atoms.name,energy))
    return energy,orbe,orbs


if __name__ == '__main__':
    from Molecule import Molecule
    h2 = Molecule('h2o',
                  [( 1,(0,0,0.35)),
                   ( 1,(0,0,-0.35))],
                   units='Angs')
    en = dft(h2)
    print en

