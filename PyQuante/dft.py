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
#from MolecularGrid import MolecularGrid
from MG2 import MG2 as MolecularGrid
from LA2 import geigh,mkdens,mkdens_spinavg,trace2
from fermi_dirac import get_efermi, get_fermi_occs,mkdens_occs, get_entropy
from NumWrap import zeros,dot,ravel,transpose,sum
from DFunctionals import XC,need_gradients
from time import time
from Convergence import DIIS
from PyQuante.cints import dist
import logging

from Convergence import SimpleAverager
from LA2 import SymOrth

def getXC(gr,nel,**kwargs):
    "Form the exchange-correlation matrix"

    functional = kwargs.get('functional','SVWN')
    do_grad_dens = need_gradients[functional]
    do_spin_polarized = kwargs.get('do_spin_polarized')
    
    gr.floor_density()  # Insure that the values of the density don't underflow
    gr.renormalize(nel) # Renormalize to the proper # electrons

    dens = gr.dens()
    weight = gr.weights()
    gamma = gr.get_gamma()
    npts = len(dens)

    if gr.version == 1:
        amdens = zeros((2,npts),'d')
        amgamma = zeros((3,npts),'d')
        amdens[0,:] = amdens[1,:] = 0.5*dens
        if gamma is not None:
            amgamma[0,:] = amgamma[1,:] = amgamma[2,:] = 0.25*gamma
    elif gr.version == 2:
        amdens = gr.density.T
        amgamma = gr.gamma.T

    fxc,dfxcdna,dfxcdnb,dfxcdgaa,dfxcdgab,dfxcdgbb = XC(amdens,amgamma,**kwargs)

    Exc = dot(weight,fxc)

    wva = weight*dfxcdna  # Combine w*v in a vector for multiplication by bfs

    # First do the part that doesn't depend upon gamma
    nbf = gr.get_nbf()
    Fxca = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        wva_i = wva*gr.bfgrid[:,i] 
        for j in xrange(nbf):
            Fxca[i,j] = dot(wva_i,gr.bfgrid[:,j])

    # Now do the gamma-dependent part.
    # Fxc_a += dot(2 dfxcdgaa*graddensa + dfxcdgab*graddensb,grad(chia*chib))
    # We can do the dot product as
    #  sum_grid sum_xyz A[grid,xyz]*B[grid,xyz]
    # or a 2d trace product
    # Here A contains the dfxcdgaa stuff
    #      B contains the grad(chia*chib)

    # Possible errors: gr.grad() here should be the grad of the b part?
    if do_grad_dens:
        # A,B are dimensioned (npts,3)
        A = transpose(0.5*transpose(gr.grad())*(weight*(2*dfxcdgaa+dfxcdgab)))
        for a in xrange(nbf):
            for b in xrange(a+1):
                B = gr.grad_bf_prod(a,b)
                Fxca[a,b] += sum(ravel(A*B))
                Fxca[b,a] = Fxca[a,b]
    if not do_spin_polarized: return Exc,Fxca

    wvb = weight*dfxcdnb  # Combine w*v in a vector for multiplication by bfs

    # First do the part that doesn't depend upon gamma
    Fxcb = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        wvb_i = wvb*gr.bfgrid[:,i] 
        for j in xrange(nbf):
            Fxcb[i,j] = dot(wvb_i,gr.bfgrid[:,j])

    # Now do the gamma-dependent part.
    # Fxc_b += dot(2 dfxcdgbb*graddensb + dfxcdgab*graddensa,grad(chia*chib))
    # We can do the dot product as
    #  sum_grid sum_xyz A[grid,xyz]*B[grid,xyz]
    # or a 2d trace product
    # Here A contains the dfxcdgaa stuff
    #      B contains the grad(chia*chib)

    # Possible errors: gr.grad() here should be the grad of the b part?
    if do_grad_dens:
        # A,B are dimensioned (npts,3)
        A = transpose(0.5*transpose(gr.grad())*(weight*(2*dfxcdgbb+dfxcdgab)))
        for a in xrange(nbf):
            for b in xrange(a+1):
                B = gr.grad_bf_prod(a,b)
                Fxcb[a,b] += sum(ravel(A*B))
                Fxcb[b,a] = Fxcb[a,b]
    return Exc,Fxca,Fxcb

def dft(atoms,**kwargs):
    """\
    dft(atoms,**kwargs) - DFT driving routine

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
    verbose = kwargs.get('verbose',False) 
    ConvCriteria = kwargs.get('ConvCriteria',1e-4)
    MaxIter = kwargs.get('MaxIter',20)
    DoAveraging = kwargs.get('DoAveraging',True)
    ETemp = kwargs.get('ETemp',False)
    functional = kwargs.get('functional','SVWN')
    kwargs['do_grad_dens'] = need_gradients[functional]

    bfs = kwargs.get('bfs',None)
    if not bfs:
        basis_data = kwargs.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = kwargs.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    enuke = atoms.get_enuke()

    # default medium mesh
    grid_nrad = kwargs.get('grid_nrad',32)
    grid_fineness = kwargs.get('grid_fineness',1)

    gr = MolecularGrid(atoms,grid_nrad,grid_fineness,**kwargs) 
    gr.set_bf_amps(bfs)

    # It would be nice to have a more intelligent treatment of the guess
    # so that I could pass in a density rather than a set of orbs.
    orbs = kwargs.get('orbs',None)
    if orbs is None: orbe,orbs = geigh(h,S)

    nclosed,nopen = atoms.get_closedopen()

    if verbose:
        print "DFT calculation on %s using functional %s" % (
            atoms.name,functional)
        print "Nbf = %d" % len(bfs)
        print "Nclosed = %d" % nclosed
        print "Nopen = %d" % nclosed
        if nopen: print "Using spin-averaged dft for open shell calculation"

    eold = 0.
    if DoAveraging:
        if verbose: print"Using DIIS averaging"
        avg=DIIS(S)

    # Converge the LDA density for the system:
    if verbose: print "Optimization of DFT density"
    for i in xrange(MaxIter):
        if ETemp:
            efermi = get_efermi(nel,orbe,ETemp)
            occs = get_fermi_occs(efermi,orbe,ETemp)
            D = mkdens_occs(orbs,occs)
            entropy = get_entropy(occs,ETemp)
        else:
            D = mkdens_spinavg(orbs,nclosed,nopen)
    
        gr.setdens(D)

        J = getJ(Ints,D)

        Exc,XC = getXC(gr,nel,**kwargs)
            
        F = h+2*J+XC
        if DoAveraging: F = avg.getF(F,D)
        
        orbe,orbs = geigh(F,S)
        
        Ej = 2*trace2(D,J)
        Eone = 2*trace2(D,h)
        energy = Eone + Ej + Exc + enuke
        if ETemp: energy += entropy
        if verbose: print "%d %10.4f %10.4f %10.4f %10.4f %10.4f" % \
           (i,energy,Eone,Ej,Exc,enuke)
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    print "Final %s energy for system %s is %f" % (functional,atoms.name,energy)
    return energy,orbe,orbs

def udft(atoms,**kwargs):
    """\
    udft(atoms,**kwargs) - Unrestricted spin DFT driving routine

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
    verbose = kwargs.get('verbose',False) 
    ConvCriteria = kwargs.get('ConvCriteria',1e-4)
    MaxIter = kwargs.get('MaxIter',20)
    DoAveraging = kwargs.get('DoAveraging',True)
    ETemp = kwargs.get('ETemp',False)
    functional = kwargs.get('functional','SVWN')
    kwargs['do_grad_dens'] = need_gradients[functional]
    kwargs['do_spin_polarized'] = True

    bfs = kwargs.get('bfs',None)
    if not bfs:
        basis_data = kwargs.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = kwargs.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    enuke = atoms.get_enuke()

    # default medium mesh
    grid_nrad = kwargs.get('grid_nrad',32)
    grid_fineness = kwargs.get('grid_fineness',1)

    gr = MolecularGrid(atoms,grid_nrad,grid_fineness,**kwargs) 
    gr.set_bf_amps(bfs)

    # It would be nice to have a more intelligent treatment of the guess
    # so that I could pass in a density rather than a set of orbs.
    orbs = kwargs.get('orbs',None)
    if not orbs: orbe,orbs = geigh(h,S)
    orbsa = orbsb = orbs

    nalpha,nbeta = atoms.get_alphabeta()

    if verbose: 
        print "UDFT calculation on %s using functional %s" \
              % (atoms.name,functional)
        print "Nbf = %d" % len(bfs)
        print "Nalpha = %d" % nalpha
        print "Nbeta = %d" % nbeta
        
    eold = 0.

    # Converge the LDA density for the system:
    if verbose: print "Optimization of DFT density"
    for i in xrange(MaxIter):
        Da = mkdens(orbsa,0,nalpha)
        Db = mkdens(orbsb,0,nbeta)

        gr.setdens(Da,Db)

        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)

        Exc,XCa,XCb = getXC(gr,nel,**kwargs)
            
        Fa = h+Ja+Jb-Ka
        Fb = h+Ja+Jb-Kb
        
        orbea,orbsa = geigh(Fa,S)
        orbeb,orbsb = geigh(Fb,S)
        
        Eja = trace2(D,Ja)
        Ejb = trace2(D,Jb)
        Eone = 2*trace2(D,h)
        energy = Eone + Eja + Ejb + Exc + enuke
        if verbose:
            print "%d %10.4f %10.4f %10.4f %10.4f %10.4f" % (
                i,energy,Eone,Eja+Ejb,Exc,enuke)
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    print "Final U%s energy for system %s is %f" % (
        functional,atoms.name,energy)
    return energy,orbe,orbs


def mk_auger_dens(c, occ):
    "Forms a density matrix from a coef matrix c and occupations in occ"
    #count how many states we were given
    nstates = occ.shape[0]
    D = 0.0
    for i in xrange(nstates):
        D += occ[i]*dot( c[:,i:i+1], transpose(c[:,i:i+1]))
    #pad_out(D)
    return D
    
    
def dft_fixed_occ(atoms,occs,**kwargs):
    """\
    dft(atoms,**kwargs) - DFT driving routine

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
    verbose = kwargs.get('verbose',False) 
    ConvCriteria = kwargs.get('ConvCriteria',1e-4)
    MaxIter = kwargs.get('MaxIter',20)
    DoAveraging = kwargs.get('DoAveraging',True)
    ETemp = kwargs.get('ETemp',False)
    functional = kwargs.get('functional','SVWN')
    kwargs['do_grad_dens'] = need_gradients[functional]

    bfs = kwargs.get('bfs',None)
    if not bfs:
        basis_data = kwargs.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = kwargs.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    enuke = atoms.get_enuke()

    # default medium mesh
    grid_nrad = kwargs.get('grid_nrad',32)
    grid_fineness = kwargs.get('grid_fineness',1)

    gr = MolecularGrid(atoms,grid_nrad,grid_fineness,**kwargs) 
    gr.set_bf_amps(bfs)

    # It would be nice to have a more intelligent treatment of the guess
    # so that I could pass in a density rather than a set of orbs.
    orbs = kwargs.get('orbs',None)
    if orbs is None: orbe,orbs = geigh(h,S)

    nclosed,nopen = atoms.get_closedopen()

    if verbose: 
        print "DFT calculation on %s using functional %s" % (
            atoms.name,functional)
        print "Nbf = %d" % len(bfs)
        print "Nclosed = %d" % nclosed
        print "Nopen = %d" % nclosed
        if nopen: print "Using spin-averaged dft for open shell calculation"

    eold = 0.
    if DoAveraging:
        print "Using DIIS averaging"
        avg=DIIS(S)

    # Converge the LDA density for the system:
    if verbose: print "Optimization of DFT density"
    for i in xrange(MaxIter):
        #print "SCF Iteration:",i,"Starting Energy:",eold
        #save the starting orbitals
        oldorbs=orbs
        if ETemp:
            efermi = get_efermi(nel,orbe,ETemp)
            occs = get_fermi_occs(efermi,orbe,ETemp)
            D = mkdens_occs(orbs,occs)
            entropy = get_entropy(occs,ETemp)
        else:
            D = mk_auger_dens(orbs, occs)
            #D = mkdens_spinavg(orbs,nclosed,nopen)
    
        gr.setdens(D)

        J = getJ(Ints,D)

        Exc,XC = getXC(gr,nel,**kwargs)
            
        F = h+2*J+XC
        if DoAveraging: F = avg.getF(F,D)
        
        orbe,orbs = geigh(F,S)
        #pad_out(orbs)
        #save the new eigenstates of the fock operator F
        neworbs=orbs
        #send oldorbs and neworbs to get biorthogonalized method
        #orbs = biorthog(neworbs, oldorbs, S, nel)
        
        Ej = 2*trace2(D,J)
        Eone = 2*trace2(D,h)
        energy = Eone + Ej + Exc + enuke
        print i+1,"   ",energy,"    ",Ej,"  ",Exc," ",enuke
        if ETemp: energy += entropy
        if verbose:
            print "%d %10.4f %10.4f %10.4f %10.4f %10.4f" % (
                i,energy,Eone,Ej,Exc,enuke)
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
        
    print "Final %s energy for system %s is %f" % (functional,atoms.name,energy)
    return energy,orbe,orbs
    
def udft_fixed_occ(atoms,occa, occb, **kwargs):
    """\
    occa and occb represent the orbital occupation arrays for 
    the calculating spin orbitals with holes
    
    udft(atoms,**kwargs) - Unrestricted spin DFT driving routine

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
    verbose = kwargs.get('verbose',False) 
    ConvCriteria = kwargs.get('ConvCriteria',1e-4)
    MaxIter = kwargs.get('MaxIter',20)
    DoAveraging = kwargs.get('DoAveraging',True)
    averaging = kwargs.get('averaging',0.5)
    ETemp = kwargs.get('ETemp',False)
    functional = kwargs.get('functional','LDA') 
    #default to LDA which has no correlation since that is easier

    kwargs['do_grad_dens'] = need_gradients[functional]
    kwargs['do_spin_polarized'] = True

    bfs = kwargs.get('bfs',None)
    if not bfs:
        basis_data = kwargs.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = kwargs.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    enuke = atoms.get_enuke()

    # default medium mesh
    grid_nrad = kwargs.get('grid_nrad',32)
    grid_fineness = kwargs.get('grid_fineness',1)

    gr = MolecularGrid(atoms,grid_nrad,grid_fineness,**kwargs) 
    gr.set_bf_amps(bfs)

    # It would be nice to have a more intelligent treatment of the guess
    # so that I could pass in a density rather than a set of orbs.
    orbs = kwargs.get('orbs',None)
    if not orbs: orbe,orbs = geigh(h,S)
    orbsa = orbsb = orbs

    nalpha,nbeta = atoms.get_alphabeta()

    if verbose:
        print "UDFT calculation on %s using functional %s" % (
            atoms.name,functional)
        print "Nbf = %d" % len(bfs)
        print "Nalpha = %d" % nalpha
        print "Nbeta = %d" % nbeta
        
    eold = 0.

    # Converge the LDA density for the system:
    print "Optimization of DFT density"
    for i in xrange(MaxIter):
        Da = mk_auger_dens(orbsa, occa)
        Db = mk_auger_dens(orbsb, occb)
        Dab = Da + Db
        
        if DoAveraging:
            if i: 
                Da = averaging*Da + (1-averaging)*Da0
                Db = averaging*Db + (1-averaging)*Db0
            Da0 = Da
            Db0 = Db

        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)

        #remember we must use a functional that has no correlation energy
        gr.setdens(Da,Db)
        exca,XCa,XCb = getXC(gr,nel,**kwargs)
        
        Fa = h+Ja+Jb+XCa
        Fb = h+Ja+Jb+XCb
        
        orbea,orbsa = geigh(Fa,S)
        orbeb,orbsb = geigh(Fb,S)
        
        Eone = trace2(Dab,h)
        Ej   = 0.5*trace2(Dab,Ja+Jb)
        Exc = 0.5*exca + 0.5*excb
        
        energy = Eone + Ej + Exc + enuke
        
        print i+1,"   ",energy,"    ",Ej,"  ",Exc
        if verbose: 
            print "%d %10.4f %10.4f %10.4f" % (i,energy,Ej,Exc)
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    print "Final U%s energy for system %s is %f" % (
        functional,atoms.name,energy)
    return energy,(orbea,orbeb),(orbsa,orbsb)

if __name__ == '__main__':
    from Molecule import Molecule
    h2 = Molecule('h2o',
                  [( 1,(0,0,0.35)),
                   ( 1,(0,0,-0.35))],
                   units='Angs')
    en = dft(h2)
    print en

    
