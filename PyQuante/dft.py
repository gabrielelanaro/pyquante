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
from LA2 import geigh,mkdens,mkdens_spinavg,trace2
from fermi_dirac import get_efermi, get_fermi_occs,mkdens_occs, get_entropy
from NumWrap import zeros,dot,array,ravel,transpose
from DFunctionals import XC,need_gradients
from time import time
from Convergence import DIIS
from PyQuante.cints import dist
import logging

from Convergence import SimpleAverager
from LA2 import SymOrth
from numpy import resize
from numpy.oldnumeric.linear_algebra import Heigenvectors, inverse
from numpy.linalg import svd, eig, eigh
from NumWrap import array2string, abs
tolerance = 0.001

# This is the version before Ann Mattsson made her changes. I'm keeping
#  it around for old time's sake. Putting this comment in 2007-02;
#  should delete the function after 2007-08.
def getXCold(gr,nel,bfgrid,**opts):
    "Form the exchange-correlation matrix"
    # Needs to be rewritten with a more intelligent
    # handling of the functionals
    from DFunctionals import XCold
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

def getXC(gr,nel,bfgrid,**opts):
    "Form the exchange-correlation matrix"

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

    fxc,dfxcdna,dfxcdnb,dfxcdgaa,dfxcdgab,dfxcdgbb = XC(amdens,amgamma,
                                                        **opts)

    # Renormalize to the proper # electrons
    renorm_factor = nel/dot(weight,dens)
    weight *= renorm_factor 

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
        # A is dimensioned (npts,3)
        A = transpose(0.5*transpose(gr.grad())*(weight*(2*dfxcdgaa+dfxcdgab)))
        B = zeros(A.shape,'d')
        for a in range(nbf):
            agrid = bfgrid[:,a]
            agradgrid = gr.bfgrad(a)
            for b in range(a+1):
                bgrid = bfgrid[:,b]
                bgradgrid = gr.bfgrad(b)
                for i in range(3):
                    B[:,i] = agrid[:]*bgradgrid[:,i] + bgrid[:]*agradgrid[:,i]
                # This is the original method, which was painfully slow:
                #B = gr.gradbfab(a,b)
                Fxc[a,b] += sum(ravel(A*B))
                Fxc[b,a] = Fxc[a,b]
    return Exc,Fxc
    

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
    if orbs is None: orbe,orbs = geigh(h,S)

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
        
        orbe,orbs = geigh(F,S)
        
        Ej = 2*trace2(D,J)
        Eone = 2*trace2(D,h)
        energy = Eone + Ej + Exc + enuke
        if ETemp: energy += entropy
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
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
    if not orbs: orbe,orbs = geigh(h,S)
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
        
        orbea,orbsa = geigh(Fa,S)
        orbeb,orbsb = geigh(Fb,S)
        
        Eja = trace2(D,Ja)
        Ejb = trace2(D,Jb)
        Eone = 2*trace2(D,h)
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


#modified by Hatem Helal hhh23@cam.ac.uk
def pad_out(matrix):
    #this will make debugging matrix operations easier by getting rid of 
    #matrix elements which are ridiculously tiny
    rows = matrix.shape[0]
    cols = matrix.shape[1]
    
    for i in range(rows):
        for j in range(cols):
                if abs(matrix[i][j]) < tolerance:
                    matrix[i][j]=0.0
    
    print array2string(matrix,200,precision=4);    print "\n\n"

    return 0
    
def diagonalize(vector):
    #takes a vector with N components and returns an NXN matrix whose diagonal elements
    #are the components of the vector, used for SVD below
    len = vector.shape[0]
    matrix = array(0.0)
    matrix = resize(matrix,(len,len))
    
    i=0
    j=0
    for i in range(len):
        for j in range(len):
            if i==j:
                matrix[i][j]= vector[i]
            else:
                matrix[i][j] = 0.0
    return matrix

def Ortho_Orbs(orbs,S):
    #Given an orbital coef matrix 'orbs' and the basis set overlap matrix
    #S this function will return the orbital coefs in the symmetrically orthogonalized
    #basis
    X = SymOrth(S)
    #in the transformed basis the basis functions are orthonormal, the following tests this
    #pad_out( dot(transpose(X),dot(sigma, X)) )
    
    #we now transform the orbital coef matrices C to C' = inverse(X) C
    #this puts the orbitals in the orthonormal basis defined by the X transformation
    orth_orbs = dot( inverse(X), orbs )

    #test that the orbitals are orthonormal...
    #pad_out( dot(transpose(orth_orbs), orth_orbs) )
    
    return orth_orbs

def biorthog(neworbs, oldorbs, sigma, nelec):
    """
    this function runs some tests on the orbital coefs and makes two sets of orbitals
    biorthogonal...
    """
    #first lets have a look at the orbital coeff matrices we were passed
    #pad_out(oldorbs)
    #pad_out(dot(inverse(SymOrth(sigma)) , oldorbs) )
    #holes = oldorbs.shape[1] - nelec
    #print "nelec,holes",nelec, holes
    oldorbs = oldorbs[...,0:nelec]

    #print "Occupied Old orbitals and all the new orbitals"
    #pad_out(oldorbs)
    #pad_out(neworbs)
    
    #pad_out(sigma)  
    #sigma is the basis set overlap matrix we now want to find the transformation X
    #that gives X^t sigma X = identity ie an orthogonalization transformation
    #pad_out(sigma)
    #calculate orthogonality transformation X using symmetric orthogonalization
    #X = sigma^-1/2 => sigma^-1/2 sigma sigma^-1/2 = sigma^-1/2 sigma^1/2 = identity
    #see Szabo and Ostlund pg 143 for more
    X = SymOrth(sigma)
    #in the transformed basis the basis functions are orthonormal, the following tests this
    #pad_out( dot(transpose(X),dot(sigma, X)) )
    
    #we now transform the orbital coef matrices C to C' = inverse(X) C
    #this puts the orbitals in the orthonormal basis defined by the X transformation
    orth_oldorbs = dot( inverse(X), oldorbs )
    orth_neworbs = dot( inverse(X), neworbs )
    
    #print "Orbitals in orthogonalized basis, old first then new"
    pad_out(orth_oldorbs)
    pad_out(orth_neworbs)

    #now we test if the two sets of orbitals are internally orthogonal
    #this is done by computing transpose(coef_matrix)*coef_matrix
    #since we are now working in an orthogonal basis
    #print "\n\nTesting internal orthogonality of our orbital sets\n"
    #print "Orbital overlap matrix Sij = < i old | j old>"
    #pad_out( dot( transpose(orth_oldorbs), orth_oldorbs ) )  
    #print "Orbital overlap matrix Sij = < i new | j new>"
    #pad_out( dot( transpose(orth_neworbs), orth_neworbs ) ) 
    
    #now we compute the overlap matrix between the new orbitals and the old ones
    S = dot( transpose(orth_oldorbs), orth_neworbs ) 
    Sdag = transpose(S)
    
    print "S" ; pad_out(S)
    #print "Sdag" ; pad_out(Sdag)
    
    print "Sdag S" ; pad_out(dot(Sdag, S))
    print "S Sdag" ; pad_out(dot(S, Sdag))
    """    
    #now we perform singular value decomposition to diagonlize the inter-scf
    #orbital overlap matrix S

    [u,sd,v] = svd(S,compute_uv=1,full_matrices=1)
    sd = diagonalize(sd)
    v=transpose(v)
    print "u" ; pad_out(u)
    print "sd"; pad_out(sd)
    print "v" ; pad_out(v)
    
    udag = transpose(u)
    vdag = transpose(v)
    
    #print "u sd vdag = S" ; pad_out(dot(u,dot(sd,vdag)))
    
    #print "udag S v  = sd" ; pad_out(dot(udag,dot(S,v)))

    #print "Cortho" ; pad_out( orth_oldorbs )
    #print "udag Cortho v" ; pad_out( dot( udag, dot(orth_oldorbs, v) ) )
    #print "Dortho" ; pad_out( orth_neworbs )
    #print "vdag Dortho u" ; pad_out( dot( vdag, dot(orth_neworbs, u) ) )

    
    #trans_oldorbs = dot( udag, dot(orth_oldorbs, v) )
    trans_neworbs = dot( orth_neworbs, v)
    pad_out(trans_neworbs)
    """

    val,v = eig(dot(Sdag, S))
    print val
    print v    
    
    pad_out(dot(X, dot(orth_neworbs, v)))
    pad_out(oldorbs)
    
    pad_out(dot(S,v))

    final_orbs = dot(X, dot(orth_neworbs, v))
    return final_orbs
    
def mk_auger_dens(c, occ):
    "Forms a density matrix from a coef matrix c and occupations in occ"
    #count how many states we were given
    nstates = occ.shape[0]
    D = 0.0
    for i in range(nstates):
        D += occ[i]*dot( c[:,i:i+1], transpose(c[:,i:i+1]))
    #pad_out(D)
    return D
    
def dft_hh(atoms,**opts):
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
    if orbs is None: orbe,orbs = geigh(h,S)

    #orthogonalize orbitals
    orbs = Ortho_Orbs(orbs, S)

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
        print "SCF Iteration:",i,"Starting Energy:",eold
        #save the starting orbitals
        oldorbs=orbs
        
        #orthogonalize orbitals   
        #orbs = Ortho_Orbs(orbs, S)
        
        if ETemp:
            efermi = get_efermi(nel,orbe,ETemp)
            occs = get_fermi_occs(efermi,orbe,ETemp)
            D = mkdens_occs(orbs,occs)
            entropy = get_entropy(occs,ETemp)
        else:
            #D = mk_auger_dens(orbs, occ)
            D = mkdens_spinavg(orbs,nclosed,nopen)
            #pad_out(orbs)
            #pad_out(D)
            #pad_out(S)
            #pad_out(0.5*dot(S,dot(D,S)))
    
        gr.setdens(D)

        J = getJ(Ints,D)

        Exc,XC = getXC(gr,nel,bfgrid,**opts)
            
        F = h+2*J+XC
        if DoAveraging: F = avg.getF(F,D)
        
        orbe,orbs = geigh(F,S)
       
        #save the new eigenstates of the fock operator F
        neworbs=orbs
        #send oldorbs and neworbs to get biorthogonalized method
        #orbs = biorthog(neworbs, oldorbs, S, nel)
            
        Ej = 2*trace2(D,J)
        Eone = 2*trace2(D,h)
        energy = Eone + Ej + Exc + enuke
        print "Ending Energy:",energy
        if ETemp: energy += entropy
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
                  (i,energy,Eone,Ej,Exc,enuke))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final %s energy for system %s is %f"
                 % (functional,atoms.name,energy))
    return energy,orbe,orbs
    
def dft_lin_mix(atoms,**opts):
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
    DoAveraging = opts.get('DoAveraging',False)
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
    if orbs is None: orbe,orbs = geigh(h,S)

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

    #added by Hatem Helal hhh23@cam.ac.uk
    logging.debug("Using linear mixing of density matrices")
    avg=SimpleAverager(0.2)

    # Converge the LDA density for the system:
    logging.debug("Optimization of DFT density")
    for i in range(MaxIter):
        print "SCF Iteration:",i,"Starting Energy:",eold
        #save the starting orbitals
        oldorbs=orbs
        if ETemp:
            efermi = get_efermi(nel,orbe,ETemp)
            occs = get_fermi_occs(efermi,orbe,ETemp)
            D = mkdens_occs(orbs,occs)
            entropy = get_entropy(occs,ETemp)
        else:
            D = mkdens_spinavg(orbs,nclosed,nopen)
    
        gr.setdens(D)
        D = avg.getD(D)
        J = getJ(Ints,D)

        Exc,XC = getXC(gr,nel,bfgrid,**opts)
            
        F = h+2*J+XC
        #if DoAveraging: F = avg.getF(F,D)
        
        orbe,orbs = geigh(F,S)
        pad_out(orbs)
        #save the new eigenstates of the fock operator F
        neworbs=orbs
        #send oldorbs and neworbs to get biorthogonalized method
        #orbs = biorthog(neworbs, oldorbs, S, nel)
        
        Ej = 2*trace2(D,J)
        Eone = 2*trace2(D,h)
        energy = Eone + Ej + Exc + enuke
        print "Ending Energy:",energy
        if ETemp: energy += entropy
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
                  (i,energy,Eone,Ej,Exc,enuke))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final %s energy for system %s is %f"
                 % (functional,atoms.name,energy))
    return energy,orbe,orbs
    
    
def dft_fixed_occ(atoms,occs,**opts):
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
    if orbs is None: orbe,orbs = geigh(h,S)

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

    #added by Hatem Helal hhh23@cam.ac.uk
    #logging.debug("Using linear mixing of density matrices")
    #avg=SimpleAverager(0.5)

    print "SCF Iteration    Total Energy    Coulomb Exchange/Correlation    Nuclear"
    # Converge the LDA density for the system:
    logging.debug("Optimization of DFT density")
    for i in range(MaxIter):
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

        Exc,XC = getXC(gr,nel,bfgrid,**opts)
            
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
        logging.debug("%d %10.4f %10.4f %10.4f %10.4f %10.4f" %
                  (i,energy,Eone,Ej,Exc,enuke))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
        
    logging.info("Final %s energy for system %s is %f"
                 % (functional,atoms.name,energy))
    return energy,orbe,orbs
    
def udft_fixed_occ(atoms,occa, occb, **opts):
    """\
    occa and occb represent the orbital occupation arrays for 
    the calculating spin orbitals with holes
    
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
    averaging = opts.get('averaging',0.5)
    ETemp = opts.get('ETemp',False)
    functional = opts.get('functional','LDA') 
    #default to LDA which has no correlation since that is easier
    
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
    if not orbs: orbe,orbs = geigh(h,S)
    orbsa = orbsb = orbs

    nalpha,nbeta = atoms.get_alphabeta()

    logging.debug("UDFT calculation on %s using functional %s"
                 % (atoms.name,functional))
    logging.debug("Nbf = %d" % len(bfs))
    logging.debug("Nalpha = %d" % nalpha)
    logging.debug("Nbeta = %d" % nbeta)
        
    eold = 0.

    print "SCF Iteration    Total Energy    Coulomb Exchange/Correlation"

    # Converge the LDA density for the system:
    logging.debug("Optimization of DFT density")
    for i in range(MaxIter):
        #Da_std = mkdens(orbsa,0,nalpha)
        #Db_std = mkdens(orbsb,0,nbeta)
        
        Da = mk_auger_dens(orbsa, occa)
        Db = mk_auger_dens(orbsb, occb)
        Dab = Da + Db
        #these should have 0.0 for all elements if everything worked alright 
        #pad_out(Da - Da_std)
        #pad_out(Db - Db_std)
        
        if DoAveraging:
            if i: 
                Da = averaging*Da + (1-averaging)*Da0
                Db = averaging*Db + (1-averaging)*Db0
            Da0 = Da
            Db0 = Db

        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)

        #remember we must use a functional that has no correlation energy
        gr.setdens(Da)
        exca,XCa = getXC(gr,nel,bfgrid,**opts)
        
        gr.setdens(Db)
        excb,XCb = getXC(gr,nel,bfgrid,**opts)
        
        Fa = h+Ja+Jb+XCa
        Fb = h+Ja+Jb+XCb
        
        orbea,orbsa = geigh(Fa,S)
        orbeb,orbsb = geigh(Fb,S)
        
        Eone = trace2(Dab,h)
        Ej   = 0.5*trace2(Dab,Ja+Jb)
        Exc = 0.5*exca + 0.5*excb
        
        energy = Eone + Ej + Exc + enuke
        
        print i+1,"   ",energy,"    ",Ej,"  ",Exc
        logging.debug("%d %10.4f %10.4f %10.4f" %
                  (i,energy,Ej,Exc))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final U%s energy for system %s is %f"
                 % (functional,atoms.name,energy))
    return energy,(orbea,orbeb),(orbsa,orbsb)