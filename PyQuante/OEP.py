#!/usr/bin/env python

# Yang/Wu's OEP implementation, in PyQuante.

from math import sqrt
from PyQuante.NumWrap import zeros,matrixmultiply,transpose,dot,identity,\
     array,solve_linear_equations
from PyQuante.Ints import getbasis, getints, getJ,get2JmK
from PyQuante.LA2 import GHeigenvectors,mkdens,TraceProperty
from PyQuante.hartree_fock import get_fock
from PyQuante.CGBF import three_center
from PyQuante.optimize import fminBFGS
from PyQuante.fermi_dirac import get_efermi, get_fermi_occs,mkdens_occs,\
     get_entropy
from PyQuante import logging

gradcall=0

def exx(atoms,orbs,**opts):
    return oep_hf(atoms,orbs,**opts)

def oep_hf(atoms,orbs,**opts):
    """oep_hf - Form the optimized effective potential for HF exchange.
       See notes on options and other args in oep routine.
    """
    return oep(atoms,orbs,get_exx_energy,get_exx_gradient,**opts)

def oep(atoms,orbs,energy_func,grad_func=None,**opts):
    """oep - Form the optimized effective potential for a given energy expression

    oep(atoms,orbs,energy_func,grad_func=None,**opts)

    atoms       A Molecule object containing a list of the atoms
    orbs        A matrix of guess orbitals
    energy_func The function that returns the energy for the given method
    grad_func   The function that returns the force for the given method

    Options
    -------
    verbose       False   Output terse information to stdout (default)
                  True    Print out additional information 
    ETemp         False   Use ETemp value for finite temperature DFT (default)
                  float   Use (float) for the electron temperature
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    """
    verbose = opts.get('verbose',False)
    ETemp = opts.get('ETemp',False)
    opt_method = opts.get('opt_method','BFGS')

    bfs = opts.get('bfs',None)
    if not bfs:
        basis = opts.get('basis',None)
        bfs = getbasis(atoms,basis)

    # The basis set for the potential can be set different from
    #  that used for the wave function
    pbfs = opts.get('pbfs',None) 
    if not pbfs: pbfs = bfs
    npbf = len(pbfs)

    integrals = opts.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    nocc,nopen = atoms.get_closedopen()

    Enuke = atoms.get_enuke()

    # Form the OEP using Yang/Wu, PRL 89 143002 (2002)
    nbf = len(bfs)
    norb = nbf
    bp = zeros(nbf,'d')

    bvec = opts.get('bvec',None)
    if bvec:
        assert len(bvec) == npbf
        b = array(bvec)
    else:
        b = zeros(npbf,'d')


    # Form and store all of the three-center integrals
    # we're going to need.
    # These are <ibf|gbf|jbf> (where 'bf' indicates basis func,
    #                          as opposed to MO)
    # N^3 storage -- obviously you don't want to do this for
    #  very large systems
    Gij = []
    for g in range(npbf):
        gmat = zeros((nbf,nbf),'d')
        Gij.append(gmat)
        gbf = pbfs[g]
        for i in range(nbf):
            ibf = bfs[i]
            for j in range(i+1):
                jbf = bfs[j]
                gij = three_center(ibf,gbf,jbf)
                gmat[i,j] = gij
                gmat[j,i] = gij

    # Compute the Fermi-Amaldi potential based on the LDA density.
    # We're going to form this matrix from the Coulombic matrix that
    # arises from the input orbitals. D0 and J0 refer to the density
    # matrix and corresponding Coulomb matrix
    
    D0 = mkdens(orbs,0,nocc)
    J0 = getJ(Ints,D0)
    Vfa = (2*(nel-1.)/nel)*J0
    H0 = h + Vfa

    b = fminBFGS(energy_func,b,grad_func,
                 (nbf,nel,nocc,ETemp,Enuke,S,h,Ints,H0,Gij),
                 logger=logging)

    energy,orbe,orbs = energy_func(b,nbf,nel,nocc,ETemp,Enuke,
                                   S,h,Ints,H0,Gij,return_flag=1)
    return energy,orbe,orbs


def get_exx_energy(b,nbf,nel,nocc,ETemp,Enuke,S,h,Ints,H0,Gij,**opts):
    """Computes the energy for the OEP/HF functional

    Options:
    return_flag    0   Just return the energy
                   1   Return energy, orbe, orbs
                   2   Return energy, orbe, orbs, F
    """
    return_flag = opts.get('return_flag',0)
    Hoep = get_Hoep(b,H0,Gij)
    orbe,orbs = GHeigenvectors(Hoep,S)
        
    if ETemp:
        efermi = get_efermi(nel,orbe,ETemp)
        occs = get_fermi_occs(efermi,orbe,ETemp)
        D = mkdens_occs(orbs,occs)
        entropy = get_entropy(occs,ETemp)
    else:
        D = mkdens(orbs,0,nocc)
        
    F = get_fock(D,Ints,h)
    energy = TraceProperty(h+F,D)+Enuke
    if ETemp: energy += entropy
    iref = nel/2
    gap = 627.51*(orbe[iref]-orbe[iref-1])

    logging.debug("EXX Energy, B, Gap: %10.5f %10.5f %10.5f"
                  % (energy,sqrt(dot(b,b)),gap))
    #logging.debug("%s" % orbe)
    if return_flag == 1:
        return energy,orbe,orbs
    elif return_flag == 2:
        return energy,orbe,orbs,F
    return energy

def get_exx_gradient(b,nbf,nel,nocc,ETemp,Enuke,S,h,Ints,H0,Gij,**opts):
    """Computes the gradient for the OEP/HF functional.

    return_flag    0   Just return gradient
                   1   Return energy,gradient 
                   2   Return energy,gradient,orbe,orbs 
    """
    # Dump the gradient every 10 steps so we can restart...
    global gradcall
    gradcall += 1
    #if gradcall % 5 == 0: logging.debug("B vector:\n%s" % b)

    # Form the new potential and the new orbitals
    energy,orbe,orbs,F = get_exx_energy(b,nbf,nel,nocc,ETemp,Enuke,
                                        S,h,Ints,H0,Gij,return_flag=2)

    Fmo = matrixmultiply(orbs,matrixmultiply(F,transpose(orbs)))

    norb = nbf
    bp = zeros(nbf,'d') # dE/db

    for g in range(nbf):
        # Transform Gij[g] to MOs. This is done over the whole
        #  space rather than just the parts we need. I can speed
        #  this up later by only forming the i,a elements required
        Gmo = matrixmultiply(orbs,matrixmultiply(Gij[g],
                                                 transpose(orbs)))

        # Now sum the appropriate terms to get the b gradient
        for i in range(nocc):
            for a in range(nocc,norb):
                bp[g] = bp[g] + Fmo[i,a]*Gmo[i,a]/(orbe[i]-orbe[a])

    #logging.debug("EXX  Grad: %10.5f" % (sqrt(dot(bp,bp))))
    return_flag = opts.get('return_flag',0)
    if return_flag == 1:
        return energy,bp
    elif return_flag == 2:
        return energy,bp,orbe,orbs
    return bp

def get_Hoep(b,H0,Gij):
    Hoep = H0
    # Add the contributions from the gaussian potential functions
    # H[ij] += b[g]*<ibf|g|jbf>
    for g in range(len(b)):
        Hoep = Hoep + b[g]*Gij[g]
    return Hoep

# Here's a much faster way to do this. Haven't figured out how to
#  do it for more generic functions like OEP-GVB
def oep_hf_an(atoms,orbs,**opts):
    """oep_hf - Form the optimized effective potential for HF exchange.

    Implementation of Wu and Yang's Approximate Newton Scheme
    from J. Theor. Comp. Chem. 2, 627 (2003).

    oep_hf(atoms,orbs,**opts)

    atoms       A Molecule object containing a list of the atoms
    orbs        A matrix of guess orbitals

    Options
    -------
    bfs           None    The basis functions to use for the wfn
    pbfs          None    The basis functions to use for the pot
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    """
    maxiter = opts.get('maxiter',100)
    tol = opts.get('tol',1e-5)
    bfs = opts.get('bfs',None)
    if not bfs:
        basis = opts.get('basis',None)
        bfs = getbasis(atoms,basis)

    # The basis set for the potential can be set different from
    #  that used for the wave function
    pbfs = opts.get('pbfs',None) 
    if not pbfs: pbfs = bfs
    npbf = len(pbfs)

    integrals = opts.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    nocc,nopen = atoms.get_closedopen()

    Enuke = atoms.get_enuke()

    # Form the OEP using Yang/Wu, PRL 89 143002 (2002)
    nbf = len(bfs)
    norb = nbf
    bp = zeros(nbf,'d')

    bvec = opts.get('bvec',None)
    if bvec:
        assert len(bvec) == npbf
        b = array(bvec)
    else:
        b = zeros(npbf,'d')


    # Form and store all of the three-center integrals
    # we're going to need.
    # These are <ibf|gbf|jbf> (where 'bf' indicates basis func,
    #                          as opposed to MO)
    # N^3 storage -- obviously you don't want to do this for
    #  very large systems
    Gij = []
    for g in range(npbf):
        gmat = zeros((nbf,nbf),'d')
        Gij.append(gmat)
        gbf = pbfs[g]
        for i in range(nbf):
            ibf = bfs[i]
            for j in range(i+1):
                jbf = bfs[j]
                gij = three_center(ibf,gbf,jbf)
                gmat[i,j] = gij
                gmat[j,i] = gij

    # Compute the Fermi-Amaldi potential based on the LDA density.
    # We're going to form this matrix from the Coulombic matrix that
    # arises from the input orbitals. D0 and J0 refer to the density
    # matrix and corresponding Coulomb matrix
    
    D0 = mkdens(orbs,0,nocc)
    J0 = getJ(Ints,D0)
    Vfa = (2*(nel-1.)/nel)*J0
    H0 = h + Vfa

    b = zeros(nbf,'d')
    eold = 0

    for iter in range(maxiter):
        Hoep = get_Hoep(b,H0,Gij)
        orbe,orbs = GHeigenvectors(Hoep,S)
        
        D = mkdens(orbs,0,nocc)
        Vhf = get2JmK(Ints,D)

        energy = TraceProperty(2*h+Vhf,D)+Enuke
        if abs(energy-eold) < tol:
            break
        else:
            eold = energy
        
        logging.debug("OEP AN Opt: %d %f" % (iter,energy))
        dV_ao = Vhf-Vfa
        dV = matrixmultiply(orbs,matrixmultiply(dV_ao,transpose(orbs)))

        X = zeros((nbf,nbf),'d')
        c = zeros(nbf,'d')
        Gkt = zeros((nbf,nbf),'d')

        for k in range(nbf):
            Gk = matrixmultiply(orbs,matrixmultiply(Gij[k],
                                                    transpose(orbs)))
            for i in range(nocc):
                for a in range(nocc,norb):
                    c[k] += dV[i,a]*Gk[i,a]/(orbe[i]-orbe[a])
                    
            for l in range(nbf):
                Gl = matrixmultiply(orbs,matrixmultiply(Gij[l],
                                                        transpose(orbs)))
                for i in range(nocc):
                    for a in range(nocc,norb):
                        X[k,l] += Gk[i,a]*Gl[i,a]/(orbe[i]-orbe[a])
        # This should actually be a pseudoinverse...
        b = solve_linear_equations(X,c)

    logging.info("Final OEP energy = %f" % energy)
    return energy,orbe,orbs

def oep_uhf_an(atoms,orbsa,orbsb,**opts):
    """oep_hf - Form the optimized effective potential for HF exchange.

    Implementation of Wu and Yang's Approximate Newton Scheme
    from J. Theor. Comp. Chem. 2, 627 (2003).

    oep_uhf(atoms,orbs,**opts)

    atoms       A Molecule object containing a list of the atoms
    orbs        A matrix of guess orbitals

    Options
    -------
    bfs           None    The basis functions to use for the wfn
    pbfs          None    The basis functions to use for the pot
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    """
    maxiter = opts.get('maxiter',100)
    tol = opts.get('tol',1e-5)
    ETemp = opts.get('ETemp',False)
    bfs = opts.get('bfs',None)
    if not bfs:
        basis = opts.get('basis',None)
        bfs = getbasis(atoms,basis)

    # The basis set for the potential can be set different from
    #  that used for the wave function
    pbfs = opts.get('pbfs',None) 
    if not pbfs: pbfs = bfs
    npbf = len(pbfs)

    integrals = opts.get('integrals',None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    nclosed,nopen = atoms.get_closedopen()
    nalpha,nbeta = nclosed+nopen,nclosed

    Enuke = atoms.get_enuke()

    # Form the OEP using Yang/Wu, PRL 89 143002 (2002)
    nbf = len(bfs)
    norb = nbf

    ba = zeros(npbf,'d')
    bb = zeros(npbf,'d')

    # Form and store all of the three-center integrals
    # we're going to need.
    # These are <ibf|gbf|jbf> (where 'bf' indicates basis func,
    #                          as opposed to MO)
    # N^3 storage -- obviously you don't want to do this for
    #  very large systems
    Gij = []
    for g in range(npbf):
        gmat = zeros((nbf,nbf),'d')
        Gij.append(gmat)
        gbf = pbfs[g]
        for i in range(nbf):
            ibf = bfs[i]
            for j in range(i+1):
                jbf = bfs[j]
                gij = three_center(ibf,gbf,jbf)
                gmat[i,j] = gij
                gmat[j,i] = gij

    # Compute the Fermi-Amaldi potential based on the LDA density.
    # We're going to form this matrix from the Coulombic matrix that
    # arises from the input orbitals. D0 and J0 refer to the density
    # matrix and corresponding Coulomb matrix
    
    D0 = mkdens(orbsa,0,nalpha)+mkdens(orbsb,0,nbeta)
    J0 = getJ(Ints,D0)
    Vfa = ((nel-1.)/nel)*J0
    H0 = h + Vfa

    eold = 0

    for iter in range(maxiter):
        Hoepa = get_Hoep(ba,H0,Gij)
        Hoepb = get_Hoep(ba,H0,Gij)

        orbea,orbsa = GHeigenvectors(Hoepa,S)
        orbeb,orbsb = GHeigenvectors(Hoepb,S)

        if ETemp:
            efermia = get_efermi(nalpha,orbea,ETemp)
            occsa = get_fermi_occs(efermia,orbea,ETemp)
            Da = mkdens_occs(orbsa,occsa)
            efermib = get_efermi(nbeta,orbeb,ETemp)
            occsb = get_fermi_occs(efermib,orbeb,ETemp)
            Db = mkdens_occs(orbsb,occsb)
            entropy = 0.5*(get_entropy(occsa,ETemp)+get_entropy(occsb,ETemp))
        else:
            Da = mkdens(orbsa,0,nalpha)
            Db = mkdens(orbsb,0,nbeta)

        J = getJ(Ints,Da) + getJ(Ints,Db)
        Ka = getK(Ints,Da)
        Kb = getK(Ints,Db)

        energy = (TraceProperty(2*h+J-Ka,Da)+TraceProperty(2*h+J-Kb,Db))/2\
                 +Enuke
        if ETemp: energy += entropy
        
        if abs(energy-eold) < tol:
            break
        else:
            eold = energy
        
        logging.debug("OEP AN Opt: %d %f" % (iter,energy))

        # Do alpha and beta separately
        # Alphas
        dV_ao = J-Ka-Vfa
        dV = matrixmultiply(orbsa,matrixmultiply(dV_ao,transpose(orbsa)))
        X = zeros((nbf,nbf),'d')
        c = zeros(nbf,'d')
        for k in range(nbf):
            Gk = matrixmultiply(orbsa,matrixmultiply(Gij[k],
                                                    transpose(orbsa)))
            for i in range(nalpha):
                for a in range(nalpha,norb):
                    c[k] += dV[i,a]*Gk[i,a]/(orbea[i]-orbea[a])
            for l in range(nbf):
                Gl = matrixmultiply(orbsa,matrixmultiply(Gij[l],
                                                        transpose(orbsa)))
                for i in range(nalpha):
                    for a in range(nalpha,norb):
                        X[k,l] += Gk[i,a]*Gl[i,a]/(orbea[i]-orbea[a])
        # This should actually be a pseudoinverse...
        ba = solve_linear_equations(X,c)
        # Betas
        dV_ao = J-Kb-Vfa
        dV = matrixmultiply(orbsb,matrixmultiply(dV_ao,transpose(orbsb)))
        X = zeros((nbf,nbf),'d')
        c = zeros(nbf,'d')
        for k in range(nbf):
            Gk = matrixmultiply(orbsb,matrixmultiply(Gij[k],
                                                    transpose(orbsb)))
            for i in range(nbeta):
                for a in range(nbeta,norb):
                    c[k] += dV[i,a]*Gk[i,a]/(orbeb[i]-orbeb[a])
            for l in range(nbf):
                Gl = matrixmultiply(orbsb,matrixmultiply(Gij[l],
                                                        transpose(orbsb)))
                for i in range(nbeta):
                    for a in range(nbeta,norb):
                        X[k,l] += Gk[i,a]*Gl[i,a]/(orbeb[i]-orbeb[a])
        # This should actually be a pseudoinverse...
        bb = solve_linear_equations(X,c)

    logging.info("Final OEP energy = %f" % energy)
    return energy,(orbea,orbeb),(orbsa,orbsb)

if __name__ == '__main__':
    from PyQuante.Molecule import Molecule
    from PyQuante.Ints import getbasis,getints
    from PyQuante.hartree_fock import rhf
    from PyQuante import logging
    logging.basicConfig(level=logging.DEBUG,format="%(message)s")

    #mol = Molecule('HF',[('H',(0.,0.,0.)),('F',(0.,0.,0.898369))],
    #              units='Angstrom')
    mol = Molecule('LiH',[(1,(0,0,1.5)),(3,(0,0,-1.5))],units = 'Bohr')
    
    bfs = getbasis(mol)
    S,h,Ints = getints(bfs,mol)
    print "after integrals"
    E_hf,orbe_hf,orbs_hf = rhf(mol,bfs=bfs,integrals=(S,h,Ints),DoAveraging=True)
    print "RHF energy = ",E_hf
    E_exx,orbe_exx,orbs_exx = exx(mol,orbs_hf,bfs=bfs,integrals=(S,h,Ints))
