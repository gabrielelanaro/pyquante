#!/usr/bin/env python
"""\
 MINDO3.py: Dewar's MINDO/3 Semiempirical Method

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from Constants import bohr2ang,e2,ev2kcal
from MINDO3_Parameters import axy,Bxy
from math import sqrt,exp,pow
from NumWrap import zeros
from NumWrap import Heigenvectors
from LA2 import mkdens,TraceProperty
A0 = bohr2ang

def get_beta0(atnoi,atnoj):
    "Resonanace integral for coupling between different atoms"
    return Bxy[(min(atnoi,atnoj),max(atnoi,atnoj))]

def get_alpha(atnoi,atnoj):
    "Part of the scale factor for the nuclear repulsion"
    return axy[(min(atnoi,atnoj),max(atnoi,atnoj))]

def get_gamma(atomi,atomj):
    "Coulomb repulsion that goes to the proper limit at R=0"
    R2 = atomi.dist2(atomj)
    return e2/sqrt(R2+0.25*pow(atomi.rho+atomj.rho,2))

def get_g(bfi,bfj):
    "Coulomb-like term for orbitals on the same atom"
    i,j = bfi.type,bfj.type
    assert bfi.atom is bfj.atom, "Incorrect call to get_g"
    if i==0 and j==0:
        return bfi.atom.gss
    elif i==0 or j==0:
        return bfi.atom.gsp
    elif i==j:
        return bfi.atom.gpp
    return bfi.atom.gppp

def get_h(bfi,bfj):
    "Exchange-like term for orbitals on the same atom"
    i,j = bfi.type,bfj.type
    assert bfi.atom is bfj.atom, "Incorrect call to get_h"
    if i==0 or j==0:
        return bfi.atom.hsp
    return bfi.atom.hppp

def get_nbf(atoms):
    "Number of basis functions in an atom list"
    nbf = 0
    for atom in atoms: nbf += atom.nbf
    return nbf

def get_F0(atoms):
    "Form the zero-iteration (density matrix independent) Fock matrix"
    nbf = get_nbf(atoms)
    nat = len(atoms)
    
    F0 = zeros((nbf,nbf),'d')

    ibf = 0 # bf number of the first bfn on iat
    for iat in range(nat):
        atomi = atoms[iat]
        for i in range(atomi.nbf):
            bfi = atomi.basis[i]
            F0[ibf+i,ibf+i] = bfi.u
            
            jbf = 0
            for jat in range(nat):
                atomj = atoms[jat]
                if iat != jat:
                    gammaij = get_gamma(atomi,atomj)
                    betaij = get_beta0(atomi.atno,atomj.atno)
                    F0[ibf+i,ibf+i] -= gammaij*atomj.Z
                    for j in range(atomj.nbf):
                        bfj = atomj.basis[j]
                        Sij = bfi.cgbf.overlap(bfj.cgbf)
                        #Sij = mopac_overlap(bfi,bfj)
                        IPij = bfi.ip+bfj.ip
                        F0[ibf+i,jbf+j] = betaij*IPij*Sij
                        F0[jbf+j,ibf+i] = F0[ibf+i,jbf+j]
                jbf += atomj.nbf
        ibf += atomi.nbf
    return F0

def get_F1(atoms,D):
    "One-center corrections to the core fock matrix"
    nbf = get_nbf(atoms)
    nat = len(atoms)
    
    F1 = zeros((nbf,nbf),'d')

    ibf = 0 # bf number of the first bfn on iat
    for iat in range(nat):
        atomi = atoms[iat]
        for i in range(atomi.nbf):
            bfi = atomi.basis[i]
            gii = get_g(bfi,bfi)
            qi =  D[ibf+i,ibf+i]
            F1[ibf+i,ibf+i] = 0.5*qi*gii
            
            for j in range(atomi.nbf):  # ij on same atom
                if j != i:
                    bfj = atomi.basis[j]
                    qj = D[ibf+j,ibf+j]
                    gij = get_g(bfi,bfj)
                    pij = D[ibf+i,ibf+j]
                    hij = get_h(bfi,bfj)
                    # the following 0.5 is something of a kludge to match
                    #  the mopac results.
                    F1[ibf+i,ibf+i] += qj*gij - 0.5*qj*hij
                    F1[ibf+i,ibf+j] += 0.5*pij*(3*hij-gij)
        ibf += atomi.nbf
    return F1


def get_F1_open(atoms,Da,Db):
    "One-center corrections to the core fock matrix"
    nbf = get_nbf(atoms)
    nat = len(atoms)
    
    F1 = zeros((nbf,nbf),'d')

    ibf = 0 # bf number of the first bfn on iat
    for iat in range(nat):
        atomi = atoms[iat]
        for i in range(atomi.nbf):
            gii = get_g(atomi.basis[i],atomi.basis[i])
            qib =  Db[ibf+i,ibf+i]
            #electron only interacts with the other electron in orb,
            # not with itself
            F1[ibf+i,ibf+i] = qib*gii 
            
            for j in range(atomi.nbf):  # ij on same atom
                if j != i:
                    qja = Da[ibf+j,ibf+j]
                    qjb = Db[ibf+j,ibf+j]
                    qj = qja+qjb
                    gij = get_g(atomi.basis[i],atomi.basis[j])
                    pija = Da[ibf+i,ibf+j] 
                    pijb = Db[ibf+i,ibf+j] 
                    pij = pija + pijb
                    hij = get_h(atomi.basis[i],atomi.basis[j])
                    # the following 0.5 is something of a kludge to match
                    #  the mopac results.
                    F1[ibf+i,ibf+i] += qj*gij - qja*hij
                    F1[ibf+i,ibf+j] += 2*pij*hij - pija*(hij+gij)
        ibf += atomi.nbf
    return F1

def get_F2(atoms,D):
    "Two-center corrections to the core fock matrix"
    nbf = get_nbf(atoms)
    nat = len(atoms)
    
    F2 = zeros((nbf,nbf),'d')

    ibf = 0 # bf number of the first bfn on iat
    for iat in range(nat):
        atomi = atoms[iat]
        jbf = 0
        for jat in range(nat):
            atomj = atoms[jat]
            if iat != jat:
                gammaij = get_gamma(atomi,atomj)
                for i in range(atomi.nbf):
                    for j in range(atomj.nbf):
                        pij = D[ibf+i,jbf+j]
                        qj = D[jbf+j,jbf+j]
                        qi = D[ibf+i,ibf+i]
                        F2[ibf+i,jbf+j] -= 0.25*pij*gammaij
                        F2[jbf+j,ibf+i] = F2[ibf+i,jbf+j]
                        # The following 0.5 is a kludge
                        F2[ibf+i,ibf+i] += 0.5*qj*gammaij
                        F2[jbf+j,jbf+j] += 0.5*qi*gammaij
            jbf += atomj.nbf
        ibf += atomi.nbf
    return F2

def get_F2_open(atoms,Da,Db):
    "Two-center corrections to the core fock matrix"
    nbf = get_nbf(atoms)
    nat = len(atoms)
    
    F2 = zeros((nbf,nbf),'d')

    ibf = 0 # bf number of the first bfn on iat
    for iat in range(nat):
        atomi = atoms[iat]
        jbf = 0
        for jat in range(nat):
            atomj = atoms[jat]
            if iat != jat:
                gammaij = get_gamma(atomi,atomj)
                for i in range(atomi.nbf):
                    for j in range(atomj.nbf):
                        pija = Da[ibf+i,jbf+j] 
                        pijb = Db[ibf+i,jbf+j] 
                        pij = pija+pijb
                        qja = Da[jbf+j,jbf+j]
                        qjb = Db[jbf+j,jbf+j]
                        qj = qja+qjb
                        qia = Da[ibf+i,ibf+i]
                        qib = Db[ibf+i,ibf+i]
                        qi = qia+qib
                        F2[ibf+i,jbf+j] -= 0.25*pij*gammaij
                        F2[jbf+j,ibf+i] = F2[ibf+i,jbf+j]
                        # The following 0.5 is a kludge
                        F2[ibf+i,ibf+i] += 0.5*qj*gammaij
                        F2[jbf+j,jbf+j] += 0.5*qi*gammaij
            jbf += atomj.nbf
        ibf += atomi.nbf
    return F2

def get_nel(atoms,charge=0):
    "Number of electrons in an atoms. Can be dependent on the charge"
    nel = 0
    for atom in atoms: nel += atom.Z
    return nel-charge

def get_enuke(atoms):
    "Compute the nuclear repulsion energy"
    enuke = 0
    for i in range(len(atoms)):
        atomi = atoms[i]
        for j in range(i):
            atomj = atoms[j]
            R2 = atomi.dist2(atomj)
            R = sqrt(R2)
            scale = get_scale(atomi.atno,atomj.atno,R)
            gammaij = get_gamma(atomi,atomj)
            enuke += atomi.Z*atomj.Z*gammaij \
                     + abs(atomi.Z*atomj.Z*(e2/R-gammaij)*scale)
    return enuke

def get_scale(atnoi,atnoj,R):
    "Prefactor from the nuclear repulsion term"
    alpha = get_alpha(atnoi,atnoj)
    if atnoi == 1:
        if atnoj == 7 or atnoj == 8:
            return alpha*exp(-R)
    elif atnoj == 1:
        if atnoi == 7 or atnoi == 8:
            return alpha*exp(-R)
    return exp(-alpha*R)

def get_guess_D(atoms):
    "Average occupation density matrix"
    nbf = get_nbf(atoms)
    D = zeros((nbf,nbf),'d')
    ibf = 0
    for atom in atoms:
        atno = atom.atno
        for i in range(atom.nbf):
            if atno == 1:
                D[ibf+i,ibf+i] = atom.Z/1.
            else:                      
                D[ibf+i,ibf+i] = atom.Z/4.
        ibf += atom.nbf
    return D    

def get_reference_energy(atoms):
    "Ref = heat of formation - energy of atomization"
    eat = 0
    hfat = 0
    for atom in atoms:
        eat += atom.Eref
        hfat += atom.Hf
    return hfat-eat*ev2kcal

def get_open_closed(nel,mult=None):
    "Get the number of open/closed orbitals based on nel & multiplicity"
    nclosed,nopen = divmod(nel,2)
    if mult: #test the multiplicity
        nopen = mult-1
        nclosed,ntest = divmod(nel-nopen,2)
        if ntest: raise "Impossible nel, multiplicity %d %d " % (nel,mult)
    return nclosed,nopen

def scfclosed(atoms,F0,nclosed,**opts):
    "SCF procedure for closed-shell molecules"
    verbose = opts.get('verbose',False)
    D = get_guess_D(atoms)
    Eold = 0
    for i in range(10):
        F1 = get_F1(atoms,D)
        F2 = get_F2(atoms,D)
        F = F0+F1+F2
        Eel = 0.5*TraceProperty(D,F0+F)
        if verbose: print i,Eel
        if abs(Eel-Eold) < 0.001: break
        Eold = Eel
        orbe,orbs = Heigenvectors(F)
        D = 2*mkdens(orbs,0,nclosed)
    return Eel

def scfopen(atoms,F0,nalpha,nbeta,**opts):
    "SCF procedure for open-shell molecules"
    verbose = opts.get('verbose',False)
    D = get_guess_D(atoms)
    Da = 0.5*D
    Db = 0.5*D
    Eold = 0
    for i in range(10):
        F1a = get_F1_open(atoms,Da,Db)
        F1b = get_F1_open(atoms,Db,Da)
        F2a = get_F2_open(atoms,Da,Db)
        F2b = get_F2_open(atoms,Db,Da)
        Fa = F0+F1a+F2a
        Fb = F0+F1b+F2b
        Eel = 0.5*TraceProperty(Da,F0+Fa)+0.5*TraceProperty(Db,F0+Fb)
        if verbose: print i,Eel
        if abs(Eel-Eold) < 0.001: break
        Eold = Eel
        orbea,orbsa = Heigenvectors(Fa)
        orbeb,orbsb = Heigenvectors(Fb)
        Da = mkdens(orbsa,0,nalpha)
        Db = mkdens(orbsb,0,nbeta)
    return Eel

def initialize(atoms):
    "Assign parameters for the rest of the calculation"
    from Slater import gauss_powers,gexps,gcoefs,s_or_p
    from MINDO3_Parameters import Uss,Upp,IPs,IPp,CoreQ,f03,nbfat,\
         zetas,zetap,Eat,Hfat,gss,gsp,gpp,gppp,hsp,hppp,NQN
    from CGBF import CGBF
    from Bunch import Bunch # Generic object to hold basis functions
    ibf = 0 # Counter to overall basis function count
    for atom in atoms:
        x,y,z = atom.pos()
        xbohr,ybohr,zbohr = x/bohr2ang,y/bohr2ang,z/bohr2ang
        atom.Z = CoreQ[atom.atno]
        atom.basis = []
        atom.rho = e2/f03[atom.atno]
        atom.nbf = nbfat[atom.atno]
        atom.Eref = Eat[atom.atno]
        atom.Hf = Hfat[atom.atno]
        atom.gss = gss[atom.atno]
        atom.gsp = gsp[atom.atno]
        atom.gpp = gpp[atom.atno]
        atom.gppp = gppp[atom.atno]
        atom.hsp = hsp[atom.atno]
        atom.hppp = hppp[atom.atno]
        for i in range(atom.nbf):
            bfunc = Bunch()
            atom.basis.append(bfunc)
            bfunc.index = ibf # pointer to overall basis function index
            ibf += 1
            bfunc.type = i # s,x,y,z
            bfunc.atom = atom # pointer to parent atom
            bfunc.cgbf = CGBF((xbohr,ybohr,zbohr),gauss_powers[i])
            zi = gexps[(NQN[atom.atno],s_or_p[i])]
            ci = gcoefs[(NQN[atom.atno],s_or_p[i])]
            if i:
                zeta = zetap[atom.atno]
                bfunc.u = Upp[atom.atno]
                bfunc.ip = IPp[atom.atno]
            else:
                zeta = zetas[atom.atno]
                bfunc.u = Uss[atom.atno]
                bfunc.ip = IPs[atom.atno]
            for j in range(len(zi)):
                bfunc.cgbf.add_primitive(zi[j]*zeta*zeta,ci[j])
            bfunc.cgbf.normalize()
    return atoms

def get_fock(atoms):
    "Just return the 0th iteration fock matrix"
    atoms = initialize(atoms)
    F0 = get_F0(atoms)
    D = get_guess_D(atoms)
    F1 = get_F1(atoms,D)
    F2 = get_F2(atoms,D)
    return F0+F1+F2

def scf(atoms,**opts):
    "Driver routine for energy calculations"
    chg = opts.get('chg',0)
    mult = opts.get('mult',None)
    verbose = opts.get('verbose',False)
    atoms = initialize(atoms)
    
    nel = get_nel(atoms)-int(chg)
    nclosed,nopen = get_open_closed(nel,mult)

    Enuke = get_enuke(atoms)
    nbf = get_nbf(atoms)
    eref = get_reference_energy(atoms)
    if verbose:
        print "Nel = %d, Nclosed = %d, Nopen = %d," % (nel,nclosed,nopen), \
              "Enuke = %10.4f, Nbf = %d" % (Enuke,nbf)
    F0 = get_F0(atoms)
    if nopen:
        Eel = scfopen(atoms,F0,nclosed+nopen,nclosed,**opts)
    else:
        Eel = scfclosed(atoms,F0,nclosed,**opts)
    Etot = Eel+Enuke
    Hf = Etot*ev2kcal+eref
    if verbose: print "Final Heat of Formation = ",Hf
    return Hf

def get_energy_forces(atoms,**opts):
    from Convergence import SimpleAverager

    chg = opts.get('chg',0)
    averaging = opts.get('averaging',True)
    verbose = opts.get('verbose',False)
    
    atoms = initialize(atoms)
    nel = get_nel(atoms)-int(chg)
    nclosed,nopen = get_open_closed(nel,None)
    assert nopen==0, "Forces only for closed-shell now"
    Enuke = get_enuke(atoms)
    nbf = get_nbf(atoms)
    eref = get_reference_energy(atoms)
    F0 = get_F0(atoms)
    D = get_guess_D(atoms)
    Eold = 0
    energies = []
    converger = SimpleAverager(averaging)
    for i in range(50):
        D = converger.getD(D)
        F1 = get_F1(atoms,D)
        F2 = get_F2(atoms,D)
        F = F0+F1+F2
        Eel = 0.5*TraceProperty(D,F0+F)
        if abs(Eel-Eold) < 0.001: break
        if verbose: print i,Eel
        energies.append(Eel)
        Eold = Eel
        orbe,orbs = Heigenvectors(F)
        D = 2*mkdens(orbs,0,nclosed)
    else:
        raise "SCF Not Converged: exiting"
    Etot = Eel+Enuke
    Hf = Etot*ev2kcal+eref
    Forces = forces(atoms,D)
    return Hf,Forces

def forces(atoms,D):
    "Compute analytic forces on list of atoms"
    nat = len(atoms)
    Forces = zeros((nat,3),'d')
    # Loop over all pairs of atoms and compute the force between them
    for iat in range(nat):
        atomi = atoms[iat]
        for jat in range(iat):
            atomj = atoms[jat]
            alpha = get_alpha(atomi.atno,atomj.atno)
            beta = get_beta0(atomi.atno,atomj.atno)
            R2 = atomi.dist2(atomj)
            R = sqrt(R2)
            c2 = 0.25*pow(atomi.rho+atomj.rho,2)

            for dir in range(3):
                Fij = 0 # Force between atoms iat and jat in direction dir
                # initialize some constants
                delta = atomi.r[dir]-atomj.r[dir]
                c1 = delta*atomi.Z*atomj.Z*e2/R
                dr1 = e2*delta*pow(R2+c2,-1.5)

                # Nuclear repulsion terms
                if ( (atomi.atno == 1
                      and (atomj.atno == 7 or atomj.atno == 8))
                     or (atomj.atno == 1
                         and (atomi.atno == 7 or atomi.atno == 8))):
                    # Special case of NH or OH bonds
                    Fij += -c1*alpha*(1/R2 - R*pow(R2+c2,-1.5)
                                      + 1/R - 1/sqrt(R2+c2))*exp(-R) \
                                      - c1*R*pow(R2+c2,-1.5)
                else:
                    Fij += -c1*(1/R2 - R*pow(R2+c2,-1.5) + alpha/R 
                                - alpha/sqrt(R2+c2))*exp(-alpha*R) \
                                - c1*R*pow(R2+c2,-1.5)

                # Overlap terms
                for bfi in atomi.basis:
                    for bfj in atomj.basis:
                        Dij = D[bfi.index,bfj.index]
                        #dSij = mopac_doverlap(bfi,bfj,dir)
                        dSij = -bfi.cgbf.doverlap(bfj.cgbf,dir)/bohr2ang
                        #dSij = -bfi.cgbf.doverlap_num(bfj.cgbf,dir)/bohr2ang
                        Fij += 2*beta*(bfi.ip+bfj.ip)*Dij*dSij

                # Core attraction terms
                for bfj in atomj.basis:
                    Fij += atomi.Z*D[bfj.index,bfj.index]*dr1
                for bfi in atomi.basis:
                    Fij += atomj.Z*D[bfi.index,bfi.index]*dr1

                # Two-electron terms
                for bfi in atomi.basis:
                    for bfj in atomj.basis:
                        Dii = D[bfi.index,bfi.index]
                        Djj = D[bfj.index,bfj.index]
                        Dij = D[bfi.index,bfj.index]
                        
                        # exchange is the first term, coulomb is second:
                        Fij += 0.5*dr1*pow(Dij,2)-dr1*Dii*Djj

                # Now sum total forces and convert to kcal/mol
                Forces[iat][dir] += ev2kcal*Fij
                Forces[jat][dir] -= ev2kcal*Fij
    return Forces

def mopac_overlap(bfi,bfj): # from the routine gover.f
    cgbfi,cgbfj = bfi.cgbf,bfj.cgbf
    ri = cgbfi.origin() # distance in bohr
    rj = cgbfj.origin()
    RR = pow(ri[0]-rj[0],2)+pow(ri[1]-rj[1],2)+pow(ri[2]-rj[2],2)
    itype = bfi.type
    jtype = bfj.type
    Sij = 0
    for primi in cgbfi.prims():
        for primj in cgbfj.prims():

            amb = primi.exp()+primj.exp()
            apb = primi.exp()*primj.exp()
            adb = apb/amb

            if itype > 0 and jtype > 0:
                #is = 4
                tomb = (ri[itype-1]-rj[itype-1])*(ri[jtype-1]-rj[jtype-1])
                abn = -adb*tomb
                if itype == jtype: abn = abn + 0.5
                abn = 4*abn*sqrt(apb)/amb
            elif itype > 0:
                #is = 3
                tomb = (ri[itype-1]-rj[itype-1])
                abn = -2*tomb*primj.exp()*sqrt(primi.exp())/amb
            elif jtype > 0:
                #is = 2
                tomb = (ri[jtype-1]-rj[jtype-1])
                abn = 2*tomb*primi.exp()*sqrt(primj.exp())/amb
            else:
                #is = 1
                abn = 1.0
                
            if adb*RR < 90:
                Sij += primi.coef()*primj.coef()*\
                       pow(2*sqrt(apb)/amb,1.5)*exp(-adb*RR)*abn
    return Sij

def mopac_doverlap(bfi,bfj,dir): # from the routine dcart.f
    cgbfi,cgbfj = bfi.cgbf,bfj.cgbf
    ri = cgbfi.origin() # distance in bohr
    rj = cgbfj.origin()
    RR = pow(ri[0]-rj[0],2)+pow(ri[1]-rj[1],2)+pow(ri[2]-rj[2],2)
    del1 = ri[dir] - rj[dir]
    itype = bfi.type
    jtype = bfj.type
    DS = 0
    for primi in cgbfi.prims():
        for primj in cgbfj.prims():
            del2 = del3 = 0
            SS = 0
            apb = primi.exp()*primj.exp()
            amb = primi.exp()+primj.exp()
            adb = apb/amb
            adr = min(adb*RR,35.0)
            if itype == 0 and jtype == 0: # ss
                # is=1
                abn = -2.*adb*del1/A0
            elif itype == 0 and jtype > 0: # sp
                if jtype-1 == dir: 
                    #is = 3
                    abn = 2*adb/sqrt(primj.exp())*(1-2*adb*del1*del1)/A0
                else:
                    #is = 2
                    del2 = ri[jtype-1]-rj[jtype-1]
                    abn = -4*adb*adb*del1*del2/sqrt(primj.exp())/A0
            elif itype > 0 and jtype == 0: # ps
                if itype-1 == dir: 
                    #is = 5
                    abn = -2*adb/sqrt(primi.exp())*(1-2*adb*del1*del1)/A0
                else:
                    #is = 4
                    del2 = ri[itype-1]-rj[itype-1]
                    abn = 4*adb*adb*del1*del2/sqrt(primi.exp())/A0
            elif itype == jtype: 
                if dir == itype-1:
                    #is = 9 (p|p)
                    abn=-8*adb*adb*del1/sqrt(apb)*(1.5-adb*del1*del1)/A0
                else:
                    #is = 8 (p'|p')
                    del2 = ri[jtype-1]-rj[jtype-1]
                    abn=-8*pow(adb,2)*del1/sqrt(apb)*(0.5-adb*del2*del2)/A0
            elif (dir != itype-1) and (dir != jtype-1):
                #is = 7(p'|p")
                del2 = ri[itype-1] - rj[itype-1]
                del3 = ri[jtype-1] - rj[jtype-1]
                abn=8*pow(adb,3)*del1*del2*del3/sqrt(apb)/A0
            else:
                #is = 6 (p|p') or (p'|p)
                del2 = ri[itype+jtype-dir-2]-rj[itype+jtype-dir-2]
                abn=-4*adb*adb*del2/sqrt(apb)*(1-2*adb*del1*del1)/A0
            SS = pow(2*sqrt(apb)/amb,1.5)*exp(-adr)*abn
            DS += SS*primi.coef()*primj.coef()
    return DS

def test_olap():
    # Test function to compare results of my CGBF overlap routines to those
    #  of mopacs. The issue is that the derivative gives different results.
    from math import sin,cos
    from copy import deepcopy

    delta = 0.001
    for theta in [0.,10.,20.,30.,45.,55.214134,90.]:
        at1 = (1,(0,0,0))
        at2 = (6,(cos(theta),sin(theta),0.1))
        atoms = initialize([at1,at2])
        bfi = atoms[0].basis[0]
        bfj = atoms[1].basis[2]

        dSijx = mopac_doverlap(bfi,bfj,0)
        dSijy = mopac_doverlap(bfi,bfj,1)
        dSijz = mopac_doverlap(bfi,bfj,2)

        dSijx2 = -bfi.cgbf.doverlap(bfj.cgbf,0)/bohr2ang
        dSijy2 = -bfi.cgbf.doverlap(bfj.cgbf,1)/bohr2ang
        dSijz2 = -bfi.cgbf.doverlap(bfj.cgbf,2)/bohr2ang
                  
        dSijx4 = -bfi.cgbf.doverlap_num(bfj.cgbf,0)/bohr2ang
        dSijy4 = -bfi.cgbf.doverlap_num(bfj.cgbf,1)/bohr2ang
        dSijz4 = -bfi.cgbf.doverlap_num(bfj.cgbf,2)/bohr2ang
        print "%2d %6.3f %6.3f %6.3f  %6.3f %6.3f %6.3f  " %\
              (theta,dSijx,dSijy,dSijz,dSijx2,dSijy2,dSijz2)
    return

def write_mopac_input(atoms,fname=None):
    from PyQuante.Element import symbol
    if not fname: fname = atoms.name + ".dat"
    lines = ['MINDO3',atoms.name,'Input file written by PyQuante']
    for atom in atoms:
        atno = atom.atno
        sym = symbol[atno]
        x,y,z = atom.r
        lines.append('%s  %10.4f  0  %10.4f  0  %10.4f  0'
                     % (sym,x,y,z))
    open(fname,'w').write('\n'.join(lines))
    return

if __name__ == '__main__':
    from Molecule import Molecule
    h2o = Molecule('H2O',atomlist=[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))])
    oh = Molecule('OH',atomlist=[(8,(0,0,0)),(1,(1.,0,0))])
    ch4 = Molecule('Methane', atomlist =
                   [(6,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0)),
                    (1,(0,0,1.)),(1,(0,0,-1.))])
    print scf(h2o)
    print scf(oh)
    print scf(ch4)
    #E,F = get_energy_forces(ch4)
    #for Fi in F: print Fi
    #import profile,pstats
    #profile.run('get_energy_forces(ch4)','prof')
    #prof = pstats.Stats('prof')
    #prof.strip_dirs().sort_stats('time').print_stats(15)
    #test_olap()
    
