"""\

NAME
      force.py

SYNOPSIS
      Force module for PyQuante electronic structure calculations. Code is
      loosely based upon Szabo and Ostlund's appendix C describing geometry
      optimization and calculating analytic derivatives.  


DESCRIPTION
	  
AUTHOR
      Hatem H. Helal, hhh23@cam.ac.uk

REPORT BUGS
      Report bugs to hhh23@cam.ac.uk

COPYRIGHT

"""

from NumWrap import array,array2string,zeros,reshape,dot
from Ints import getbasis
from LA2 import trace2
from math import sqrt
from PyQuante.cints import ijkl2intindex
from AnalyticDerivatives import der_Hcore_element,der_overlap_element,der_Jints

def hf_force(mol,wf,bname):
# calculates Hartree-Fock derived atomic forces through
# analytic derivatives of the HF energy.  Stores the forces 
# the atom class which can later be accessed through 
# atomlist[j].forces[i] which would give you component i
# of the force on atom j
    bset = getbasis(mol.atoms,bname)
     
    if wf.restricted:
        rhf_force(mol,wf,bset)
    if wf.unrestricted:
        uhf_force(mol,wf,bset)
    if wf.fixedocc:
        fixedocc_uhf_force(mol,wf,bset)
        
    return

def rhf_force(mol,wf,bset):
    #need to check if this still works for restricted
    #open shell calculations
    Dmat = wf.mkdens()
    Qmat = wf.mkQmatrix()
    
    #compute the force on each atom
    for atom in mol.atoms:
        #equation C.12 on page 441 of Szabo
        #Label the terms in this equation as
        #dE/dRa = d(one electron)/dRa + d(two electron)/dRa 
        #        + d(density matrix)/dRa + d(nuclear repulsion)/dRa
        #the names for these terms are probably open for dispute...
        dE_dRa =   2*der_oneeE(atom.atid,Dmat,bset,mol.atoms) \
                 + der_twoeE(atom.atid,Dmat,bset) \
                 - 2*der_dmat(atom.atid,Qmat,bset) \
                 + der_enuke(atom.atid,mol.atoms)

        fa = -dE_dRa

        atom.set_force(fa)

    return
    
def uhf_force(mol,wf,bset):
    Da,Db = wf.mkdens()
    Qa,Qb = wf.mkQmatrix()
    
    for atom in mol.atoms:
        dEa_dR =   der_oneeE(atom.atid,Da,bset,mol.atoms) \
                 - der_dmat(atom.atid,Qa,bset)
                 
        dEb_dR =   der_oneeE(atom.atid,Db,bset,mol.atoms) \
                 - der_dmat(atom.atid,Qb,bset) 

        dtwoe  =  der_twoeE_uhf(atom.atid,Da,Db,bset)

        denuke = der_enuke(atom.atid,mol.atoms)
        
        f = -(dEa_dR + dEb_dR + dtwoe + denuke)
        
        atom.set_force(f)
    return

def fixedocc_uhf_force(mol,wf,bset):
    Da,Db = wf.mk_auger_dens()
    Qa,Qb = wf.mk_auger_Qmatrix()
    
    for atom in mol.atoms:
        dEa_dR =   der_oneeE(atom.atid,Da,bset,mol.atoms) \
                 - der_dmat(atom.atid,Qa,bset)
                 
        dEb_dR =   der_oneeE(atom.atid,Db,bset,mol.atoms) \
                 - der_dmat(atom.atid,Qb,bset) 

        dtwoe  =  der_twoeE_uhf(atom.atid,Da,Db,bset)

        denuke = der_enuke(atom.atid,mol.atoms)
        
        f = -(dEa_dR + dEb_dR + dtwoe + denuke)
        
        atom.set_force(f)
    return
        
def der_oneeE(a,D,bset,atoms):

    dH_dXa,dH_dYa,dH_dZa = der_Hcore_matrix(a,bset,atoms)
    
    doneE_Xa = trace2(D,dH_dXa)
    doneE_Ya = trace2(D,dH_dYa)
    doneE_Za = trace2(D,dH_dZa)
    
    #print doneE_Xa,doneE_Ya,doneE_Za
    return array([doneE_Xa,doneE_Ya,doneE_Za],'d')

def der_twoeE(a,D,bset):
    d2Ints_dXa,d2Ints_dYa,d2Ints_dZa  = der2Ints(a,bset)

    Gx,Gy,Gz = der2JmK(D,d2Ints_dXa,d2Ints_dYa,d2Ints_dZa)
    
    dtwoeE_Xa = trace2(D,Gx)
    dtwoeE_Ya = trace2(D,Gy)
    dtwoeE_Za = trace2(D,Gz)

    #print dtwoeE_Xa,dtwoeE_Ya,dtwoeE_Za
    return array([dtwoeE_Xa,dtwoeE_Ya,dtwoeE_Za],'d')
    
def der_twoeE_uhf(a,Da,Db,bset):
    d2Ints_dXa,d2Ints_dYa,d2Ints_dZa  = der2Ints(a,bset)

    dJax,dJay,dJaz = derJ(Da,d2Ints_dXa,d2Ints_dYa,d2Ints_dZa)
    dJbx,dJby,dJbz = derJ(Db,d2Ints_dXa,d2Ints_dYa,d2Ints_dZa)
    
    dKax,dKay,dKaz = derK(Da,d2Ints_dXa,d2Ints_dYa,d2Ints_dZa)
    dKbx,dKby,dKbz = derK(Db,d2Ints_dXa,d2Ints_dYa,d2Ints_dZa)

    Dab = Da + Db
    
    dtwoeE_Xa = 0.5*trace2(Dab,dJax+dJbx) - 0.5*(trace2(Da,dKax)+trace2(Db,dKbx))
    dtwoeE_Ya = 0.5*trace2(Dab,dJay+dJby) - 0.5*(trace2(Da,dKay)+trace2(Db,dKby))
    dtwoeE_Za = 0.5*trace2(Dab,dJaz+dJbz) - 0.5*(trace2(Da,dKaz)+trace2(Db,dKbz))
    
    return array([dtwoeE_Xa,dtwoeE_Ya,dtwoeE_Za],'d')

def der_dmat(a,Qmat,bset):
    """
    Looking at Szabo's equation C.12 you can see that this term results 
    from adding up the terms that involve derivatives of the density matrix
    elements P_mu,nu.  Following Szabo we compute this term as
      sum_mu,nu Q_mu,nu dS_mu,nu / dRa
    where the Q matrix is essentially a density matrix that is weighted by 
    the orbital eigenvalues.  Computing the derivative of the overlap integrals
    is done in the der_overlap_matrix function later in this module.
    """
    dS_dXa,dS_dYa,dS_dZa = der_overlap_matrix(a,bset)
    
    dDmat_Xa = trace2(Qmat,dS_dXa)
    dDmat_Ya = trace2(Qmat,dS_dYa)
    dDmat_Za = trace2(Qmat,dS_dZa)
    
    #print dDmat_Xa,dDmat_Ya,dDmat_Za
    return array([dDmat_Xa,dDmat_Ya,dDmat_Za],'d')


def der_enuke(a,atoms):
    """
    returns the derivative of the nuclear repulsion energy
    with respect to the atomic coordinates of atom 'a'
    """
    natoms = len(atoms)
    at_a = atoms[a]
    
    denuke_Xa,denuke_Ya,denuke_Za = 0.0,0.0,0.0
    
    for b in xrange(natoms):
        if b!=a:
            at_b = atoms[b]
            coef = at_a.get_nuke_chg()*at_b.get_nuke_chg()/(at_a.dist(at_b))**3
            
            denuke_Xa += coef*(at_b.r[0]-at_a.r[0])
            denuke_Ya += coef*(at_b.r[1]-at_a.r[1])
            denuke_Za += coef*(at_b.r[2]-at_a.r[2])

    #print denuke_Xa,denuke_Ya,denuke_Za

    return array([denuke_Xa,denuke_Ya,denuke_Za],'d')

def der_Hcore_matrix(a,bset,atoms):
    """
    Evaluates the derivative of the single particle Hamiltonian matrix elements.
    This amounts to evaluating the derivatives of the kinetic energy and the 
    nuclear-electron attraction.
    """
    #initialize dHcore/dR matrices
    nbf = len(bset)
    dHcore_dXa = zeros((nbf,nbf),'d')
    dHcore_dYa = zeros((nbf,nbf),'d')
    dHcore_dZa = zeros((nbf,nbf),'d')

    for i in xrange(nbf): #rows
        for j in xrange(nbf): #columns
            dHcore_dXa[i][j],dHcore_dYa[i][j],dHcore_dZa[i][j] = der_Hcore_element(a,bset[i],bset[j],atoms)
    
    #if a==1: print "dH_dXa"; pad_out(dHcore_dXa); print "dH_dYa"; pad_out(dHcore_dYa); print "dH_dZa"; pad_out(dHcore_dZa);
    return dHcore_dXa,dHcore_dYa,dHcore_dZa 

def der_overlap_matrix(a,bset):
    """
    Evaluates the derivative of the overlap integrals
    with respect to atomic coordinate
    """
    #initialize dS/dR matrices
    nbf = len(bset)
    dS_dXa = zeros((nbf,nbf),'d')
    dS_dYa = zeros((nbf,nbf),'d')
    dS_dZa = zeros((nbf,nbf),'d')
    
    for i in xrange(nbf): #rows
        for j in xrange(nbf): #columns
            dS_dXa[i][j],dS_dYa[i][j],dS_dZa[i][j] = der_overlap_element(a,bset[i],bset[j])

    #if a==1: print "dS_dXa"; pad_out(dS_dXa); print "dS_dYa"; pad_out(dS_dYa); print "dS_dZa"; pad_out(dS_dZa);
    return dS_dXa,dS_dYa,dS_dZa

def der2Ints(a,bset):
    #modified from Ints.py -> get2ints
    """Store integrals in a long array in the form (ij|kl) (chemists
    notation. We only need i>=j, k>=l, and ij <= kl"""
    from array import array
    nbf = len(bset)
    totlen = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
    d2Ints_dXa = array('d',[0]*totlen)
    d2Ints_dYa = array('d',[0]*totlen)
    d2Ints_dZa = array('d',[0]*totlen)
    for i in xrange(nbf):
        for j in xrange(i+1):
            ij = i*(i+1)/2+j
            for k in xrange(nbf):
                for l in xrange(k+1):
                    kl = k*(k+1)/2+l
                    if ij >= kl:
                        ijkl = ijkl2intindex(i,j,k,l)
                        d2Ints_dXa[ijkl],d2Ints_dYa[ijkl],d2Ints_dZa[ijkl] =\
                                 der_Jints(a,bset[i],bset[j],bset[k],bset[l])
    return d2Ints_dXa,d2Ints_dYa,d2Ints_dZa

def derJ(D,d2Ints_dXa,d2Ints_dYa,d2Ints_dZa):
    #modified from Ints.py -> getJ
    "Form the Coulomb operator corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    dJx = zeros((nbf,nbf),'d')
    dJy = zeros((nbf,nbf),'d')
    dJz = zeros((nbf,nbf),'d')

    for i in xrange(nbf):
        for j in xrange(i+1):
            xtemp = zeros(nbf*nbf,'d')
            ytemp = zeros(nbf*nbf,'d')
            ztemp = zeros(nbf*nbf,'d')
            kl = 0
            for k in xrange(nbf):
                for l in xrange(nbf):
                    index = ijkl2intindex(i,j,k,l)
                    xtemp[kl] = d2Ints_dXa[index]
                    ytemp[kl] = d2Ints_dYa[index]
                    ztemp[kl] = d2Ints_dZa[index]
                    kl += 1
            dJx[i,j] = dot(xtemp,D1d)
            dJx[j,i] = dJx[i,j]
            
            dJy[i,j] = dot(ytemp,D1d)
            dJy[j,i] = dJy[i,j]
            
            dJz[i,j] = dot(ztemp,D1d)
            dJz[j,i] = dJz[i,j]
    return dJx,dJy,dJz

def derK(D,d2Ints_dXa,d2Ints_dYa,d2Ints_dZa):
    #modified from Ints.py -> getK
    "Form the exchange operator corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    dKx = zeros((nbf,nbf),'d')
    dKy = zeros((nbf,nbf),'d')
    dKz = zeros((nbf,nbf),'d')
    for i in xrange(nbf):
        for j in xrange(i+1):
            xtemp = zeros(nbf*nbf,'d')
            ytemp = zeros(nbf*nbf,'d')
            ztemp = zeros(nbf*nbf,'d')
            kl = 0
            for k in xrange(nbf):
                for l in xrange(nbf):
                    index_k1 = ijkl2intindex(i,k,j,l)
                    index_k2 = ijkl2intindex(i,l,k,j)
                    xtemp[kl] = 0.5*(d2Ints_dXa[index_k1]+d2Ints_dXa[index_k2])
                    ytemp[kl] = 0.5*(d2Ints_dYa[index_k1]+d2Ints_dYa[index_k2])
                    ztemp[kl] = 0.5*(d2Ints_dZa[index_k1]+d2Ints_dZa[index_k2])
                    kl += 1
            dKx[i,j] = dot(xtemp,D1d)
            dKx[j,i] = dKx[i,j]
            
            dKy[i,j] = dot(ytemp,D1d)
            dKy[j,i] = dKy[i,j]
            
            dKz[i,j] = dot(ztemp,D1d)
            dKz[j,i] = dKz[i,j]
    return dKx,dKy,dKz
    
def der2JmK(D,d2Ints_dXa,d2Ints_dYa,d2Ints_dZa):
    #modified from Ints.py -> get2Jmk
    "Form the 2J-K integrals corresponding to a density matrix D"
    nbf = D.shape[0]
    D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
    Gx = zeros((nbf,nbf),'d')
    Gy = zeros((nbf,nbf),'d')
    Gz = zeros((nbf,nbf),'d')
    
    for i in xrange(nbf):
        for j in xrange(i+1):
            xtemp = zeros(nbf*nbf,'d')
            ytemp = zeros(nbf*nbf,'d')
            ztemp = zeros(nbf*nbf,'d')
            kl = 0
            for k in xrange(nbf):
                for l in xrange(nbf):
                    index_j = ijkl2intindex(i,j,k,l)
                    index_k1 = ijkl2intindex(i,k,j,l)
                    index_k2 = ijkl2intindex(i,l,k,j)
                    xtemp[kl] = 2.*d2Ints_dXa[index_j]-0.5*d2Ints_dXa[index_k1]\
                               -0.5*d2Ints_dXa[index_k2]
                               
                    ytemp[kl] = 2.*d2Ints_dYa[index_j]-0.5*d2Ints_dYa[index_k1]\
                               -0.5*d2Ints_dYa[index_k2]
                               
                    ztemp[kl] = 2.*d2Ints_dZa[index_j]-0.5*d2Ints_dZa[index_k1]\
                               -0.5*d2Ints_dZa[index_k2]                    
                    kl += 1

            Gx[i,j] = dot(xtemp,D1d)
            Gx[j,i] = Gx[i,j]
            
            Gy[i,j] = dot(ytemp,D1d)
            Gy[j,i] = Gy[i,j]            
            
            Gz[i,j] = dot(ztemp,D1d)
            Gz[j,i] = Gz[i,j]

    return Gx,Gy,Gz
    
