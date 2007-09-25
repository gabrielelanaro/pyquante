#!/usr/bin/env python
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

from NumWrap import array,array2string,zeros
from Ints import getbasis
from LA2 import trace2
from math import sqrt

from AnalyticDerivatives import der_Hcore_element,der_overlap_element

def hf_forces(mol,wf,bname):
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
        dE_dRa =   der_oneeE(atom.atid,Dmat,bset,mol.atoms) \
                 + der_twoeE(atom.atid,Dmat,bset) \
                 + der_dmat(atom.atid,Qmat,bset) \
                 + der_enuke(atom.atid,mol.atoms)

        fa = -dE_dRa
        
        atom.set_force(fa)

    return
    
def uhf_force(mol,wf,bset):

    return

        
def der_oneeE(a,D,bset,atoms):

    dH_dXa,dH_dYa,dH_dZa = der_Hcore_matrix(a,bset,atoms)
    
    doneE_Xa = 2*trace2(D,dH_dXa)
    doneE_Ya = 2*trace2(D,dH_dYa)
    doneE_Za = 2*trace2(D,dH_dZa)
    
    return array([doneE_Xa,doneE_Ya,doneE_Za],'d')

def der_twoeE(a,D,bset):
    
    dtwoeE_Xa = 0.0
    dtwoeE_Ya = 0.0
    dtwoeE_Za = 0.0
    
    return array([dtwoeE_Xa,dtwoeE_Ya,dtwoeE_Za],'d')

def der_dmat(a,Qmat,bset):
    """
    Looking at Szabo's equation C.12 you can see that this term results 
    from adding up the terms that involve derivatives of the density matrix
    elements P_mu,nu.  Following Szabo we compute this term as
     - sum_mu,nu Q_mu,nu dS_mu,nu / dRa
    where the Q matrix is essentially a density matrix that is weighted by 
    the orbital eigenvalues.  Computing the derivative of the overlap integrals
    is done in the der_overlap_matrix function later in this module.
    """
    dS_dXa,dS_dYa,dS_dZa = der_overlap_matrix(a,bset)
    
    dDmat_Xa = -2*trace2(Qmat,dS_dXa)
    dDmat_Ya = -2*trace2(Qmat,dS_dYa)
    dDmat_Za = -2*trace2(Qmat,dS_dZa)
    
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
    
    for b in range(natoms):
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

    for i in range(nbf): #rows
        for j in range(nbf): #columns
            dHcore_dXa[i][j],dHcore_dYa[i][j],dHcore_dZa[i][j] = der_Hcore_element(a,bset[i],bset[j],atoms)
    
    if a==1: print "dH_dXa"; pad_out(dHcore_dXa); print "dH_dYa"; pad_out(dHcore_dYa); print "dH_dZa"; pad_out(dHcore_dZa);
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
    
    for i in range(nbf): #rows
        for j in range(nbf): #columns
            dS_dXa[i][j],dS_dYa[i][j],dS_dZa[i][j] = der_overlap_element(a,bset[i],bset[j])

    #if a==1: print "dS_dXa"; pad_out(dS_dXa); print "dS_dYa"; pad_out(dS_dYa); print "dS_dZa"; pad_out(dS_dZa);
    return dS_dXa,dS_dYa,dS_dZa

    
    
#below are some helper functions for debugging and testing that are used above and probably 
#should be moved somewhere else or outright removed eventually
def pad_out(matrix):
    #this will make debugging matrix operations easier by getting rid of 
    #matrix elements which are ridiculously tiny
    rows = matrix.shape[0]
    cols = matrix.shape[1]
    
    tolerance = 0.001
    
    for i in range(rows):
        for j in range(cols):
                if abs(matrix[i][j]) < tolerance:
                    matrix[i][j]=0.0
    
    print array2string(matrix,200,precision=4);    print "\n\n"

    return