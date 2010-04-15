"""\
 Wavefunction.py: trial module defining an orbital class
 The idea is to make handling computations over orbitals slightly easier

 Started by Hatem H Helal hhh23@cam.ac.uk on 17/07/2007

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from NumWrap import matrixmultiply,transpose,dot
from LA2 import diagonal_mat

class Wavefunction:
    """\
    Data class to hold all relevant wavefunction data in PyQuante

    Example of how to use this in a restricted HF (or DFT or ....) calculation:
    Wavefunction(orbs=orbital_coef_matrix,orbe=orbital_eigenvalues, \
                nclosed=closed_shell_orbitals, nopen=open_shell_orbitals) 

    The unrestricted analog:
    Wavefunction(orbs_a=spin_up_coef_matrix, orbs_b=spin_down_coef_matrix, \
            orbe_a=spin_up_orbital_eigenvalues, orbe_b=spin_down_orbital_eigenvalues, \
            nalpha=spin_up_electrons, nbeta=spin_down_electrons)
    
    This class also allows for a container for orbital occupation numbers which
    could be useful for caclulations involving non-standard occupations, excitations,
    or thermal smearing through a Fermi-Dirac function.

    Options:       Description
    --------       -----------
Restricted
    orbs            Orbital coefficient matrix, columns are eigenstates of 
                    the Fock matrix arranged in ascending order by eigenvalue

    orbe            Orbital eigenvalues in list form
    nclosed         Number of closed shell orbitals
    nopen           Number of open shell orbitals
    occs            Orbital occupation numbers, useful for non-standard occupations

Unrestricted 
    orbs_a          Spin up orbital coefficient matrix, columns are the eigenstates
                    of the spin up Fock matrix arranged in ascending order by eigenvalue
    orbs_b          Spin down analog of orbs_a
    
    orbe_a          Spin up orbital eigenvalues in list form
    orbe_b          Spin down analog of orbe_a

    nalpha          Number of spin up electrons
    nbeta           Number of spin down electrons

    occs_a          Spin up occupation numbers
    occs_b          Spin down analog of occs_a
    """
    def __init__(self,**opts):
        #orbital coefs, energies, and occupations for restricted calculations
        self.orbs    = opts.get('orbs',None)
        self.orbe    = opts.get('orbe',None)
        self.nclosed = opts.get('nclosed',None)
        self.nopen   = opts.get('nopen',None)
        self.occs    = opts.get('occs',None)
        
        #orbital coefs, energies, and occupations for unrestricted calculations        
        self.orbs_a = opts.get('orbs_a',None)
        self.orbs_b = opts.get('orbs_b',None)
        
        self.orbe_a = opts.get('orbe_a',None)
        self.orbe_b = opts.get('orbe_b',None)
        
        self.nalpha = opts.get('nalpha',None)
        self.nbeta  = opts.get('nbeta',None)
        
        self.occs_a = opts.get('occs_a',None)
        self.occs_b = opts.get('occs_b',None)
        
        self.restricted   = opts.get('restricted',None)
        self.unrestricted = opts.get('unrestricted',None)
        self.fixedocc    = opts.get('fixedocc',None)
        return
        
    def update_wf(self,**opts):
        #use this function to update coef and eigenvalues within SCF loop
        #restricted
        self.orbs    = opts.get(orbs,None)
        self.orbe    = opts.get(orbe,None)

        #unrestricted
        self.orbs_a = opts.get(orbs_a,None)
        self.orbs_b = opts.get(orbs_b,None)
                
        self.orbe_a = opts.get(orbe_a,None)
        self.orbe_b = opts.get(orbe_b,None)
        
        return
    
    def set_occs(self,**opts):
        #sets orbital occupation arrays
        self.occs   = opts.get(occs,None)
        self.occs_a = opts.get(occs_a,None)
        self.occs_b = opts.get(occs_b,None)
        
        return
        
    
    def mkdens(self):
        #form density matrix D for restricted wavefuntions
        #for unrestricted returns Da and Db matrices
        if self.restricted:
            #makes spin avg density matrix for open shell systems
            D = self.mk_R_Dmat(0,self.nclosed) #
            if self.nopen>=1:
                D += 0.5*self.mk_R_Dmat(0,self.nclosed+self.nopen)
            return D
        if self.unrestricted:
            Da = self.mk_A_Dmat(0,self.nalpha)
            Db = self.mk_B_Dmat(0,self.nbeta)
            return Da,Db
        
    def mk_R_Dmat(self,nstart,nstop):
        "Form a density matrix C*Ct given eigenvectors C[nstart:nstop,:]"
        d = self.orbs[:,nstart:nstop]
        Dmat = matrixmultiply(d,transpose(d))
        return Dmat

    def mk_A_Dmat(self,nstart,nstop):
        "Form a density matrix C*Ct given eigenvectors C[nstart:nstop,:]"
        d = self.orbs_a[:,nstart:nstop]
        Dmat = matrixmultiply(d,transpose(d))
        return Dmat

    def mk_B_Dmat(self,nstart,nstop):
        "Form a density matrix C*Ct given eigenvectors C[nstart:nstop,:]"
        d = self.orbs_b[:,nstart:nstop]
        Dmat = matrixmultiply(d,transpose(d))
        return Dmat

    def mkQmatrix(self):
        #computes Q matrix which is basically density matrix weighted by the orbital eigenvalues
        if self.restricted:
            occ_orbe = self.orbe[0:self.nclosed]
            occ_orbs = self.orbs[:,0:self.nclosed]
            Qmat = matrixmultiply( occ_orbs, matrixmultiply( diagonal_mat(occ_orbe), transpose(occ_orbs) ) )
            return Qmat
            
        if self.unrestricted:
            occ_orbs_A = self.orbs_a[:,0:self.nalpha]
            occ_orbs_B = self.orbs_b[:,0:self.nbeta]
            occ_orbe_A = self.orbe_a[0:self.nalpha]
            occ_orbe_B = self.orbe_b[0:self.nbeta]
            
            Qa = matrixmultiply( occ_orbs_A, matrixmultiply( diagonal_mat(occ_orbe_A), transpose(occ_orbs_A) ) )
            Qb = matrixmultiply( occ_orbs_B, matrixmultiply( diagonal_mat(occ_orbe_B), transpose(occ_orbs_B) ) )
            return Qa,Qb
            
    def mk_auger_dens(self):
        "Forms a density matrix from a coef matrix c and occupations in occ"
        if self.restricted:
            nstates = self.occs.shape[0]
            D = 0.0
            for i in xrange(nstates):
                D += self.occs[i]*dot( self.orbs[:,i:i+1], transpose(self.orbs[:,i:i+1]))
            #pad_out(D)
            return D
        if self.fixedocc:
            nastates = self.occs_a.shape[0]
            nbstates = self.occs_b.shape[0]
            Da,Db = 0.0,0.0
            for i in xrange(nastates):
                Da += self.occs_a[i]*dot( self.orbs_a[:,i:i+1], transpose(self.orbs_a[:,i:i+1]))
            for i in xrange(nbstates):
                Db += self.occs_b[i]*dot( self.orbs_b[:,i:i+1], transpose(self.orbs_b[:,i:i+1]))
            
            return Da,Db
        
    def mk_auger_Qmatrix(self):
        #computes Q matrix using occ arrays, which is basically density matrix weighted by the orbital eigenvalues
        if self.restricted:
            nstates = self.occs.shape[0]
            Qmat = 0.0
            for i in xrange(nstates):
                Qmat += self.orbe[i]*self.occs[i]*dot( self.orbs[:,i:i+1], transpose(self.orbs[:,i:i+1]))
            return Qmat
            
        if self.fixedocc:
            nastates = self.occs_a.shape[0]
            nbstates = self.occs_b.shape[0]
            Qa,Qb = 0.0,0.0
            for i in xrange(nastates):
                Qa += self.orbe_a[i]*self.occs_a[i]*dot( self.orbs_a[:,i:i+1], transpose(self.orbs_a[:,i:i+1]))
            for i in xrange(nbstates):
                Qb += self.orbe_a[i]*self.occs_b[i]*dot( self.orbs_b[:,i:i+1], transpose(self.orbs_b[:,i:i+1]))
            
            return Qa,Qb
