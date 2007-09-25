#!/usr/bin/env python
"""\

NAME
      AnalyticDerivatives.py

SYNOPSIS
      Workhorse of the force.py module.


DESCRIPTION
	  
AUTHOR
      Hatem H. Helal, hhh23@cam.ac.uk

REPORT BUGS
      Report bugs to hhh23@cam.ac.uk

COPYRIGHT

"""

from NumWrap import array2string
from math import sqrt
from PyQuante.cints import overlap
from PGBF import PGBF

def der_Hcore_element(a,bfi,bfj,atoms):
    """
    Finds the derivative of the core-Hamiltonian matrix elements, which can
    be written as
    
    H_ij = T_ij + VNe_ij
    
    Where T_ij is the kinetic energy integral and VNe_ij is the nuclear
    attraction integral.
    """
    
    dTij_dXa,dTij_dYa,dTij_dZa = der_kinetic_integral(a,bfi,bfj)
    dVij_dXa,dVij_dYa,dVij_dZa = der_nuc_att(a,bfi,bfj,atoms)
    #if a==1: print "\ndTij_dXa,dTij_dYa,dTij_dZa",dTij_dXa,dTij_dYa,dTij_dZa
    if a==1: print "dVij_dXa,dVij_dYa,dVij_dZa",dVij_dXa,dVij_dYa,dVij_dZa,"\n"
    
    dHij_dXa = dTij_dXa + dVij_dXa
    dHij_dYa = dTij_dYa + dVij_dYa
    dHij_dZa = dTij_dZa + dVij_dZa
    
    return dHij_dXa,dHij_dYa,dHij_dZa 

def der_kinetic_integral(a,bfi,bfj):
    """
    The kinetic energy operator does not depend on the atomic position so we only
    have to consider differentiating the Gaussian functions.  There are 4 possible
    cases we have to evaluate
    
    Case 1: Neither of the basis functions depends on the position of atom A which gives:
        dT_ij/dXa = 0
    
    Cases 2 and 3: Only one of the basis functions depends the position of atom A which 
    gives us either of the following possible integrals to evaluate
        dT_ij/dXa = integral{dr dg_i/dXa T g_j }
        
        dT_ij/dXa = integral{dr g_i T dg_j/dXa  }
    
    Case 4: Both of the basis functions depend on the position of atom A which gives the
    following integral to evaluate
        dT_ij/dXa = integral{dr dg_i/dXa T g_j + g_i T dg_j/dXa  }
    """
    dTij_dXa,dTij_dYa,dTij_dZa = 0.0,0.0,0.0
 
    #we use atom ids on the CGBFs to evaluate which of the 4 above case we have
    #bfi is centered on atom a
    if bfi.atid==a:
        for upbf in bfj.prims():
            for vpbf in bfi.prims():
                alpha = vpbf.exp()
                l,m,n = vpbf.powers()
                origin = vpbf.origin()
                coefs  = upbf.coef()*vpbf.coef()
                                
                #x component
                v = PGBF(alpha,origin,(l+1,m,n))
                v.normalize()
                
                terma = sqrt(alpha*(2.0*l+1.0))*coefs*v.kinetic(upbf)
                
                if l>0:
                    v.reset_powers(l-1,m,n)
                    v.normalize()
                    termb = -2*l*sqrt(alpha/(2.0*l-1.0))*coefs*v.kinetic(upbf)
                else: termb = 0.0
                
                dTij_dXa += terma + termb

                #y component
                v.reset_powers(l,m+1,n)
                v.normalize()
                terma = sqrt(alpha*(2.0*m+1.0))*coefs*v.kinetic(upbf)
                
                if m>0:
                    v.reset_powers(l,m-1,n)
                    v.normalize()
                    termb = -2*m*sqrt(alpha/(2.0*m-1.0))*coefs*v.kinetic(upbf)
                else: termb = 0.0
                
                dTij_dYa += terma + termb
                
                #z component
                v.reset_powers(l,m,n+1)
                v.normalize()
                terma = sqrt(alpha*(2.0*n+1.0))*coefs*v.kinetic(upbf)
                
                if n>0:
                    v.reset_powers(l,m,n-1)
                    v.normalize()
                    termb = -2*n*sqrt(alpha/(2.0*n-1.0))*coefs*v.kinetic(upbf)
                else: termb = 0.0
                
                dTij_dZa += terma + termb
                

    #bfj is centered on atom a
    if bfj.atid==a:
        for upbf in bfi.prims():
            for vpbf in bfj.prims():
                alpha = vpbf.exp()
                l,m,n = vpbf.powers()
                origin = vpbf.origin()
                coefs  = upbf.coef()*vpbf.coef()
                                
                #x component
                v = PGBF(alpha,origin,(l+1,m,n))
                v.normalize()
                
                terma = sqrt(alpha*(2.0*l+1.0))*coefs*v.kinetic(upbf)
                
                if l>0:
                    v.reset_powers(l-1,m,n)
                    v.normalize()
                    termb = -2*l*sqrt(alpha/(2.0*l-1.0))*coefs*v.kinetic(upbf)
                else: termb = 0.0
                
                dTij_dXa += terma + termb
                
                #y component
                v.reset_powers(l,m+1,n)
                v.normalize()
                terma = sqrt(alpha*(2.0*m+1.0))*coefs*v.kinetic(upbf)
                
                if m>0:
                    v.reset_powers(l,m-1,n)
                    v.normalize()
                    termb = -2*m*sqrt(alpha/(2.0*m-1.0))*coefs*v.kinetic(upbf)
                else: termb = 0.0
                
                dTij_dYa += terma + termb
                
                #z component
                v.reset_powers(l,m,n+1)
                v.normalize()
                terma = sqrt(alpha*(2.0*n+1.0))*coefs*v.kinetic(upbf)
                
                if n>0:
                    v.reset_powers(l,m,n-1)
                    v.normalize()
                    termb = -2*n*sqrt(alpha/(2.0*n-1.0))*coefs*v.kinetic(upbf)
                else: termb = 0.0
                
                dTij_dZa += terma + termb

    #print  dTij_dXa,dTij_dYa,dTij_dZa
    return dTij_dXa,dTij_dYa,dTij_dZa
    
def der_nuclear_attraction_integral_test(a,bfi,bfj,atoms):
    """
    The operator describing the attraction between the nuclei and the electrons depends of 
    course on the atomic positions.  This leads to an additional integral involving the 
    derivative of this operator with respect to the atomic position of the form
    
        integral{dr dg_i dV_Ne/dXa g_j } = integral{dr dg_i (x-Xa)/|r-Ra|^3 g_j }
        
    We can show that this integral evaluates to zero by writing it in spherical coordinates.
    So we can proceed by only evaluating the derivatives of the Gaussian basis functions which 
    leads to 4 possible cases.
    
    Case 1: Neither of the basis functions depends on the position of atom A which gives:
        dV_ij/dXa = 0 
    
    Cases 2 and 3: Only one of the basis functions depends on the position of atom A which gives
    us either of the following integrals to evaluate
        dV_ij/dXa = integral{dr dg_i/dXa V g_j }
        
        dV_ij/dXa = integral{dr g_i V dg_j/dXa  }
        
    Case 4: Both basis functions depend on the position of atom A which gives the the following
    integral to evaluate.
        dV_ij/dXa = integral{dr dg_i/dXa V g_j + dr g_i V dg_j/dXa  }
    
    Remember that the nuclear attraction operator is written as a sum over all atoms so this is implemented
    by looping through the atomlist passed in to evaluate the necessary integrals.
    """

    dVij_dXa,dVij_dYa,dVij_dZa = 0.0,0.0,0.0
    #print bfi,bfj
    if bfi.atid==a: #bfi is centered on atom a
        for upbf in bfj.prims():
            for vpbf in bfi.prims():
                alpha = vpbf.exp()
                l,m,n = vpbf.powers()
                origin = vpbf.origin()
                coefs  = upbf.coef()*vpbf.coef()
                                
                #x component
                v = PGBF(alpha,origin,(l+1,m,n))
                v.normalize()
                
                for atom in atoms:
                    dVij_dXa += atom.atno*sqrt(alpha*(2.0*l+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if l>0:
                    v.reset_powers(l-1,m,n)
                    v.normalize()
                    for atom in atoms:
                        dVij_dXa += atom.atno*-2*l*sqrt(alpha/(2.0*l-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                #y component
                v.reset_powers(l,m+1,n)
                v.normalize()
                
                for atom in atoms:
                    dVij_dYa += atom.atno*sqrt(alpha*(2.0*m+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if m>0:
                    v.reset_powers(l,m-1,n)
                    v.normalize()
                    for atom in atoms:
                        dVij_dYa += atom.atno*-2*m*sqrt(alpha/(2.0*m-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                
                #z component
                v.reset_powers(l,m,n+1)
                v.normalize()
                
                for atom in atoms:
                    dVij_dZa += atom.atno*sqrt(alpha*(2.0*n+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if n>0:
                    v.reset_powers(l,m,n-1)
                    v.normalize()
                    for atom in atoms:
                        dVij_dZa += atom.atno*-2*n*sqrt(alpha/(2.0*n-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
        
        
    if bfj.atid==a: #bfj is centered on atom a
        for upbf in bfi.prims():
            for vpbf in bfj.prims():
                alpha = vpbf.exp()
                l,m,n = vpbf.powers()
                origin = vpbf.origin()
                coefs  = upbf.coef()*vpbf.coef()
                                
                #x component
                v = PGBF(alpha,origin,(l+1,m,n))
                v.normalize()
                
                for atom in atoms:
                    dVij_dXa += atom.atno*sqrt(alpha*(2.0*l+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if l>0:
                    v.reset_powers(l-1,m,n)
                    v.normalize()
                    for atom in atoms:
                        dVij_dXa += atom.atno*-2*l*sqrt(alpha/(2.0*l-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                #y component
                v.reset_powers(l,m+1,n)
                v.normalize()
                for atom in atoms:
                    dVij_dYa += atom.atno*sqrt(alpha*(2.0*m+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if m>0:
                    v.reset_powers(l,m-1,n)
                    v.normalize()
                    for atom in atoms:
                        dVij_dYa += atom.atno*-2*m*sqrt(alpha/(2.0*m-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                #z component
                v.reset_powers(l,m,n+1)
                v.normalize()
                for atom in atoms:
                    dVij_dZa += atom.atno*sqrt(alpha*(2.0*n+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if n>0:
                    v.reset_powers(l,m,n-1)
                    v.normalize()
                    for atom in atoms:
                        dVij_dZa += atom.atno*-2*n*sqrt(alpha/(2.0*n-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
    print "\n",bfi.atid,bfj.atid,"\n",dVij_dXa,dVij_dYa,dVij_dZa
    return dVij_dXa,dVij_dYa,dVij_dZa

    
def der_nuclear_attraction_integral(a,bfi,bfj,atoms):
    """
    The operator describing the attraction between the nuclei and the electrons depends of 
    course on the atomic positions.  This leads to an additional integral involving the 
    derivative of this operator with respect to the atomic position of the form
    
        integral{dr dg_i dV_Ne/dXa g_j } = integral{dr dg_i (x-Xa)/|r-Ra|^3 g_j }
        
    We can show that this integral evaluates to zero by writing it in spherical coordinates.
    So we can proceed by only evaluating the derivatives of the Gaussian basis functions which 
    leads to 4 possible cases.
    
    Case 1: Neither of the basis functions depends on the position of atom A which gives:
        dV_ij/dXa = 0 
    
    Cases 2 and 3: Only one of the basis functions depends on the position of atom A which gives
    us either of the following integrals to evaluate
        dV_ij/dXa = integral{dr dg_i/dXa V g_j }
        
        dV_ij/dXa = integral{dr g_i V dg_j/dXa  }
        
    Case 4: Both basis functions depend on the position of atom A which gives the the following
    integral to evaluate.
        dV_ij/dXa = integral{dr dg_i/dXa V g_j + dr g_i V dg_j/dXa  }
    
    Remember that the nuclear attraction operator is written as a sum over all atoms so this is implemented
    by looping through the atomlist passed in to evaluate the necessary integrals.
    """

    dVij_dXa,dVij_dYa,dVij_dZa = 0.0,0.0,0.0

    for atom in atoms: #we loop over all atoms to evaluate the attraction integrals
        if a==1: print atom,"\n",dVij_dXa,dVij_dYa,dVij_dZa
        if bfi.atid==a: #bfi is centered on atom a
            for upbf in bfj.prims():
                for vpbf in bfi.prims():
                    alpha = vpbf.exp()
                    l,m,n = vpbf.powers()
                    origin = vpbf.origin()
                    coefs  = upbf.coef()*vpbf.coef()
                                    
                    #x component
                    v = PGBF(alpha,origin,(l+1,m,n))
                    v.normalize()
                    
                    terma = sqrt(alpha*(2.0*l+1.0))*coefs*v.nuclear(upbf,atom.pos())
                    
                    if l>0:
                        v.reset_powers(l-1,m,n)
                        v.normalize()
                        termb = -2*l*sqrt(alpha/(2.0*l-1.0))*coefs*v.nuclear(upbf,atom.pos())
                    else: termb = 0.0
                    
                    dVij_dXa += atom.atno*(terma + termb)
                    
                    #y component
                    v.reset_powers(l,m+1,n)
                    v.normalize()
                    terma = sqrt(alpha*(2.0*m+1.0))*coefs*v.nuclear(upbf,atom.pos())
                    
                    if m>0:
                        v.reset_powers(l,m-1,n)
                        v.normalize()
                        termb = -2*m*sqrt(alpha/(2.0*m-1.0))*coefs*v.nuclear(upbf,atom.pos())
                    else: termb = 0.0
                    
                    dVij_dYa += atom.atno*(terma + termb)
                    
                    #z component
                    v.reset_powers(l,m,n+1)
                    v.normalize()
                    terma = sqrt(alpha*(2.0*n+1.0))*coefs*v.nuclear(upbf,atom.pos())
                    
                    if n>0:
                        v.reset_powers(l,m,n-1)
                        v.normalize()
                        termb = -2*n*sqrt(alpha/(2.0*n-1.0))*coefs*v.nuclear(upbf,atom.pos())
                    else: termb = 0.0
                    
                    dVij_dZa += atom.atno*(terma + termb)
                    
                    
            
        if bfj.atid==a: #bfj is centered on atom a
            for upbf in bfi.prims():
                for vpbf in bfj.prims():
                    alpha = vpbf.exp()
                    l,m,n = vpbf.powers()
                    origin = vpbf.origin()
                    coefs  = upbf.coef()*vpbf.coef()
                                    
                    #x component
                    v = PGBF(alpha,origin,(l+1,m,n))
                    v.normalize()
                    
                    terma = sqrt(alpha*(2.0*l+1.0))*coefs*v.nuclear(upbf,atom.pos())
                    
                    if l>0:
                        v.reset_powers(l-1,m,n)
                        v.normalize()
                        termb = -2*l*sqrt(alpha/(2.0*l-1.0))*coefs*v.nuclear(upbf,atom.pos())
                    else: termb = 0.0
                    
                    dVij_dXa += atom.atno*(terma + termb)
                    
                    #y component
                    v.reset_powers(l,m+1,n)
                    v.normalize()
                    terma = sqrt(alpha*(2.0*m+1.0))*coefs*v.nuclear(upbf,atom.pos())
                    
                    if m>0:
                        v.reset_powers(l,m-1,n)
                        v.normalize()
                        termb = -2*m*sqrt(alpha/(2.0*m-1.0))*coefs*v.nuclear(upbf,atom.pos())
                    else: termb = 0.0
                    
                    dVij_dYa += atom.atno*(terma + termb)
                    
                    #z component
                    v.reset_powers(l,m,n+1)
                    v.normalize()
                    terma = sqrt(alpha*(2.0*n+1.0))*coefs*v.nuclear(upbf,atom.pos())
                    
                    if n>0:
                        v.reset_powers(l,m,n-1)
                        v.normalize()
                        termb = -2*n*sqrt(alpha/(2.0*n-1.0))*coefs*v.nuclear(upbf,atom.pos())
                    else: termb = 0.0
                    
                    dVij_dZa += atom.atno*(terma + termb)
                    
    return dVij_dXa,dVij_dYa,dVij_dZa

def der_overlap_element(a,bfi, bfj):
    """
    finds the derivative of the overlap integral with respect to the 
    atomic coordinate of atom "a".  Note there are four possible cases
    for evaluating this integral:
     1. Neither of the basis functions depend on the position of atom a
        ie. they are centered on atoms other than atom a
     2 and 3. One of the basis functions depends on the position of atom a
        so we need to evaluate the derivative of a Gaussian with the 
        recursion (right word?) relation derived on page 442 of Szabo.
     4. Both of the basis functions are centered on atom a, which through the
        recursion relation for the derivative of a Gaussian basis function will
        require the evaluation of 4 overlap integrals...

    this function will return a 3 element list with the derivatives of the overlap
    integrals with respect to the atomic coordinates Xa,Ya,Za.
    """
    dSij_dXa,dSij_dYa,dSij_dZa = 0.0,0.0,0.0
    
    #we use atom ids on the CGBFs to evaluate which of the 4 above case we have
    if bfi.atid==a: #bfi is centered on atom a
        for upbf in bfj.prims():
            for vpbf in bfi.prims():
                alpha = vpbf.exp()
                l,m,n = vpbf.powers()
                origin = vpbf.origin()
                coefs  = upbf.coef()*vpbf.coef()
                                
                #x component
                v = PGBF(alpha,origin,(l+1,m,n))
                v.normalize()
                
                terma = sqrt(alpha*(2.0*l+1.0))*coefs*v.overlap(upbf)
                
                if l>0:
                    v.reset_powers(l-1,m,n)
                    v.normalize()
                    termb = -2*l*sqrt(alpha/(2.0*l-1.0))*coefs*v.overlap(upbf)
                else: termb = 0.0
                
                dSij_dXa += terma + termb
                
                
                #y component
                v.reset_powers(l,m+1,n)
                v.normalize()
                terma = sqrt(alpha*(2.0*m+1.0))*coefs*v.overlap(upbf)
                
                if m>0:
                    v.reset_powers(l,m-1,n)
                    v.normalize()
                    termb = -2*m*sqrt(alpha/(2.0*m-1.0))*coefs*v.overlap(upbf)
                else: termb = 0.0
                
                dSij_dYa += terma + termb
                
                #z component
                v.reset_powers(l,m,n+1)
                v.normalize()
                terma = sqrt(alpha*(2.0*n+1.0))*coefs*v.overlap(upbf)
                
                if n>0:
                    v.reset_powers(l,m,n-1)
                    v.normalize()
                    termb = -2*n*sqrt(alpha/(2.0*n-1.0))*coefs*v.overlap(upbf)
                else: termb = 0.0
                
                dSij_dZa += terma + termb

                
    #bfj is centered on atom a
    if bfj.atid==a:
        for upbf in bfi.prims():
            for vpbf in bfj.prims():
                alpha = vpbf.exp()
                l,m,n = vpbf.powers()
                origin = vpbf.origin()
                coefs = upbf.coef()*vpbf.coef()
                
                #x component 
                v = PGBF(alpha,origin,(l+1,m,n))
                v.normalize()
                
                terma = sqrt(alpha*(2.0*l+1.0))*coefs*v.overlap(upbf)
                
                if l>0:
                    v.reset_powers(l-1,m,n)
                    v.normalize()
                    termb = -2*l*sqrt(alpha/(2.0*l-1.0))*coefs*v.overlap(upbf)
                else: termb = 0.0
                
                dSij_dXa += terma + termb
                

                #y component
                v.reset_powers(l,m+1,n)
                v.normalize()
                terma = sqrt(alpha*(2.0*m+1.0))*coefs*v.overlap(upbf)
                
                if m>0:
                    v.reset_powers(l,m-1,n)
                    v.normalize()
                    termb = -2*m*sqrt(alpha/(2.0*m-1.0))*coefs*v.overlap(upbf)
                else: termb = 0.0
                
                dSij_dYa += terma + termb
                
                #z component
                v.reset_powers(l,m,n+1)
                v.normalize()
                terma = sqrt(alpha*(2.0*n+1.0))*coefs*v.overlap(upbf)
                
                if n>0:
                    v.reset_powers(l,m,n-1)
                    v.normalize()
                    termb = -2*n*sqrt(alpha/(2.0*n-1.0))*coefs*v.overlap(upbf)
                else: termb = 0.0
                
                dSij_dZa += terma + termb

    return dSij_dXa,dSij_dYa,dSij_dZa

def der_nuc_att(a,bfi,bfj,atoms):
    dVij_dXa,dVij_dYa,dVij_dZa = 0.0,0.0,0.0
    #we use atom ids on the CGBFs to evaluate which of the 4 above case we have
    if bfi.atid==a: #bfi is centered on atom a
        for upbf in bfj.prims():
            for vpbf in bfi.prims():
                alpha = vpbf.exp()
                l,m,n = vpbf.powers()
                origin = vpbf.origin()
                coefs  = upbf.coef()*vpbf.coef()
                                
                #x component
                v = PGBF(alpha,origin,(l+1,m,n))
                v.normalize()
                
                terma=0.0
                for atom in atoms:                    
                    terma += atom.atno*sqrt(alpha*(2.0*l+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if l>0:
                    v.reset_powers(l-1,m,n)
                    v.normalize()
                    
                    termb=0.0
                    for atom in atoms:
                        termb += -2*l*atom.atno*sqrt(alpha/(2.0*l-1.0))*coefs*v.nuclear(upbf,atom.pos())
                        
                else: termb = 0.0
                
                dVij_dXa += terma + termb
                
                
                #y component
                v.reset_powers(l,m+1,n)
                v.normalize()

                terma=0.0
                for atom in atoms:
                    terma += atom.atno*sqrt(alpha*(2.0*m+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if m>0:
                    v.reset_powers(l,m-1,n)
                    v.normalize()
                    
                    termb=0.0
                    for atom in atoms:
                        termb += -2*m*atom.atno*sqrt(alpha/(2.0*m-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                dVij_dYa += terma + termb
                
                #z component
                v.reset_powers(l,m,n+1)
                v.normalize()
                
                terma=0.0
                for atom in atoms:
                    terma += atom.atno*sqrt(alpha*(2.0*n+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if n>0:
                    v.reset_powers(l,m,n-1)
                    v.normalize()
                    
                    termb=0.0
                    for atom in atoms:
                        termb += -2*n*atom.atno*sqrt(alpha/(2.0*n-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                dVij_dZa += terma + termb

                
    #bfj is centered on atom a
    if bfj.atid==a:
        for upbf in bfi.prims():
            for vpbf in bfj.prims():
                alpha = vpbf.exp()
                l,m,n = vpbf.powers()
                origin = vpbf.origin()
                coefs = upbf.coef()*vpbf.coef()
                
                #x component 
                v = PGBF(alpha,origin,(l+1,m,n))
                v.normalize()
                
                terma=0.0
                for atom in atoms:
                    terma += atom.atno*sqrt(alpha*(2.0*l+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if l>0:
                    v.reset_powers(l-1,m,n)
                    v.normalize()
                    
                    termb=0.0
                    for atom in atoms:
                        termb += -2*l*atom.atno*sqrt(alpha/(2.0*l-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                dVij_dXa += terma + termb
                

                #y component
                v.reset_powers(l,m+1,n)
                v.normalize()

                terma=0.0
                for atom in atoms:
                    terma += atom.atno*sqrt(alpha*(2.0*m+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if m>0:
                    v.reset_powers(l,m-1,n)
                    v.normalize()
                    
                    termb=0.0
                    for atom in atoms:
                        termb += -2*m*atom.atno*sqrt(alpha/(2.0*m-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                dVij_dYa += terma + termb
                
                #z component
                v.reset_powers(l,m,n+1)
                v.normalize()
                
                terma=0.0
                for atom in atoms:
                    terma += atom.atno*sqrt(alpha*(2.0*n+1.0))*coefs*v.nuclear(upbf,atom.pos())
                
                if n>0:
                    v.reset_powers(l,m,n-1)
                    v.normalize()
                    
                    termb = 0.0
                    for atom in atoms:
                        termb += -2*n*atom.atno*sqrt(alpha/(2.0*n-1.0))*coefs*v.nuclear(upbf,atom.pos())
                else: termb = 0.0
                
                dVij_dZa += terma + termb

    return dVij_dXa,dVij_dYa,dVij_dZa
    
def num_der_nuc_att(a,bfi,bfj,atoms):
    dVij_dXa,dVij_dYa,dVij_dZa = 0.0,0.0,0.0
    
    delta = 0.001
    x0,y0,z0 = a.r
    
    for atom in atoms:
        if atom.atid == a:
            bfi.nuclear(bfj,atom.pos())
    
    a.update_coords((x0+delta,y0,z0))

    
    return dVij_dXa,dVij_dYa,dVij_dZa

    
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
