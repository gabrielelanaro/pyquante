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
from PGBF import PGBF,coulomb
from pyints import grad_nuc_att

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
        for upbf in bfj.prims:
            for vpbf in bfi.prims:
                alpha = vpbf.exp
                l,m,n = vpbf.powers
                origin = vpbf.origin
                coefs  = upbf.coef*vpbf.coef
                                
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
        for upbf in bfi.prims:
            for vpbf in bfj.prims:
                alpha = vpbf.exp
                l,m,n = vpbf.powers
                origin = vpbf.origin
                coefs  = upbf.coef*vpbf.coef
                                
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

    return dTij_dXa,dTij_dYa,dTij_dZa

def der_nuc_att(a,bfi,bfj,atoms):
    """
    This function finds the atomic gradient of the nuclear attraction integrals. Since the
    nuclear attraction operator explicitly depends on the atomic coordinates we find 
    
    grad <i|V|j> = <grad i|V|j> + <i|V|grad j> + <i| grad V |j>
    
    The first two terms are straightforward to evaluate using the recursion relation for the
    derivative of a Gaussian basis function.  The last term found through the nuclear_gradient
    function in the primitive Gaussian class.  
    """
    
    dVij_dXa,dVij_dYa,dVij_dZa = 0.0,0.0,0.0
    if bfi.atid==a: #bfi is centered on atom a
        for upbf in bfj.prims:
            for vpbf in bfi.prims:
                alpha = vpbf.exp
                l,m,n = vpbf.powers
                origin = vpbf.origin
                coefs  = upbf.coef*vpbf.coef
                                
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
        for upbf in bfi.prims:
            for vpbf in bfj.prims:
                alpha = vpbf.exp
                l,m,n = vpbf.powers
                origin = vpbf.origin
                coefs = upbf.coef*vpbf.coef
                
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
    
    #finally evaluate <i| grad V |j>
    for atom in atoms:
        if atom.atid==a:
            for upbf in bfi.prims:
                for vpbf in bfj.prims:
                    prefactor = upbf.coef*vpbf.coef*atom.atno
                    list = upbf.nuclear_gradient(vpbf,atom.pos())
                    
                    dVij_dXa+=prefactor*list[0]
                    dVij_dYa+=prefactor*list[1]
                    dVij_dZa+=prefactor*list[2]
        

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
        for upbf in bfj.prims:
            for vpbf in bfi.prims:
                alpha = vpbf.exp
                l,m,n = vpbf.powers
                origin = vpbf.origin
                coefs  = upbf.coef*vpbf.coef
                                
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
        for upbf in bfi.prims:
            for vpbf in bfj.prims:
                alpha = vpbf.exp
                l,m,n = vpbf.powers
                origin = vpbf.origin
                coefs = upbf.coef*vpbf.coef
                
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
    
def der_Jints(a, bfi,bfj,bfk,bfl):
    """
    This function will find the atomic gradient of the Coloumb integral over
    basis functions i,j,k, and l as in
    
    grad_a <ij|kl> = <gi j|kl> + <i gj|kl> + <ij|gk l> + <ij|k gl>
    """
    dJint_dXa,dJint_dYa,dJint_dZa = 0.0,0.0,0.0
    
    if bfi.atid==a: #bfi is centered on atom a
        for tpbf in bfi.prims:
            for upbf in bfj.prims:
                for vpbf in bfk.prims:
                    for wpbf in bfl.prims:
                        alpha  = tpbf.exp
                        l,m,n  = tpbf.powers
                        origin = tpbf.origin
                        coefs  = tpbf.coef*upbf.coef*vpbf.coef*wpbf.coef
                        
                        #x component
                        tmp = PGBF(alpha, origin,(l+1,m,n)) #temp pgbf
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*l+1.0))*coefs*coulomb(tmp,upbf,vpbf,wpbf)
                        
                        if l>0:
                            tmp.reset_powers(l-1,m,n)
                            tmp.normalize()
                            termb = -2*l*sqrt(alpha/(2.*l-1))*coefs*coulomb(tmp,upbf,vpbf,wpbf)
                        else: termb = 0.0
                        
                        dJint_dXa += terma+termb
                        
                        #y component
                        tmp.reset_powers(l,m+1,n)
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*m+1.0))*coefs*coulomb(tmp,upbf,vpbf,wpbf)
                        
                        if m>0:
                            tmp.reset_powers(l,m-1,n)
                            tmp.normalize()
                            termb = -2*m*sqrt(alpha/(2.*m-1))*coefs*coulomb(tmp,upbf,vpbf,wpbf)
                        else: termb=0.0
                        
                        dJint_dYa += terma + termb

                        #z component
                        tmp.reset_powers(l,m,n+1)
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*n+1.0))*coefs*coulomb(tmp,upbf,vpbf,wpbf)
                        
                        if n>0:
                            tmp.reset_powers(l,m,n-1)
                            tmp.normalize()
                            termb = -2*n*sqrt(alpha/(2.*n-1))*coefs*coulomb(tmp,upbf,vpbf,wpbf)
                        else: termb=0.0
                        
                        dJint_dZa += terma + termb
                        
    if bfj.atid==a: #bfj is centered on atom a
        for tpbf in bfi.prims:
            for upbf in bfj.prims:
                for vpbf in bfk.prims:
                    for wpbf in bfl.prims:
                        alpha  = upbf.exp
                        l,m,n  = upbf.powers
                        origin = upbf.origin
                        coefs  = tpbf.coef*upbf.coef*vpbf.coef*wpbf.coef
                        
                        #x component
                        tmp = PGBF(alpha, origin,(l+1,m,n)) #temp pgbf
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*l+1.0))*coefs*coulomb(tpbf,tmp,vpbf,wpbf)
                        
                        if l>0:
                            tmp.reset_powers(l-1,m,n)
                            tmp.normalize()
                            termb = -2*l*sqrt(alpha/(2.*l-1))*coefs*coulomb(tpbf,tmp,vpbf,wpbf)
                        else: termb = 0.0
                        
                        dJint_dXa += terma+termb
                        
                        #y component
                        tmp.reset_powers(l,m+1,n)
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*m+1.0))*coefs*coulomb(tpbf,tmp,vpbf,wpbf)
                        
                        if m>0:
                            tmp.reset_powers(l,m-1,n)
                            tmp.normalize()
                            termb = -2*m*sqrt(alpha/(2.*m-1))*coefs*coulomb(tpbf,tmp,vpbf,wpbf)
                        else: termb=0.0
                        
                        dJint_dYa += terma + termb

                        #z component
                        tmp.reset_powers(l,m,n+1)
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*n+1.0))*coefs*coulomb(tpbf,tmp,vpbf,wpbf)
                        
                        if n>0:
                            tmp.reset_powers(l,m,n-1)
                            tmp.normalize()
                            termb = -2*n*sqrt(alpha/(2.*n-1))*coefs*coulomb(tpbf,tmp,vpbf,wpbf)
                        else: termb=0.0
                        
                        dJint_dZa += terma + termb
    
    if bfk.atid==a: #bfk is centered on atom a
        for tpbf in bfi.prims:
            for upbf in bfj.prims:
                for vpbf in bfk.prims:
                    for wpbf in bfl.prims:
                        alpha  = vpbf.exp
                        l,m,n  = vpbf.powers
                        origin = vpbf.origin
                        coefs  = tpbf.coef*upbf.coef*vpbf.coef*wpbf.coef
                        
                        #x component
                        tmp = PGBF(alpha, origin,(l+1,m,n)) #temp pgbf
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*l+1.0))*coefs*coulomb(tpbf,upbf,tmp,wpbf)
                        
                        if l>0:
                            tmp.reset_powers(l-1,m,n)
                            tmp.normalize()
                            termb = -2*l*sqrt(alpha/(2.*l-1))*coefs*coulomb(tpbf,upbf,tmp,wpbf)
                        else: termb = 0.0
                        
                        dJint_dXa += terma+termb
                        
                        #y component
                        tmp.reset_powers(l,m+1,n)
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*m+1.0))*coefs*coulomb(tpbf,upbf,tmp,wpbf)
                        
                        if m>0:
                            tmp.reset_powers(l,m-1,n)
                            tmp.normalize()
                            termb = -2*m*sqrt(alpha/(2.*m-1))*coefs*coulomb(tpbf,upbf,tmp,wpbf)
                        else: termb=0.0
                        
                        dJint_dYa += terma + termb

                        #z component
                        tmp.reset_powers(l,m,n+1)
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*n+1.0))*coefs*coulomb(tpbf,upbf,tmp,wpbf)
                        
                        if n>0:
                            tmp.reset_powers(l,m,n-1)
                            tmp.normalize()
                            termb = -2*n*sqrt(alpha/(2.*n-1))*coefs*coulomb(tpbf,upbf,tmp,wpbf)
                        else: termb=0.0
                        
                        dJint_dZa += terma + termb
    
    if bfl.atid==a: #bfl is centered on atom a
        for tpbf in bfi.prims:
            for upbf in bfj.prims:
                for vpbf in bfk.prims:
                    for wpbf in bfl.prims:
                        alpha  = wpbf.exp
                        l,m,n  = wpbf.powers
                        origin = wpbf.origin
                        coefs  = tpbf.coef*upbf.coef*vpbf.coef*wpbf.coef
                        
                        #x component
                        tmp = PGBF(alpha, origin,(l+1,m,n)) #temp pgbf
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*l+1.0))*coefs*coulomb(tpbf,upbf,vpbf,tmp)
                        
                        if l>0:
                            tmp.reset_powers(l-1,m,n)
                            tmp.normalize()
                            termb = -2*l*sqrt(alpha/(2.*l-1))*coefs*coulomb(tpbf,upbf,vpbf,tmp)
                        else: termb = 0.0
                        
                        dJint_dXa += terma+termb
                        
                        #y component
                        tmp.reset_powers(l,m+1,n)
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*m+1.0))*coefs*coulomb(tpbf,upbf,vpbf,tmp)
                        
                        if m>0:
                            tmp.reset_powers(l,m-1,n)
                            tmp.normalize()
                            termb = -2*m*sqrt(alpha/(2.*m-1))*coefs*coulomb(tpbf,upbf,vpbf,tmp)
                        else: termb=0.0
                        
                        dJint_dYa += terma + termb

                        #z component
                        tmp.reset_powers(l,m,n+1)
                        tmp.normalize()
                        
                        terma = sqrt(alpha*(2.0*n+1.0))*coefs*coulomb(tpbf,upbf,vpbf,tmp)
                        
                        if n>0:
                            tmp.reset_powers(l,m,n-1)
                            tmp.normalize()
                            termb = -2*n*sqrt(alpha/(2.*n-1))*coefs*coulomb(tpbf,upbf,vpbf,tmp)
                        else: termb=0.0
                        
                        dJint_dZa += terma + termb
    
    return dJint_dXa,dJint_dYa,dJint_dZa
    
