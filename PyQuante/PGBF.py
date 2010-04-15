"""\
  PGBF.py Perform basic operations over primitive
    gaussian basis functions. The equations herein are based upon
    'Gaussian Expansion Methods for Molecular Orbitals.' H. Taketa,
    S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
    [THO paper].

  For the purposes of this routine, a gaussian is defined as:

    g(x,y,z) = A*(x^i)*(y^j)*(z^k)*exp{-a*(r-ro)^2}

  The functions defined are:

  overlap(g'): Compute the overlap matrix element of g with g': Int(g*g')
  
  kinetic(g'): Compute the kinetic energy matrix element
    between g and g' = Int(G*lapl(G')), where lapl is the Laplacian.

  nuclear(g',r): Compute the nuclear attraction integral
    Int(g*(1/r)*g'). Only programmed for 1s gaussians.

  coulomb(g,g',g'',g'''): Compute the two-electron colombic repulsion
    integral Int(g(1)g'(1)(1/r12)g''(2)g'''(2)). 

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from math import sqrt,pi,pow,exp
from NumWrap import array

from PyQuante.cints import kinetic,overlap,nuclear_attraction,fact2,dist2
from PyQuante.cints import binomial, three_center_1D
from PyQuante.chgp import coulomb_repulsion

#added 2/8/07 by Hatem Helal hhh23@cam.ac.uk
#probably need to write the C version in cints...
from PyQuante.pyints import grad_nuc_att

class PGBF:
    "Class for Primitive Gaussian Basis Functions."

    # Constructor
    def __init__(self,exponent,origin,powers=(0,0,0),norm=1.):
        self._exponent = float(exponent)
        self._origin = tuple([float(i) for i in origin])
        self._powers= powers
        self._normalization = float(norm)
        self.normalize()
        self._coefficient = 1
        return

    # Public
    def exp(self): return self._exponent
    def origin(self): return self._origin
    def powers(self): return self._powers
    def norm(self): return self._normalization
    def coef(self): return self._coefficient

    def reset_powers(self,px,py,pz):
        self._powers = (px,py,pz)
        return

    def overlap(self,other):
        "Compute overlap element with another PGBF"
        #print self._normalization,other._normalization,\
        #       self._exponent,self._powers,self._origin,\
        #       other._exponent,other._powers,other._origin
        return self._normalization*other._normalization*\
               overlap(self._exponent,self._powers,self._origin,
                       other._exponent,other._powers,other._origin)

    def kinetic(self,other):
        "Overlap between two gaussians. THO eq. 2.14."
        return self._normalization*other._normalization*\
               kinetic(self._exponent,self._powers,self._origin,
                       other._exponent,other._powers,other._origin)

    def multipole(self,other,i,j,k):
        from pyints import multipole_ints
        return self._normalization*other._normalization*\
               multipole_ints((i,j,k),
                              self._exponent,self._powers,self._origin,
                              other._exponent,other._powers,other._origin)


    #Need to rewrite this to:
    #  1. pull the _normalization constants out front to be consistent
    #     with overlap() and kinetic()
    #  2. reorder the arguments to be in the same order as overlap()
    #     and kinetic()
    def nuclear(self,other,C):
        "THO eq. 2.17 and 3.1"
        return nuclear_attraction(self._origin,self._normalization,
                                  self._powers,self._exponent,
                                  other._origin,other._normalization,
                                  other._powers,other._exponent,
                                  C)

    def nuclear_gradient(self,other,C):
        return self._normalization*other._normalization*\
               array(grad_nuc_att(self._origin,self._powers,self._exponent,
                            other._origin,other._powers,other._exponent,
                            C))

    def amp(self,x,y,z):
        "Compute the amplitude of the PGBF at point x,y,z"
        i,j,k = self._powers
        x0,y0,z0 = self._origin
        return self._normalization*self._coefficient*\
               pow(x-x0,i)*pow(y-y0,j)*pow(z-z0,k)*\
               exp(-self._exponent*dist2((x,y,z),(x0,y0,z0)))

    def move_center(self,dx,dy,dz):
        self._origin = (self._origin[0]+dx,self._origin[1]+dy,self._origin[2]+dz)
        return

    # Private
    def normalize(self):
        "Normalize basis function. From THO eq. 2.2"
        l,m,n = self._powers
        alpha = self._exponent
        self._normalization = sqrt(pow(2,2*(l+m+n)+1.5)*
                                   pow(alpha,l+m+n+1.5)/
                                   fact2(2*l-1)/fact2(2*m-1)/
                                   fact2(2*n-1)/pow(pi,1.5))
        # This is a hack to allow zero-valued exponents, for testing
        if abs(alpha) < 1e-8: self._normalization = 1.
        return


    # Other overloads
    def __str__(self):
	return "PGBF(%.2f," % self._exponent +\
               "(%.2f,%.2f,%.2f)," % self._origin +\
               "(%d,%d,%d)," % self._powers +\
               "%.2f)" % self._normalization

    def prim_str(self,topnorm=1):
        return "    <prim exp=\"%6.4f\" coeff=\"%6.4f\" ncoeff=\"%6.4f\"/>\n" \
               % (self.exp(),self.coef(),topnorm*self.norm()*self.coef())

    def laplacian(self,pos):
        amp = self.amp(pos[0],pos[1],pos[2])
        alpha = self._exponent
        x = pos[0]-self._origin[0]
        y = pos[1]-self._origin[1]
        z = pos[2]-self._origin[2]
        x2 = x*x
        y2 = y*y
        z2 = z*z
        r2 = x2+y2+z2
        L,M,N = self._powers
        term = (L*(L-1)/x2 + M*(M-1)/y2 + N*(N-1)/z2) +\
                4*alpha*alpha*r2 - 2*alpha*(2*(L+M+N)+3)
        return self._normalization*self._coefficient*amp*term

    def grad_old(self,pos):
        amp = self.amp(pos[0],pos[1],pos[2])
        alpha = self._exponent
        L,M,N = self._powers
        x = pos[0]-self._origin[0]
        y = pos[1]-self._origin[1]
        z = pos[2]-self._origin[2]
        val = array([L/x - 2*x*alpha,M/y - 2*y*alpha,N/z-2*z*alpha])
        return self._normalization*self._coefficient*val*amp

    def grad(self,x,y,z):
        alpha = self._exponent
        I,J,K = self._powers
        C = self._normalization*self._coefficient
        x0,y0,z0 = self._origin
        fx = pow(x-x0,I)*exp(-alpha*pow(x-x0,2))
        fy = pow(y-y0,J)*exp(-alpha*pow(y-y0,2))
        fz = pow(z-z0,K)*exp(-alpha*pow(z-z0,2))
        gx = -2*alpha*(x-x0)*fx
        gy = -2*alpha*(y-y0)*fy
        gz = -2*alpha*(z-z0)*fz
        if I > 0: gx += pow(x-x0,I-1)*exp(-alpha*pow(x-x0,2))
        if J > 0: gy += pow(y-y0,J-1)*exp(-alpha*pow(y-y0,2))
        if K > 0: gz += pow(z-z0,K-1)*exp(-alpha*pow(z-z0,2))
        return array([C*gx*fy*fz,C*fx*gy*fz,C*fx*fy*gz])
        

# Friend functions
def coulomb(gA,gB,gC,gD):
    """Coulomb interaction between four cartesian Gaussians; THO eq. 2.22"""
    return coulomb_repulsion(gA._origin,gA._normalization,gA._powers,
                             gA._exponent,gB._origin,gB._normalization,
                             gB._powers,gB._exponent,gC._origin,
                             gC._normalization,gC._powers,gC._exponent,
                             gD._origin,gD._normalization,gD._powers,
                             gD._exponent)

def three_center(gA,gB,gC):
    "Three-center integral between Gaussians"
    na = gA._normalization
    nb = gB._normalization
    nc = gC._normalization
    ix = three_center_1D(gA._origin[0],gA._powers[0],gA._exponent,
                         gB._origin[0],gB._powers[0],gB._exponent,
                         gC._origin[0],gC._powers[0],gC._exponent)
    iy = three_center_1D(gA._origin[1],gA._powers[1],gA._exponent,
                        gB._origin[1],gB._powers[1],gB._exponent,
                        gC._origin[1],gC._powers[1],gC._exponent)
    iz = three_center_1D(gA._origin[2],gA._powers[2],gA._exponent,
                        gB._origin[2],gB._powers[2],gB._exponent,
                        gC._origin[2],gC._powers[2],gC._exponent)
    return na*nb*nc*ix*iy*iz

def test_3cent():
    gA = PGBF(1,(0,0,0))
    gB = PGBF(1,(1,0,0))
    gC = PGBF(2.,(0,0,0))
    gD = PGBF(1,(0,0,0),(1,0,0))
    # Here we construct a "fake" function that has zero exponent, so that
    #  the three-center integrals that follow will be the same as a
    #  normal overlap integral
    g0 = PGBF(0,(0,0,0))
    for a,b in [(gA,gA),(gB,gB),(gA,gB),(gC,gC),(gA,gC),
                (gD,gD),(gA,gD),(gB,gD)]: tester(a,b)

def tester(gA,gB):
    # insure that the overlap integrals <ga|gb> and <gb|ga> are equal
    # to the various three center integrals <ga|gb|g0> and its
    # various permutations
    g0 = PGBF(0,(0,0,0))
    olab = gA.overlap(gB)
    olba = gB.overlap(gA)
    tcab0 = three_center(gA,gB,g0)
    tcba0 = three_center(gB,gA,g0)
    tca0b = three_center(gA,g0,gB)
    tcb0a = three_center(gB,g0,gA)
    tc0ab = three_center(g0,gA,gB)
    tc0ba = three_center(g0,gB,gA)
    diff = 0
    vals = [olab,olba,tcab0,tcba0,tca0b,tcb0a,tc0ab,tc0ba]
    for i in xrange(len(vals)):
        for j in xrange(i):
            diff = max(diff,abs(vals[i]-vals[j]))
    print "For ints %s %s max diff is %f" % (gA,gB,diff)

if __name__ == '__main__':  test_3cent()
