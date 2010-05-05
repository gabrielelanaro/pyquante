"""\
 CGBF.py Perform basic operations over contracted gaussian basis
  functions. Uses the functions in PGBF.py.

 References:
  OHT = K. O-ohata, H. Taketa, S. Huzinaga. J. Phys. Soc. Jap. 21, 2306 (1966).
  THO = Taketa, Huzinaga, O-ohata, J. Phys. Soc. Jap. 21,2313 (1966).

 This program is part of the PyQuante quantum chemistry program suite

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

import PGBF
from NumWrap import zeros,array
from math import sqrt

from PyQuante.cints import overlap
#from PyQuante.chgp import contr_coulomb
from PyQuante.crys import contr_coulomb
from PyQuante.contracted_gto import ContractedGTO

class CGBF(ContractedGTO):
    "Class for a contracted Gaussian basis function"
    def __init__(self,origin,powers=(0,0,0),atid=0):
        super(CGBF, self).__init__(origin, powers, atid)
        
        self.origin = tuple([float(i) for i in origin])
        self.powers = powers
        self.norm = 1.
        self.prims = []
        self.pnorms = []
        self.pexps = []
        self.pcoefs = []

        self.ang_mom = sum(powers)
        #added by Hatem H Helal hhh23@cam.ac.uk
        #stores atom id number for easy method to identify which atom
        #a particular bf is centered on
        self.atid = atid
        return

    def __repr__(self):
        s = "<cgbf atomid=%d origin=\"(%f,%f,%f)\" powers=\"(%d,%d,%d)\">\n" % \
            (self.atid,self.origin[0],self.origin[1],self.origin[2],
             self.powers[0],self.powers[1],self.powers[2])
        for prim in self.prims:
            s = s + prim.prim_str(self.norm)
        s = s + "</cgbf>\n"
        return s

    def center(self,other):
        # Crude estimate to where the center is. The correct form
        #  would use gaussian_product_center, but this has multiple
        #  values for each pair of primitives
        xa,ya,za = self.origin
        xb,yb,zb = other.origin
        return 0.5*(xa+xb),0.5*(ya+yb),0.5*(za+zb)


    def add_primitive(self,exponent,coefficient):
        "Add a primitive BF to this contracted set"
        pbf = PGBF.PGBF(exponent,self.origin,self.powers)
        pbf.coef = coefficient # move into PGBF constructor
        super(CGBF,self).add_primitive(pbf,coefficient)
        
        self.pnorms.append(pbf.norm)
        self.prims.append(pbf)
        self.pexps.append(exponent)
        self.pcoefs.append(coefficient)
        return

    def reset_powers(self,px,py,pz):
        self.powers = (px,py,pz)
        for prim in self.prims:
            prim.reset_powers(px,py,pz)
        return

    # Normalize defined in the superclass...
    def overlap(self,other):
        "Overlap matrix element with another CGBF"
        Sij = 0.
        for ipbf in self.prims:
            for jpbf in other.prims:
                Sij = Sij + ipbf.coef*jpbf.coef*ipbf.overlap(jpbf)
        return self.norm*other.norm*Sij

    def kinetic(self,other):
        "KE matrix element with another CGBF"
        Tij = 0.
        for ipbf in self.prims:
            for jpbf in other.prims:
                Tij = Tij + ipbf.coef*jpbf.coef*ipbf.kinetic(jpbf)
        return self.norm*other.norm*Tij

    def multipole(self,other,i,j,k):
        "Overlap matrix element with another CGBF"
        Mij = 0.
        for ipbf in self.prims:
            for jpbf in other.prims:
                Mij += ipbf.coef*jpbf.coef*ipbf.multipole(jpbf,i,j,k)
        return self.norm*other.norm*Mij

    def nuclear(self,other,C):
        "Nuclear matrix element with another CGBF and a center C"
        Vij = 0.
        for ipbf in self.prims:
            for jpbf in other.prims:
                Vij = Vij + ipbf.coef*jpbf.coef*ipbf.nuclear(jpbf,C)
        return self.norm*other.norm*Vij

    def amp(self,x,y,z):
        "Compute the amplitude of the CGBF at point x,y,z"
        val = 0.
        for prim in self.prims: val+= prim.amp(x,y,z)
        return self.norm*val

    def move_center(self,dx,dy,dz):
        "Move the basis function to another center"
        self.origin = (self.origin[0]+dx,self.origin[1]+dy,self.origin[2]+dz)
        for prim in self.prims: prim.move_center(dx,dy,dz)
        return

    def doverlap(self,other,dir):
        "Overlap of func with derivative of another"
        dSij = 0.
        l = other.powers[dir]
        
        ijk_plus = list(other.powers)
        ijk_plus[dir] += 1
        ijk_plus = tuple(ijk_plus)
        
        for ipbf in self.prims:
            for jpbf in other.prims:
                dSij += 2*jpbf.exp*ipbf.coef*jpbf.coef*\
                        ipbf.norm*jpbf.norm*\
                        overlap(ipbf.exp,ipbf.powers,ipbf.origin,
                                jpbf.exp,ijk_plus,jpbf.origin)
        if l>0:
            ijk_minus = list(other.powers)
            ijk_minus[dir] -= 1
            ijk_minus = tuple(ijk_minus)

            for ipbf in self.prims:
                for jpbf in other.prims:
                    dSij -= l*ipbf.coef*jpbf.coef*\
                            ipbf.norm*jpbf.norm*\
                            overlap(ipbf.exp,ipbf.powers,ipbf.origin,
                                    jpbf.exp,ijk_minus,jpbf.origin)
        return self.norm*other.norm*dSij

    def doverlap_num(self,other,dir):
        "Overlap of func with derivative of another: numeric approximation"
        dSij = 0.
        delta = 0.001 # arbitrary shift amount
        origin_plus = list(other.origin)
        origin_plus[dir] += delta
        origin_plus = tuple(origin_plus)

        origin_minus = list(other.origin)
        origin_minus[dir] -= delta
        origin_minus = tuple(origin_minus)

        for ipbf in self.prims:
            for jpbf in other.prims:
                dSij += 0.5*ipbf.coef*jpbf.coef*ipbf.norm*jpbf.norm*(
                    overlap(ipbf.exp,ipbf.powers,ipbf.origin,
                            jpbf.exp,ipbf.powers,origin_plus)
                    -overlap(ipbf.exp,ipbf.powers,ipbf.origin,
                            jpbf.exp,ipbf.powers,origin_plus)
                    )/delta
        return self.norm*other.norm*dSij

    def laplacian(self,pos):
        "Evaluate the laplacian of the function at pos=x,y,z"
        val = 0.
        for prim in self.prims: val += prim.laplacian(pos)
        return self.norm*val

    def grad(self,x,y,z):
        "Evaluate the grad of the function at pos=x,y,z"
        val = zeros(3,'d')
        for prim in self.prims:
            val += prim.grad(x,y,z)
        return self.norm*val

def coulomb(a,b,c,d):
    "Coulomb interaction between 4 contracted Gaussians"

    Jij = contr_coulomb(a.pexps,a.pcoefs,a.pnorms,a.origin,a.powers,
                        b.pexps,b.pcoefs,b.pnorms,b.origin,b.powers,
                        c.pexps,c.pcoefs,c.pnorms,c.origin,c.powers,
                        d.pexps,d.pcoefs,d.pnorms,d.origin,d.powers)
    
    return a.norm*b.norm*c.norm*d.norm*Jij

def three_center(a,b,c):
    import PGBF
    sum = 0
    for ac,ap in zip(a.pcoefs,a.prims):
        for bc,bp in zip(b.pcoefs,b.prims):
            for cc,cp in zip(c.pcoefs,c.prims):
                sum += ac*bc*cc*PGBF.three_center(ap,bp,cp)
    return a.norm*b.norm*c.norm*sum

