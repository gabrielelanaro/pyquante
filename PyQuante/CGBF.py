#!/usr/bin/env python
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
from NumWrap import zeros
from math import sqrt

from PyQuante.cints import overlap
from PyQuante.crys import contr_coulomb

class CGBF:
    "Class for a contracted Gaussian basis function"
    def __init__(self,origin,powers=(0,0,0)):
        self._origin = tuple([float(i) for i in origin])
        self._powers = powers
        self._normalization = 1.
        self._prims = []
        self._pnorms = []
        self._pexps = []
        self._pcoefs = []
        return

    def __repr__(self):
        s = "<cgbf origin=\"(%f,%f,%f)\" powers=\"(%d,%d,%d)\">\n" % \
            (self._origin[0],self._origin[1],self._origin[2],
             self._powers[0],self._powers[1],self._powers[2])
        for prim in self._prims:
            s = s + prim.prim_str(self.norm())
        s = s + "</cgbf>\n"
        return s

    def norm(self): return self._normalization
    def origin(self): return self._origin
    def powers(self): return self._powers
    def prims(self): return self._prims

    def exps(self): return self._pexps
    def coefs(self): return self._pcoefs
    def pnorms(self): return self._pnorms

    def center(self,other):
        # Crude estimate to where the center is. The correct form
        #  would use gaussian_product_center, but this has multiple
        #  values for each pair of primitives
        xa,ya,za = self._origin
        xb,yb,zb = other.origin()
        return 0.5*(xa+xb),0.5*(ya+yb),0.5*(za+zb)


    def add_primitive(self,exponent,coefficient):
        "Add a primitive BF to this contracted set"
        pbf = PGBF.PGBF(exponent,self._origin,self._powers)
        pbf._coefficient = coefficient # move into PGBF constructor
        self._prims.append(pbf)
        self._pexps.append(exponent)
        self._pcoefs.append(coefficient)
        return

    def reset_powers(self,px,py,pz):
        self._powers = (px,py,pz)
        for prim in self.prims():
            prim.reset_powers(px,py,pz)
        return

    def normalize(self):
        "Normalize the current CGBF"
        olap = self.overlap(self)
        self._normalization = 1./sqrt(olap)
        for prim in self._prims: self._pnorms.append(prim.norm())

    def overlap(self,other):
        "Overlap matrix element with another CGBF"
        Sij = 0.
        for ipbf in self._prims:
            for jpbf in other._prims:
                Sij = Sij + ipbf.coef()*jpbf.coef()*ipbf.overlap(jpbf)
        return self.norm()*other.norm()*Sij

    def kinetic(self,other):
        "KE matrix element with another CGBF"
        Tij = 0.
        for ipbf in self._prims:
            for jpbf in other._prims:
                Tij = Tij + ipbf.coef()*jpbf.coef()*ipbf.kinetic(jpbf)
        return self.norm()*other.norm()*Tij

    def nuclear(self,other,C):
        "Nuclear matrix element with another CGBF and a center C"
        Vij = 0.
        for ipbf in self._prims:
            for jpbf in other._prims:
                Vij = Vij + ipbf.coef()*jpbf.coef()*ipbf.nuclear(jpbf,C)
        return self.norm()*other.norm()*Vij

    def amp(self,x,y,z):
        "Compute the amplitude of the CGBF at point x,y,z"
        val = 0.
        for prim in self._prims: val+= prim.amp(x,y,z)
        return self._normalization*val

    def move_center(self,dx,dy,dz):
        "Move the basis function to another center"
        self._origin = (self._origin[0]+dx,self._origin[1]+dy,self._origin[2]+dz)
        for prim in self._prims: prim.move_center(dx,dy,dz)
        return

    def doverlap(self,other,dir):
        "Overlap of func with derivative of another"
        dSij = 0.
        l = other.powers()[dir]
        
        ijk_plus = list(other.powers())
        ijk_plus[dir] += 1
        ijk_plus = tuple(ijk_plus)
        
        for ipbf in self._prims:
            for jpbf in other._prims:
                dSij += 2*jpbf.exp()*ipbf.coef()*jpbf.coef()*\
                        ipbf.norm()*jpbf.norm()*\
                        overlap(ipbf.exp(),ipbf.powers(),ipbf.origin(),
                                jpbf.exp(),ijk_plus,jpbf.origin())
        if l>0:
            ijk_minus = list(other.powers())
            ijk_minus[dir] -= 1
            ijk_minus = tuple(ijk_minus)

            for ipbf in self._prims:
                for jpbf in other._prims:
                    dSij -= l*ipbf.coef()*jpbf.coef()*\
                            ipbf.norm()*jpbf.norm()*\
                            overlap(ipbf.exp(),ipbf.powers(),ipbf.origin(),
                                    jpbf.exp(),ijk_minus,jpbf.origin())
        return self.norm()*other.norm()*dSij

    def doverlap_num(self,other,dir):
        "Overlap of func with derivative of another: numeric approximation"
        dSij = 0.
        delta = 0.001 # arbitrary shift amount
        origin_plus = list(other.origin())
        origin_plus[dir] += delta
        origin_plus = tuple(origin_plus)

        origin_minus = list(other.origin())
        origin_minus[dir] -= delta
        origin_minus = tuple(origin_minus)

        for ipbf in self._prims:
            for jpbf in other._prims:
                dSij += 0.5*ipbf.coef()*jpbf.coef()*ipbf.norm()*jpbf.norm()*(
                    overlap(ipbf.exp(),ipbf.powers(),ipbf.origin(),
                            jpbf.exp(),ipbf.powers(),origin_plus)
                    -overlap(ipbf.exp(),ipbf.powers(),ipbf.origin(),
                            jpbf.exp(),ipbf.powers(),origin_plus)
                    )/delta
        return self.norm()*other.norm()*dSij

    def laplacian(self,pos):
        "Evaluate the laplacian of the function at pos=x,y,z"
        val = 0.
        for prim in self._prims: val += prim.laplacian(pos)
        return self._normalization*val

    def grad(self,x,y,z):
        "Evaluate the grad of the function at pos=x,y,z"
        val = zeros(3,'d')
        for prim in self._prims:
            val += prim.grad(x,y,z)
        return self._normalization*val

def coulomb(a,b,c,d):
    "Coulomb interaction between 4 contracted Gaussians"

    #if sum(a.powers()) < sum(b.powers()): a,b = b,a
    #if sum(c.powers()) < sum(d.powers()): c,d = d,c
    
    Jij = contr_coulomb(a.exps(),a.coefs(),a.pnorms(),a.origin(),a.powers(),
                        b.exps(),b.coefs(),b.pnorms(),b.origin(),b.powers(),
                        c.exps(),c.coefs(),c.pnorms(),c.origin(),c.powers(),
                        d.exps(),d.coefs(),d.pnorms(),d.origin(),d.powers())
    return a.norm()*b.norm()*c.norm()*d.norm()*Jij

def three_center(a,b,c):
    import PGBF
    sum = 0
    for ac,ap in zip(a.coefs(),a.prims()):
        for bc,bp in zip(b.coefs(),b.prims()):
            for cc,cp in zip(c.coefs(),c.prims()):
                sum += ac*bc*cc*PGBF.three_center(ap,bp,cp)
    return a.norm()*b.norm()*c.norm()*sum

