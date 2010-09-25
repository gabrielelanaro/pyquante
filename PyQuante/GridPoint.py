"""\
 GridPoint.py: A class to hold grid point data

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from math import sqrt
from NumWrap import zeros,reshape,dot,array,matrixmultiply

class GridPoint:
    def __init__(self,x,y,z,w=1.0,**opts):
        self.do_grad_dens = opts.get('do_grad_dens',False)
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)
        self._w = float(w)
        self.xyz = array((self._x,self._y,self._z),'d')
        self._r = sqrt(self._x*self._x+self._y*self._y+self._z*self._z)
        self._gamma = None
        self._dens = 0
        self._dens0 = None
        self._grad = None
        self.bfs = []
        return

    def xyzw(self): return (self._x,self._y,self._z,self._w)
    def dens(self): return self._dens
    def get_gamma(self): return self._gamma
    def grad(self): return self._grad
    def weight(self): return self._w
    def nbf(self): return len(self.bfs)
    def r(self): return self._r

    def setweight(self,w): self._w = w
    def zeroweight(self): self._w = 0

    def set_bf_amps(self,bfs):
        x,y,z,w = self.xyzw()
        nbf = len(bfs)
        self.bfs = zeros(nbf,'d')
        for i in xrange(nbf):
            self.bfs[i] = bfs[i].amp(x,y,z)
        # This *if* statement is potentially slow. If it becomes
        #  a factor, pull it up to AtomicGrids and have two
        #  explicit cases here.
        if self.do_grad_dens:
            self.bfgrads = zeros((nbf,3),'d')
            for i in xrange(nbf):
                self.bfgrads[i,:] = bfs[i].grad(x,y,z)
        return

    def floor_density(self,tol=1e-9):
        if self._dens < tol:
            self._dens = 0
        return

    def setdens(self,D):
        self._dens = 2*dot(self.bfs,matrixmultiply(D,self.bfs))
        # This *if* statement is potentially slow. If it becomes
        #  a factor, pull it up to AtomicGrids and have two
        #  explicit cases here.
        if self.do_grad_dens:
            self._grad = 2*dot(self.bfs.T,matrixmultiply(D,self.bfgrads)) +\
                         2*dot(self.bfgrads.T,matrixmultiply(D,self.bfs))
            self._gamma = dot(self._grad,self._grad)
        return

    def scale_density(self,factor):
        self._dens *= factor
        return

    def scale_weights(self,factor):
        self._w *= factor
        return


