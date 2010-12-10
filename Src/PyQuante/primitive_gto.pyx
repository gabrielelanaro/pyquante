#		 This source code is part of the
#		 PyQuante quantum chemistry suite
#
# Written by Gabriele Lanaro, 2009-2010
# Copyright (c) 2009-2010, Gabriele Lanaro
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

cimport primitive_gto
from stdlib cimport *

cdef class PrimitiveGTO:
    def __cinit__(self, *a, **kw):
        self.this = primitive_gto.primitive_gto_new()

    def __init__(self,double alpha, origin = (0.0,0.0,0.0), powers = (0,0,0), coef=1):
        cdef double x,y,z
        cdef int l,m,n
        x,y,z = origin
        l,m,n = powers
        primitive_gto.primitive_gto_init(self.this,alpha,x,y,z,l,m,n,coef)

    def amp(self, double x, double y, double z):
        return primitive_gto_amp(self.this, x,y,z)
    
    def __dealloc__(self):
        primitive_gto_free(self.this)
        
    property norm:
        def __get__(self):
            return self.this.norm
    property alpha:
        def __get__(self):
            return self.this.alpha
