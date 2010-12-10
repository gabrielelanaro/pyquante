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

cimport contracted_gto
import contracted_gto
from primitive_gto cimport PrimitiveGTO, cPrimitiveGTO
#from primitive_gto import PrimitiveGTO

from stdlib cimport *

cdef class ContractedGTO:
    '''
    Class that represents a linear combination of PrimitiveGTOs
    '''
    def __cinit__(self, origin, powers, atid=0):
        self.this = contracted_gto.contracted_gto_new()

    def __init__(self,origin, powers, atid=0):
        x,y,z = origin
        l,m,n = powers        
        contracted_gto.contracted_gto_init(self.this, x, y, z, l, m, n, atid)

    @classmethod
    def from_primitives(cls, pgtos, coefs, atid=0):
        raise NotImplementedError()
        cdef cPrimitiveGTO **pgto_array
        cdef int n,i
        cdef PrimitiveGTO cy_pgto
        
        n = <int>len(pgtos)
        # this array will be free in the cgto destructor 
        pgto_array = <cPrimitiveGTO **>malloc(n * sizeof(cPrimitiveGTO*))
        
        # Constructing the PrimitiveGTO** array...
        for i in range(n):
            cy_pgto = <PrimitiveGTO> pgtos[i]
            cy_pgto.this.coef  = <double>coefs[i] # set the coef
            
            pgto_array[i] = <cPrimitiveGTO *> cy_pgto.this
        
        #obj = ContractedGTO(x,y,z,l,m,n,atid)
        #contracted_gto.contracted_gto_from_primitives(obj.this , pgto_array, n)
        #return obj
        
    
    def __dealloc__(self):
        contracted_gto.contracted_gto_free(self.this)
        
    def add_primitive(self, primitive_gto.PrimitiveGTO pgto, double coeff):
        contracted_gto.contracted_gto_add_primitive(self.this, pgto.this, coeff)
        
    def amp(self,double x,double y,double z):
        '''
        given a point in the space it returns the amplitude of the
        GTO
        '''
        return contracted_gto.contracted_gto_amp(self.this,x,y,z)
    
    def normalize(self):
        contracted_gto.contracted_gto_normalize(self.this)
        
    property norm:
        
        def __get__(self):
            return self.this.norm
        def __set__(self, double value):
            self.this.norm = value
        
    property primitives:
        
        def __get__(self):
            cdef cPrimitiveGTO* prim
            ret = []
            for i in range(self.this.nprims):
                prim = <cPrimitiveGTO *>self.this.primitives[i]
                ret.append(PrimitiveGTO(prim.alpha,
                                        (prim.x0,prim.y0,prim.z0),
                                        (prim.l,prim.m,prim.n),
                                        prim.coef))
            return ret
                            
'''            
cpdef double coulomb( ContractedGTO a,
                      ContractedGTO b,
                      ContractedGTO c,
                      ContractedGTO d ):
    return contracted_gto.contracted_gto_coulomb(a.this,
                                                 b.this,
                                                 c.this,
                                                 d.this)
'''
