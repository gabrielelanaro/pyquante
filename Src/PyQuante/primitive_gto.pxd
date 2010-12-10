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

cdef extern from "primitive-gto.h":
    ctypedef struct cPrimitiveGTO "PrimitiveGTO":
        double alpha
        double x0,y0,z0
        int l,m,n
        double norm
        double coef
    
    cPrimitiveGTO *primitive_gto_new()
    void primitive_gto_init( cPrimitiveGTO *gto, double alpha,
                             double x0, double y0, double z0,
                             int l, int m, int n, double coef)
    void primitive_gto_free(cPrimitiveGTO *pgto)
    void primitive_gto_normalize(cPrimitiveGTO *pgto)
    double primitive_gto_overlap(cPrimitiveGTO *pgto1, cPrimitiveGTO *pgto2)
    double primitive_gto_amp(cPrimitiveGTO *pgto, double x, double y, double z)

cdef class PrimitiveGTO:
    cdef cPrimitiveGTO *this
