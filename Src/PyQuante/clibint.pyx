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

from shell cimport Shell,cShell
cimport shell

from stdlib cimport *
from PyQuante.cints import ijkl2intindex as intindex

# My Libint library
cdef extern from "clibint.h":
    int cshell_compute_eri "shell_compute_eri" (cShell *shell1, cShell *shell2, cShell *shell3, cShell *shell4, double *output)

# Libint initialization, should be done once.
cdef extern from "libint.h":
   void init_libint_base()
init_libint_base()

def shell_compute_eri(Shell shell1,Shell shell2,Shell shell3,Shell shell4, Ints):
    '''
    Compute a Cartesian shell putting the result in the Ints array.
    The position of each ERI is given by the
    PyQuante.cints.ijkl2intindex function.

    The index of each basis function is the index that it has in the
    original basis set.
    '''
    cdef double *output
    cdef int n # Total size
    cdef int i,j,k,l # Four index
    
    n =  shell1.nfuncs*shell2.nfuncs*shell3.nfuncs*shell4.nfuncs
    
    output = <double *>malloc(n * sizeof(double))

    cshell_compute_eri(shell1.this,
                   shell2.this,
                   shell3.this,
                   shell4.this, output)
    
    basis_index1 = shell1.basis_index
    basis_index2 = shell2.basis_index
    basis_index3 = shell3.basis_index
    basis_index4 = shell4.basis_index
    
    for i in range(shell1.nfuncs):
        for j in range(shell2.nfuncs):
            for k in range(shell3.nfuncs):
                for l in range(shell4.nfuncs):                
                    basis_index = (basis_index1[i], basis_index2[j],
                                   basis_index3[k], basis_index4[l])
                    eri_index = ((i*shell2.nfuncs + j)*shell3.nfuncs+k)*shell4.nfuncs + l
                    Ints[intindex(*basis_index)] = output[eri_index]
                    
    free(output)
