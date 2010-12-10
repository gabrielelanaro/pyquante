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

from contracted_gto cimport ContractedGTO
cimport shell

sym2int = { "S" : 0,
            "P" : 1,
            "D" : 2,
            "F" : 3,
            "G" : 4,
            "H" : 5,
            "I" : 6}

cdef class Shell:
    def __cinit__(self, *a,**kw):
        self.this = shell_new()

    def __init__(self, ang_mom):
        shell_init(self.this, sym2int[ang_mom])
        self._cgtolist = [] # To hold references

    def append(self, ContractedGTO function, index):
        shell_append(self.this, function.this, index)
        self._cgtolist.append(function) # To hold references

    def __dealloc__(self):
        shell_free(self.this)

    property ang_mom:
        def __get__(self):
            return self.this.ang_mom

    property functions:
        def __get__(self):
            return

    property nfuncs:
        def __get__(self):
            return self.this.nfuncs

    property basis_index:
        def __get__(self):
            return [self.this.basis_index[i] for i in range(self.nfuncs)]
