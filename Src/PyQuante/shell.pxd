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

cdef extern from "shell.h":
    ctypedef struct cShell "Shell":
        contracted_gto.cContractedGTO *functions
        int *basis_index
        int ang_mom
        double r[3]
        int nfuncs
    cShell * shell_new()
    void shell_init(cShell *shell, int ang_mom)
    void shell_append (cShell *shell, contracted_gto.cContractedGTO *cgto, int index)
    void shell_free(cShell *shell)

cdef class Shell:
    cdef cShell *this
    cdef object _cgtolist # Just to handle reference counting
