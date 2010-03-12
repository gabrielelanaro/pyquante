"""\
PyQuante: Python Quantum Chemistry
 PyQuante is a toolkit for building quantum chemistry programs in the
 Python programming language.

 Modules:
 CGBF         Class/functions for manipulating contracted Gaussian bfns
 PGBF         Class/functions for manipulating primitive Gaussian bfns
 hartree_fock Functions for simple closed-shell Hartree-Fock calculations
 LA2          Add-ons to the LinearAlgebra package
 pyints       Basic one- and two-electron integrals written in python
 cints        C version of the functions in pyints
 rys          Two-electron integrals using Rys quadrature
 crys         C version of the functions in rys
 hgp          Two-electron ints using Head-Gordon/Pople recursion relations
 chgp         C version of the functions in hgp
 basis_*      Different basis sets
 dft          Module for closed-shell LDA density functional theory
 DFunctionals Different density functionals (right now only S & VWN)

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Starting to think about making modules available from main-level
# import. These are good candidates

from PyQuante.Molecule import Molecule
from PyQuante.PyQuante2 import SCF
from PyQuante.CGBF import CGBF
from PyQuante.logger import configure_output
