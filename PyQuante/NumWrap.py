#!/usr/bin/env python
"""
 NumWrap.py - Interface to Numeric and numpy

 An interface to the Numeric and numpy libraries that will (hopefully)
 make the transformation to numpy go as seamlessly as possible.

 Also interfaces the LinearAlgebra and numpy.linalg libraries, since
 these have to be done consistently with the Numeric/numpy choice
"""

use_numpy = True
import re
pat = re.compile('\D')

if use_numpy:
    from numpy import array,zeros,concatenate,dot,ravel,arange
    from numpy import arcsinh,diagonal,identity,choose,transpose
    from numpy import reshape,take
    from numpy import where
    from numpy import __version__ as version
    from numpy.oldnumeric import NewAxis

    words = map(int,pat.split(version))
    big_version = 100*words[0] + 10*words[1]
    if len(words) > 2: big_version += words[0]

    if big_version >= 100:
        from numpy.oldnumeric.linear_algebra import determinant
        from numpy.oldnumeric.linear_algebra import Heigenvectors
        from numpy.oldnumeric.linear_algebra import solve_linear_equations
        import numpy.oldnumeric.mlab as MLab
    elif big_version >= 98:
        from numpy.linalg.old import determinant
        from numpy.linalg.old import Heigenvectors
        from numpy.linalg.old import solve_linear_equations
        import numpy.lib.mlab as MLab
    else:
        from numpy.linalg import determinant
        from numpy.linalg import Heigenvectors
        from numpy.linalg import solve_linear_equations
        import numpy.lib.mlab as MLab
    matrixmultiply = dot
    import numpy as Numeric
    
else:
    from Numeric import array,zeros,concatenate,dot,ravel,matrixmultiply
    from Numeric import arange
    from Numeric import arcsinh,diagonal,identity,choose,transpose
    from Numeric import reshape,take
    from Numeric import where
    from Numeric import NewAxis
    from LinearAlgebra import solve_linear_equations,Heigenvectors
    from LinearAlgebra import determinant
    import Numeric
    import MLab

