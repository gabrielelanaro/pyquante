#!/usr/bin/env python
"""\
 Generate the Legendre polynomial roots and weights using the SciPy library

 This program is part of the PyQuante quantum chemistry suite.
 PyQuante is copyright (c) 2002 Richard P. Muller. All Rights Reserved.
 You may contact the author at rpm@wag.caltech.edu.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307

"""
from scipy.special.orthogonal import p_roots

r= {}
w= {}
r[20],w[20] = p_roots(20)
r[24],w[24] = p_roots(24)
r[28],w[28] = p_roots(28)
r[32],w[32] = p_roots(32)
r[36],w[36] = p_roots(36)

print "Legendre = {"
for i in [20,24,28,32,36]:
    print "    %d: [" % i
    roots,weights = r[i],w[i]
    for ii in range(i):
        if roots[ii].imag:
            raise "Error, P_%d root %d complex"
        print "        (%20.15f,%20.15f)," % (roots[ii].real,weights[ii])
    print "    ],"
print "}"

