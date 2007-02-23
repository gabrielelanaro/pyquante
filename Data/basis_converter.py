#!/usr/bin/env python
"""\
 basis_converter.py Convert a GAMESS-US basis set to PyQuante format

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
from pprint import pprint
from PyQuante.Basis.Tools import parse_gamess_basis

header = """\
 %s basis set for use with PyQuante

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

def main():
    import sys
    #fnames = ["basis_321.dat","basis_6311pp_2d_2p.dat",
    #          "basis_6311ss.dat","basis_ccpvdz.dat",
    #          "basis_ccpvtz.dat","basis_dzvp.dat",
    #          "basis_sto3g.dat","basis_sto6g.dat"]
    fnames = sys.argv[1:]
    do_write = True
    for fname in fnames:
        print "Processing ",fname
        oname = fname.replace(".dat",".py")
        basis = parse_gamess_basis(fname)
        if do_write:
            file = open(oname,"w")
            file.write('"""\\\n')
            file.write(header % fname)
            file.write('"""\n')
            file.write("basis_data = \\\n")
            pprint(basis,file)
            file.close()
        else:
            pprint(basis)
    return

if __name__ == "__main__": main()
    
