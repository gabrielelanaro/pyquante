#!/usr/bin/env python
"""\
  basis_converter.py Convert a Jaguar basis record to a native python format.

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

import sys,string

sym2no = {
    'X' : 0, 'H'  : 1, 'He' : 2,

    'Li' : 3, 'Be' : 4, 'B'  : 5, 'C'  : 6, 'N'  : 7,
    'O'  : 8, 'F'  : 9, 'Ne' : 10,
    
    'Na' : 11, 'Mg' : 12, 'Al' : 13, 'Si' : 14,
    'P'  : 15, 'S'  : 16, 'Cl' : 17, 'Ar' : 18,

    'K' : 19, 'Ca' : 20, 'Sc':21, 'Ti':22, 'V':23,'Cr':24,'Mn':25,
    'Fe' : 26, 'Co':27, 'Ni':28, 'Cu':29,'Zn':30,
    'Ga' : 31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,

    'Rb':37, 'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,
    'Ru' : 44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,
    'Sn':50,'Sb':51,'Te':52,'I':53,'Xe':54,

    'Cs':55,'Ba':56,'La':57,'Ce':58,'Pr':59,'Nd':60,'Pm':61,'Sm':62,
    'Eu':63,'Gd':64,'Tb':65,'Dy':66,'Ho':67,'Er':68,'Tm':69,'Yb':70,
    'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,
    'Au':79,'Hg':80,'Tl':81,'Pb':82,'Bi':83,'At':85,'Rn':86,
    'U' : 92,

    'x' : 0, 'h'  : 1, 'he' : 2,

    'li' : 3, 'be' : 4, 'b'  : 5, 'c'  : 6, 'n'  : 7,
    'o'  : 8, 'f'  : 9, 'ne' : 10,
    
    'na' : 11, 'mg' : 12, 'al' : 13, 'si' : 14,
    'p'  : 15, 's'  : 16, 'cl' : 17, 'ar' : 18,

    'k' : 19, 'ca' : 20, 'sc':21, 'ti':22, 'v':23,'cr':24,'mn':25,
    'fe' : 26, 'co':27, 'ni':28, 'cu':29,'zn':30,
    'ga' : 31,'ge':32,'as':33,'se':34,'br':35,'kr':36,

    'rb':37, 'sr':38,'y':39,'zr':40,'nb':41,'mo':42,'tc':43,
    'ru' : 44,'rh':45,'pd':46,'ag':47,'cd':48,'in':49,
    'sn':50,'sb':51,'te':52,'i':53,'xe':54,

    'cs':55,'ba':56,'la':57,'ce':58,'pr':59,'nd':60,'pm':61,'sm':62,
    'eu':63,'gd':64,'tb':65,'dy':66,'ho':67,'er':68,'tm':69,'yb':70,
    'lu':71,'hf':72,'ta':73,'w':74,'re':75,'os':76,'ir':77,'pt':78,
    'au':79,'hg':80,'tl':81,'pb':82,'bi':83,'at':85,'rn':86,
    'u' : 92,
} 

def main(filename="basis_631ss.dat"):
    outfilename = string.replace(filename,'.dat','.py')
    file = open(filename)
    bfs = []
    while 1:
        line = file.readline()
        if not line: break
        words = string.split(line)
        sym = words[0]
        atno = sym2no[sym]
        nat = len(bfs)
        if len(bfs) < atno+1:
            for i in range(atno+1-len(bfs)): bfs.append([])
        while 1:
            line = file.readline()
            if not line: break
            words = string.split(line)
            if len(words) < 1: break
            if words[0] == '****': break
            type,nprim = words[0],int(words[2])
            try:
                nprim2 = int(words[3])
                nprim = nprim + nprim2
            except:
                pass
            prims = []
            pprims = []
            for i in range(nprim):
                line = file.readline()
                words = string.split(line)
                expnt = float(words[0])
                coef = float(words[1])
                prims.append((expnt,coef))
                if type == 'SP':
                    coef2 = float(words[2])
                    pprims.append((expnt,coef2))
            if type == 'SP':
                bfs[atno].append(('S',prims))
                bfs[atno].append(('P',pprims))
            else:
                bfs[atno].append((type,prims))
    file.close()
    file = open(outfilename,'w')
    file.write('basis = [\n')
    for bf in bfs:
        if bf:
            file.write('    [\n')
            for type,prims in bf:
                file.write('    (\'%s\',[\n' % type)
                for expnt,coef in prims:
                    file.write('        (%f, %f),\n' % (expnt,coef))
                file.write('        ]),\n')
            file.write('    ],\n')
        else:
            file.write('    None,\n')
    file.write('    ]\n') 
                
if __name__ == '__main__':
    if len(sys.argv) < 2:
        main()
    else:
        main(sys.argv[1])

