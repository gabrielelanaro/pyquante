#!/usr/bin/env python
"""\
 Code for miscellaneous input/output formats in PyQuante.
 
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from NumWrap import zeros

def mtx2file(a,filename='bs.dat'):
    "General format for matrices: one space-delimited row per line"
    n,m = a.shape
    file = open(filename,'w')
    file.write("%d %d\n" % (n,m))
    for i in range(n):
        for j in range(m):
            file.write("%f " % a[i,j])
        file.write("\n")
    file.close()
    return

def rdmtx(filename):            
    file = open(filename)
    line = file.readline()
    n,m = map(int,line.split())
    A = zeros((n,m),'d')
    for i in range(n):
        line = file.readline()
        vals = map(float,line.split())
        for j in range(m):
            A[i,j] = vals[j]
    file.close()
    return A

def print_halfmat(A,name=None):
    "Print the lower half of a square matrix"
    if name: print "%s Matrix" % name
    for i in range(A.shape[0]):
        for j in range (i+1):
            print '%10.4f ' % A[i,j],
        print ''
    return

def read_xyz(filename):
    from Element import sym2no
    file = open(filename)
    geometries = []
    while 1:
        line = file.readline()
        if not line: break
        nat = int(line.split()[0])
        title = file.readline()
        atoms = []
        for i in range(nat):
            line = file.readline()
            words = line.split()
            atno = sym2no[words[0]]
            x,y,z = map(float,words[1:])
            atoms.append((atno,(x,y,z)))
        geometries.append(atoms)
    return geometries

def write_xyz(filename,geometries,title="File written by XYZ.py"):
    file = open(filename,'w')
    for atoms in geometries: append(file,atoms,title)
    file.close()
    return

def append_xyz(file,atoms,title="File written by XYZ.py"):
    from Element import symbol
    file.write("%d\n%s\n" % (len(atoms),title))
    for atno,(x,y,z) in atoms:
        file.write("%4s %10.4f %10.4f %10.4f\n"
                   % (symbol[atno],x,y,z))
    file.flush()
    return

def read_jaguar_restart(filename,nbf):
    "Gets the orbitals from a Jaguar restart file"
    import re
    file = open(filename)
    guess_pat = re.compile('&guess')
    end_pat = re.compile('&')
    reading = 0
    for line in file.readlines():
        if reading:
            if end_pat.search(line):
                reading = 0
                break
            guess_section.append(line)
        if not reading:
            if guess_pat.search(line):
                reading = 1
                guess_section = []
    file.close()
    orbe = []
    orbs = zeros((nbf,nbf),'d')
    orbptr = -1
    orb_pat = re.compile('Orbital Energy')
    for line in guess_section:
        if orb_pat.search(line):
            basptr = 0
            orbptr += 1
        else:
            words = map(float,line.split())
            nwords = len(words)
            for i in range(nwords):
                orbs[orbptr,basptr+i] = words[i]
            basptr += nwords
    return orbs

