"""\
 Molf - Module containing functions for reading/writing Molden MOLF files

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 

"""

import re
from PyQuante.NumWrap import zeros
from PyQuante.Util import parseline

def section(fname):
    title_lines = []
    atom_lines = []
    basis_lines = []
    orb_lines = []
    misc_lines = []
    recipient = misc_lines
        
    for line in open(fname):
        if line.startswith("[Title"): recipient = title_lines
        elif line.startswith("[Atoms"): recipient = atom_lines
        elif line.startswith("[GTO"): recipient = basis_lines
        elif line.startswith("[MO"): recipient = orb_lines
        recipient.append(line)
    return title_lines,atom_lines,basis_lines,orb_lines

def parse_atoms(atom_lines,**opts):
    from PyQuante.Molecule import Molecule
    
    atomlist = []
    for line in atom_lines[1:]:
        atno,x,y,z = parseline(line,'xxdfff')
        atomlist.append((atno,(x,y,z)))

    molopts = {}
    if re.search("Angs",atom_lines[0]): molopts['units'] = 'angs'
    atoms = Molecule('Molf',atomlist,**molopts)

    return atoms

def parse_basis(lines):
    bfs = []
    while 1:
        try:
            line = lines.pop(0)
        except:
            break
        if not line: break
        words = line.split()
        if not words: break
        iat,ibs = map(int,words)
        atbasis = []
        while 1:
            line = lines.pop(0)
            words = line.split()
            if not words: break
            type = words[0]
            nprim = int(words[1])
            coef = float(words[2])
            prims = []
            for i in xrange(nprim):
                line = lines.pop(0)
                exp,coef = map(float,line.split())
                prims.append((exp,coef))
            atbasis.append((type,prims))
        bfs.append((iat,atbasis))
    return bfs

def parse_orbs(lines,nbf):
    sympat = re.compile('^Sym\s*=\s*(\S+)\s*$')
    enepat = re.compile('^Ene\s*=\s*(\S+)$')
    spinpat = re.compile('^Spin\s*=\s*(\S+)$')
    occpat = re.compile('^Occup\s*=\s*(\S+)$')
    coefpat = re.compile('\s*(\d+)\s*(\S+)$')

    orbs = []
    sym,ene,spin,occ,coefs = None,None,None,None,[]

    for line in lines:
        if sympat.search(line):
            if sym: orbs.append((sym,ene,spin,occ,coefs))
            sym,ene,spin,occ,coefs = None,None,None,None,[]
            sym = sympat.search(line).groups()[0]
        elif enepat.search(line):
            ene = float(enepat.search(line).groups()[0])
        elif spinpat.search(line):
            spin = spinpat.search(line).groups()[0]
        elif occpat.search(line):
            occ = float(occpat.search(line).groups()[0])
        elif coefpat.search(line):
            a,b = coefpat.search(line).groups()
            iorb,coef = int(a),float(b)
            coefs.append((iorb,coef))
    #print_orbs(orbs)
    norb = len(orbs)
    orbmtx = zeros((nbf,norb),'d')
    orbe = zeros(norb,'d')
    occs = zeros(norb,'d')
    for i in xrange(norb):
        sym,ene,spin,occ,coefs = orbs[i]
        occs[i] = occ
        orbe[i] = ene
        for ibf,coef in coefs:
            orbmtx[ibf-1,i] = coef
    return occs,orbe,orbmtx

def make_pyquante_basis(atoms,input_bfs):
    from PyQuante.Ints import sym2powerlist
    from PyQuante.CGBF import CGBF
    bfs = []
    for atom_index,info in input_bfs:
        x,y,z = atoms[atom_index-1].pos()
        for type,primlist in info:
            for power in sym2powerlist[type]:
                bf = CGBF((x,y,z),power)
                for expnt,coef in primlist:
                    bf.add_primitive(expnt,coef)
                bf.normalize()
                bfs.append(bf)
    return bfs
            
def read_molf(fname):
    title_lines,atom_lines,basis_lines,orb_lines = section(fname)
    # Atoms is the only parser that takes the first line (the line with the [] label)
    atoms = parse_atoms(atom_lines)
    bfs = parse_basis(basis_lines[1:])
    bfs = make_pyquante_basis(atoms,bfs)
    nbf = len(bfs)
    occs,orbe,orbs = parse_orbs(orb_lines[1:],nbf)
    return atoms,bfs,occs,orbe,orbs

