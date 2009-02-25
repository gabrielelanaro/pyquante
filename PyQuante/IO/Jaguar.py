"""\
 Routines to read/write information from/to Jaguar files
 
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

def get_guess_lines(fname):
    lines = []
    started = False
    for line in open(fname):
        if started:
            if line.startswith('&'): break
            lines.append(line)
        if line.startswith('&guess'):
            started = True
    return lines

def split_guess_lines(guess_lines):
    from PyQuante.Util import parseline
    from PyQuante.NumWrap import transpose,array
    import re
    orbpat = re.compile('Orbital Energy')
    orb = []
    orbe = []
    occs = []
    orbs = []
    for line in guess_lines:
        if orbpat.search(line):
            orbei,occi = parseline(line,'xxxfxf')
            orbe.append(orbei)
            occs.append(occi)
            orb = []
            orbs.append(orb)
        else:
            orb.extend(map(float,line.split()))
    orbs = transpose(array(orbs))
    return occs,orbe,orbs

def orbs_from_restart(fname):
    guess_lines = get_guess_lines(fname)
    occs,orbe,orbs = split_guess_lines(guess_lines)
    return occs,orbe,orbs

def geo_from_output(fname):
    import re
    from PyQuante.Util import parseline,cleansym
    from PyQuante.Element import sym2no
    from PyQuante.Molecule import Molecule
    igeo = re.compile('Input geometry')
    sgeo = re.compile('Symmetrized geometry')
    # Double check the syntax of these last two
    ngeo = re.compile('new geometry')
    fgeo = re.compile('final geometry')
    geo = []
    file = open(fname)
    while 1:
        line = file.readline()
        if not line: break
        if igeo.search(line) or sgeo.search(line) or ngeo.search(line) \
           or fgeo.search(line):
            geo = []
            line = file.readline()
            line = file.readline()
            while 1:
                line = file.readline()
                if len(line.split()) < 4: break
                sym,x,y,z = parseline(line,'sfff')
                atno = sym2no[cleansym(sym)]
                geo.append((atno,(x,y,z)))
    return Molecule('jaguar molecule',geo)
