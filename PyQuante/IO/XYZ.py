"""\
 Code for reading/writing Xmol XYZ files.
 
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

def read_xyz(filename):
    from PyQuante.Element import sym2no
    from PyQuante.Molecule import Molecule
    file = open(filename)
    geometries = []
    igeo = 1
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
        atoms = Molecule("XYZ geometry #%d" % igeo,atoms)
        igeo += 1
        geometries.append(atoms)
    return geometries # Now a list of PyQuante molecules

def write_xyz(filename,geometries,title="File written by XYZ.py"):
    file = open(filename,'w')
    for atoms in geometries: append_xyz(file,atoms,title)
    file.close()
    return

def append_xyz(file,atoms,title="File written by XYZ.py"):
    from Element import symbol
    file.write("%d\n%s\n" % (len(atoms),title))
    for atom in atoms:
        atno,(x,y,z) = atom.atuple()
        file.write("%4s %10.4f %10.4f %10.4f\n"
                   % (symbol[atno],x,y,z))
    file.flush()
    return
