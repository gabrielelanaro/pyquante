#!/usr/bin/env python
"""\
 Molecule.py: Simple class for molecules.

 TODO: *Really* need to think of a more intelligent way of handling
 units!

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from PyQuante.Atom import Atom
from PyQuante.Element import sym2no
from PyQuante.Constants import bohr2ang,ang2bohr

allowed_units = ['bohr','angs']

class Molecule:
    """\
    Molecule(name,atomlist,**opts) - Object to hold molecular data in PyQuante

    name       The name of the molecule
    atomlist   A list of (atno,(x,y,z)) information for the molecule

    Options:      Value   Description
    --------      -----   -----------
    units         Bohr    Units for the coordinates of the molecule
    charge        0       The molecular charge
    multiplicity  1       The spin multiplicity of the molecule
    """
    def __init__(self,name='molecule',atomlist=[],**opts):
        self.name = name
        self.atoms = []
        self.basis = []
        self.grid = None
        units = opts.get('units','bohr')
        units = units.lower()[:4]
        assert units in allowed_units
        self.units = units
        if atomlist: self.add_atuples(atomlist)
        self.charge = opts.get('charge',0)
        self.multiplicity = opts.get('multiplicity',1)
        return

    def __repr__(self):
        outl = "\n%s: %s" % (self.__class__.__name__,self.name)
        for atom in self.atoms:
            outl += "\n\t%s" % str(atom)
        return outl

    def update_from_atuples(self,geo):
        for i in range(len(geo)):
            self.atoms[i].update_from_atuple(geo[i])
        return

    def translate(self,pos):
        for atom in self.atoms: atom.translate(pos)
        return

    def add_atom(self,atom): self.atoms.append(atom)
    
    def add_atuple(self,atno,xyz):
        if self.units != 'bohr': xyz = toBohr(xyz[0],xyz[1],xyz[2])
        if type(atno) == type(''): atno = sym2no[atno]
        self.atoms.append(Atom(atno,xyz[0],xyz[1],xyz[2]))

    def add_atuples(self,atoms):
        "Add a list of (atno,(x,y,z)) tuples to the atom list"
        from Atom import Atom
        for atno,xyz in atoms: self.add_atuple(atno,xyz)
        return

    def atuples(self):
        "Express molecule as a list of (atno,(x,y,z)) tuples"
        atoms = []
        for atom in self.atoms: atoms.append(atom.atuple())
        return atoms

    def atuples_angstrom(self):
        atoms = []
        for atom in self.atoms:
            atno,xyz = atom.atuple()
            atoms.append((atno,toAng(xyz[0],xyz[1],xyz[2])))
        return atoms

    # Not really used. Consider removing
    def add_xyz_file(self,filename,which_frame=-1):
        "Input atoms from xyz file. By default choose the last frame"
        from IO import read_xyz
        geos = read_xyz(filename)
        self.add_atuples(geos[which_frame])
        return

    def set_charge(self,charge): self.charge = int(charge)
    def get_charge(self): return self.charge

    def get_nel(self,charge=None):
        if charge:
            # Deprecation warning inserted 8/2005
            print "Warning: use of mol.get_nel(charge) has been deprecated"
            print "Please supply charge at construction of molecule or use"
            print "mol.set_charge(charge)"
            self.set_charge(charge)
        nel = -self.charge
        for atom in self.atoms: nel += atom.get_nel()
        return nel

    def get_enuke(self):
        enuke = 0.
        nat = len(self.atoms)
        for i in range(nat):
            ati = self.atoms[i]
            for j in range(i):
                atj = self.atoms[j]
                enuke += ati.get_nuke_chg()*atj.get_nuke_chg()/ati.dist(atj)
        return enuke

    def get_closedopen(self):
        multiplicity = self.multiplicity
        nel = self.get_nel()

        assert multiplicity > 0

        if (nel%2 == 0 and multiplicity%2 == 0) \
               or (nel%2 == 1 and multiplicity%2 == 1):
            print "Incompatible number of electrons and spin multiplicity"
            print "nel = ",nel
            print "multiplicity = ",multiplicity
            raise "Incompatible number of electrons and spin multiplicity"

        nopen = multiplicity-1
        nclosed,ierr = divmod(nel-nopen,2)
        if ierr:
            print "Error in Molecule.get_closedopen()"
            print 'nel = ',nel
            print 'multiplicity = ',multiplicity
            print 'nopen = ',nopen
            print 'nclosed = ',nclosed
            print 'ierr = ',ierr
            raise "Error in Molecule.get_closedopen()"
        return nclosed, nopen

    def get_alphabeta(self,**opts):
        nclosed,nopen = self.get_closedopen(**opts)
        return nclosed+nopen,nclosed

    # These two overloads let the molecule act as a list of atoms
    def __getitem__(self,i):return self.atoms[i]
    def __len__(self): return len(self.atoms)

    # These overloads let one create a subsystem from a list of atoms
    def subsystem(self,name,indices,**opts):
        submol = Molecule(name,None,**opts)
        for i in indices: submol.add_atom(self.atoms[i])
        return submol

def toBohr(*args):
    if len(args) == 1: return ang2bohr*args[0]
    return [ang2bohr*arg for arg in args]
def toAng(*args):
    if len(args) == 1: return bohr2ang*args[0]
    return [bohr2ang*arg for arg in args]

def cleansym(s):
    import re
    return re.split('[^a-zA-z]',s)[0]

def ParseXYZLines(name,xyz_lines,**opts):
    atoms = []
    for line in xyz_lines.splitlines():
        words = line.split()
        if not words: continue
        sym = cleansym(words[0])
        xyz = map(float,words[1:4])
        atoms.append((sym,xyz))
    return Molecule(name,atoms,**opts)

if __name__ == '__main__':
    h2o = Molecule('h2o',
                   [('O',(0.,0.,0.)),('H',(1.,0.,0.)),('H',(0.,1.,0.))],
                   units='Angstrom')
    print h2o
    print h2o.subsystem([0,1])
    
