"""\
 Atom.py: Simple class for atoms.

 This program is part of the PyQuante quantum chemistry program suite

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# My goal is to keep everything generic,i.e. not specific to
#  MINDO, HF, DFT, etc.

from NumWrap import array
from PyQuante.cints import dist2,dist
from Element import mass,symbol
from Constants import bohr2ang

# Careful about units! I'm not doing anything about them here;
#  whatever you store you get back.

class Atom:
    def __init__(self,atno,x,y,z,atid=0,fx=0.0,fy=0.0,fz=0.0,vx=0.0,vy=0.0,vz=0.0):
        self.atno = atno
        self.r = array([x,y,z],'d')
        #added by Hatem Helal hhh23@cam.ac.uk
        #atom id defaults to zero so as not to break preexisting code...
        self.atid = atid
        self.f = array([fx,fy,fz],'d')
        self.vel = array([vx,vy,vz],'d')
        return

    def __repr__(self): return "Atom ID: %d Atomic Num: %2d (%6.3f,%6.3f,%6.3f)" % \
        (self.atid,self.atno,self.r[0],self.r[1],self.r[2])
    def __getitem__(self, i):
        return self.r[i]
    def mass(self): return mass[self.atno]
    def pos(self): return tuple(self.r)
    def symbol(self): return symbol[self.atno]
    
    def force(self): return self.f
    def velocity(self): return self.vel
    
    # Could also do the following dists with numpy:
    def dist2(self,atom): return dist2(self.pos(),atom.pos())
    def dist(self,atom): return dist(self.pos(),atom.pos())
    def atuple(self): return (self.atno,self.r)
    def translate(self,pos): self.r += pos

    # The next two I've written as functions since if I ever handle
    #  pseudopotentials I'll need to do something clever, and this
    #  gives me a degree of indirection that will allow me to do this
    def get_nel(self): return self.atno
    def get_nuke_chg(self): return self.atno

    # This is set by the MINDO initialize routine. Will raise
    #  an error otherwise
    def get_nel_mindo(self): return self.Z

    def update_coords(self,xyz): self.r = array(xyz)
    def update_from_atuple(self,(atno,xyz)): self.update_coords(xyz)
    
    def set_force(self,fxfyfz): self.f = array(fxfyfz)
    def set_velocity(self,vxvyvz): self.vel = array(vxvyvz)

    def urotate(self,U):
        "Rotate molecule by the unitary matrix U"
        from PyQuante.NumWrap import matrixmultiply
        #self.r = matrixmultiply(self.r,U)
        self.r = matrixmultiply(U,self.r)
        return
    

def test():
    at1 = Atom(1,0,0,0)
    at2 = Atom(1,1,0,0)
    print at1.dist(at2)

if __name__ == '__main__': test()
