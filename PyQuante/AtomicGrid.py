"""\
 AtomicGrid.py Simple atomic grids, based on:
  A.D. Becke, 'A multicenter numerical integration scheme for
   polyatomic molecules.' J. Chem. Phys 88(4) 1988.

 The atomic grids are constructed from atomic grids that use
 Lebedev grids for the angular part, and Legendre grids for
 the radial parts.

 This program is part of the PyQuante quantum chemistry program suite

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from GridPoint import GridPoint
from Lebedev import Lebedev
from Legendre import Legendre
from math import sin,cos,pi
from NumWrap import zeros,array
from Constants import ang2bohr

# Where do these values come from? What are the units?
Bragg = [None,
     0.50, 2.00, 1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50,
     2.25, 1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 2.50,
     2.20, 1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35,
     1.35, 1.35, 1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 2.75,
     2.35, 2.00, 1.80, 1.55, 1.45, 1.45, 1.35, 1.30, 1.35,
     1.40, 1.60, 1.55, 1.55, 1.45, 1.45, 1.40, 1.40, 3.00,
     2.60, 2.15, 1.95, 1.85, 1.85, 1.85, 1.85, 1.85, 1.85,
     1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.55,
     1.45, 1.35, 1.35, 1.30, 1.35, 1.35, 1.35, 1.50, 1.90,
     1.80, 1.60, 1.90, 1.65, 3.25, 2.80, 2.15, 1.95, 1.80,
     1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75,
     1.75, 1.75, 1.75, 1.75, 1.55, 1.55
     ] 

# Pople atomic radii (bohr) for the SG-1 grid
PopleRadii = [None,
              1.0000, 0.5882, 3.0769, 2.0513, 1.5385, 1.2308,
              1.0256, 0.8791, 0.7692, 0.6839, 4.0909, 3.1579,
              2.5714, 2.1687, 1.8750, 1.6514, 1.4754, 1.3333]

class AtomicGrid:
    """\
    Put together an atomic grid consisting of a Legendre radial grid
    plus a Lebedev angular grid.

    The grid can be specified using a number of radial shells in
    [20,24,28,32,36] plus a fineness in range(4), or by using a list
    of angular momenta, where the number of radial shells is
    determined by the length of the list.
    """
    def __init__(self, atom, **opts):
        # Currently I'm keeping the grid in an array of points,
        # and in an array of shells of points. The goal is to move
        # to only having the array of shells 
        self.do_grad_dens = opts.get('do_grad_dens',False)
        radial = opts.get('radial','EulerMaclaurin')
        nrad = opts.get('nrad',32)
        fineness = opts.get('fineness',1)
        self.points = []
        self.shells = []
        self.Z = atom.atno
        self.Rmax = 0.5*Bragg[self.Z]*ang2bohr

        x,y,z = atom.pos()

        if radial == 'Legendre':
            self.grid = LegendreGrid(nrad,self.Rmax,fineness)
        else:
            self.grid = EulerMaclaurinGrid(nrad,self.Z,do_sg1=False)
            #self.grid = EulerMaclaurinGrid(50,self.Z,nang=194,do_sg1=False)
        # Could also call for EML or SG1 grids:

        for rrad,wrad,nangpts in self.grid:
            shell = []
            self.shells.append(shell)
            for (xang,yang,zang,wang) in Lebedev[nangpts]:
                weight = wrad*wang
                point = GridPoint(rrad*xang+x,rrad*yang+y,
                                  rrad*zang+z,weight,**opts)
                self.points.append(point)
                shell.append(point)
        return

    def __len__(self): return len(self.points)
    def __getitem__(self,index): return self.points[index]

    def floor_density(self,tol=1e-9):
        for point in self.points:
            point.floor_density(tol)
        return

    def set_bf_amps(self,bfs,**opts):
        "Set the basis function amplitude at each grid point"
        for point in self.points: point.set_bf_amps(bfs,**opts)
        return

    def setdens(self,D,**opts):
        "Set the density at each grid point"
        for point in self.points: point.setdens(D)
        return

    def weights(self):
        "Return a vector of weights of each point in the grid"
        weights = zeros(len(self.points),'d')
        for i in xrange(len(self.points)):
            weights[i] = self.points[i].weight()
        return weights
    
    def dens(self):
        "Return the density for each point in the grid"
        dens = zeros(len(self.points),'d')
        for i in xrange(len(self.points)):
            dens[i] = self.points[i].dens()
        return dens

    def get_gamma(self):
        "Return the density gradient gamma for each point in the grid"
        if not self.do_grad_dens: return None
        gamma = zeros(len(self.points),'d')
        for i in xrange(len(self.points)):
            gamma[i] = self.points[i].get_gamma()
        return gamma

    def bfs(self,i):
        "Return a vector of the product of two basis functions "
        " over the entire grid"
        bfs = zeros(len(self.points),'d')
        for j in xrange(len(self.points)):
            bfs[j] = self.points[j].bfs[i]
        return bfs

    def make_bfgrid(self):
        "Construct a matrix with bfs in columns over the entire grid, "
        " so that R[0] is the first basis function, R[1] is the second..."
        bfs = []
        for point in self.points:
            bfs.extend(point.bfs)
        return array(bfs)

    def nbf(self):
        return self.points[0].nbf()

    def npts(self): return len(self.points)

    def scale_weights(self,factor):
        for point in self.points:
            point.scale_weights(factor)
        return

    def scale_density(self,factor):
        for point in self.points:
            point.scale_density(factor)
        return

# The following two routines return [(ri,wi,nangi)] for nrad shells.
# The ri's are properly adjusted to go to the proper distances.
# The wi's are adjusted to only have to be multiplied by wrad from
# the lebedev shell
def EulerMaclaurinGrid(nrad,Z,**opts):
    do_sg1 = opts.get('do_sg1',True)
    nang = opts.get('nang',194)
    radial = EulerMaclaurinRadialGrid(nrad,Z)
    if do_sg1:
        grid = [(r,w,SG1Angs(r,Z)) for r,w in radial]
    else:
        grid = [(r,w,nang) for r,w in radial]
    return grid

def LegendreGrid(nrad,Rmax,fineness):
    #Rmax = 0.5*Bragg[Z]*ang2bohr

    radial = Legendre[nrad]
    grid = []
    for i in xrange(nrad):
        xrad,wrad = radial[i]
        rrad = BeckeRadMap(xrad,Rmax)
        dr = 2*Rmax/pow(1-xrad,2)
        vol = 4*pi*rrad*rrad*dr
        nangpts = ang_mesh(float(i+1)/nrad,fineness)
        grid.append((rrad,wrad*vol,nangpts))
    return grid
    
def BeckeRadMap(x,Rmax):
    return Rmax*(1.0+x)/(1.0-x)

def ang_mesh(frac,fineness,alevs = None):
    """\
    Determine the number of points in the angular mesh based on
    the fraction of the total radial grid index frac c (0,1).

    You can optionally pass in the number of points for
    the 5 different regions
    """
    if not alevs:
        ang_levels = [
            [ 6, 14, 26, 26, 14], # Coarse
            [ 50, 50,110, 50, 26], # Medium
            [ 50,110,194,110, 50], # Fine
            [194,194,194,194,194]  # ultrafine
            ]
        alevs = ang_levels[fineness]
    nang = alevs[0]
    if frac > 0.4: nang = alevs[1]
    if frac > 0.5: nang = alevs[2]
    if frac > 0.7: nang = alevs[3]
    if frac > 0.8: nang = alevs[4]
    return nang

def random_spin(angshell):
    "Spin an angular shell by a random amount"
    import random
    newshell = []

    x=2*pi*random.random()
    y=2*pi*random.random()
    z=2*pi*random.random()

    ang1=cos(z)
    ang2=sin(z)
    ang3=cos(y)
    ang4=sin(y)
    ang5=cos(x)
    ang6=sin(x)

    npts = len(angshell)
    # Rotations around the z axis affect x and y
    for i in xrange(npts):
        x,y,z,w = angshell[i]
        xtmp = ang1*x+ang2*y
        ytmp = ang1*y-ang2*x
        x,y = xtmp,ytmp

        xtmp = ang3*x+ang4*z
        ztmp = ang3*z-ang4*x
        x,z = xtmp,ztmp

        ytmp = ang5*y+ang6*z
        ztmp = ang5*z-ang6*y
        y,z = ytmp,ztmp
        newshell.append((x,y,z,w))

    return newshell

def regular_spin(angshell,i):
    "Spin an angular shell by a regular amount per shell"

    npts = len(angshell)
    deg_per_shell=5
    newshell = []
    
    imod = i % 4
    if imod == 0:
        return angshell
    elif imod == 1:
        x=2*pi*(deg_per_shell*i)/360.
        ang5=cos(x)
        ang6=sin(x)
        for i in xrange(npts):
            x,y,z,w = angshell[i]
            ytmp = ang5*y+ang6*z
            ztmp = ang5*z-ang6*y
            y,z = ytmp,ztmp
            newshell.append((x,y,z,w))
        return newshell
    elif imod == 2:
        y=2*pi*(deg_per_shell*i)/360.
        ang3=cos(y)
        ang4=sin(y)
        for i in xrange(npts):
            x,y,z,w = angshell[i]
            xtmp = ang3*x+ang4*z
            ztmp = ang3*z-ang4*x
            x,z = xtmp,ztmp
            newshell.append((x,y,z,w))
        return newshell
    elif imod == 3:
        z=2*pi*(deg_per_shell*i)/360.
        ang1=cos(z)
        ang2=sin(z)
        for i in xrange(npts):
            x,y,z,w = angshell[i]
            xtmp = ang1*x+ang2*y
            ytmp = ang1*y-ang2*x
            x,y = xtmp,ytmp
            newshell.append((x,y,z,w))
        return newshell
    return newshell

def EulerMaclaurinRadialGrid(nrad,Z):
    # Radial part of the Gill, Johnson, Pople SG-1 grid
    R = PopleRadii[Z]
    grid = []
    for i in xrange(1,nrad+1):
        # Changed to include a factor of 4pi
        #w = 2.*pow(R,3)*(nrad+1.)*pow(i,5)*pow(nrad+1-i,-7)
        w = 8.*pi*pow(R,3)*(nrad+1.)*pow(i,5)*pow(nrad+1-i,-7)
        r = R*i*i*pow(nrad+1-i,-2)
        grid.append((r,w))
    return grid

def SG1Angs(r,Z):
    # Gill, Johnson, Pople rules for SG-1 angular densities
    R = PopleRadii[Z]
    if Z in xrange(1,3): # H-He
        alphas = [0.25,0.5,1.0,4.5]
    elif Z in xrange(3,11): # Li-Ne
        alphas = [0.1667, 0.500, 0.900, 3.5]
    else: # only fit for Na-Ar
        alphas = [0.1,0.4,0.8,2.5]

    if r < alphas[0]*R: return 6
    elif r < alphas[1]*R: return 38
    elif r < alphas[2]*R: return 86
    elif r < alphas[3]*R: return 194
    return 86
