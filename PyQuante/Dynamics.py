"""\
 Dynamics.py: Module for molecular dynamics

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
from NumWrap import array,zeros
from Constants import Rgas
from math import sqrt,pow
from IO import append_xyz

# Derivation of units for Force constant: (Thanks to Alejandro Strachan)
# We have accel in (kcal/mol)/(A*g/mol) and we want it in A/ps^2
# (kcal/mol)/(A*g/mol) = kcal/(A*g) = 1000 kcal/(A*kg)
#   = 1000 kcal/(A*kg) * 4.184 J/cal = 4184 kJ/(A*kg)
#   = 4.184e6 (kg m^2/s^2)/(A*kg) = 4.184e6 m^2/s^2/A
#   = 4.184e26 A/s^2 = 418.4 A/(ps^2)
fconst = 418.4 # convert (kcal/mol)/(A*g/mol) to A/ps^2
# The inverse of this quantity transforms back from amu*(A^2/ps^2) to kcal/mol

def Dynamics(atoms,EnergyForces,nsteps=1000,Ti=298,dt=1e-3):
    xyz = open('pyqmd.xyz','w')
    dat = open('pyqmd.dat','w')
    set_boltzmann_velocities(atoms,Ti)
    Etot = 0
    for step in range(nsteps):
        append_xyz(xyz,atoms.atuples(),"PQMD %4d  E = %10.4f" % (step,Etot))
        try:
            Ev,F = EnergyForces(atoms)
        except:
            print "Using averaging to try and converge"
            Ev,F = EnergyForces(atoms,0.5)
        set_forces(atoms,F)
        LeapFrogUpdate(atoms,dt)
        #for atom in atoms: flask.bounce(atom)
        Ek = get_kinetic(atoms)
        T = get_temperature(atoms)
        #rescale_velocities(atoms,Ti) # uncomment for iso-kinetics
        Etot = Ev+Ek
        print step*dt,Etot,Ev,Ek,T
        dat.write("%10.4f %10.4f %10.4f %10.4f %10.4f\n" %
                  (step*dt,Etot,Ev,Ek,T))
        dat.flush()
    return

def get_kinetic(atoms):
    sum_mv2 = 0
    for atom in atoms: sum_mv2 += atom.mass()*atom.v0.squared()
    return 0.5*sum_mv2/fconst

# There's a disconnect here, in that the kinetic energy is being
#  computed with v0 (v(t)) and the temperature is being computed
#  at v (v(t+dt/2))
def get_temperature(atoms):
    sum_mv2 = 0
    for atom in atoms: sum_mv2 += atom.mass()*atom.v.squared()
    return 1000*sum_mv2/((3*len(atoms)-6)*Rgas*fconst) 
    
def LeapFrogUpdate(atoms,dt):
    # Leap-frog Verlet dynamics is based on the equations
    #  v(t+dt/2) = v(t-dt/2)+dt*a(t)
    #  r(t+dt) = r(t) + dt*v(t+dt/2)
    # so that the positions, calculated at dt,2dt,3dt, etc.,
    # leap-frog over the velocities, calculated at dt/2,3dt/2,5dt/2...
    for atom in atoms:
        m = atom.mass()
        a = -atom.F*fconst/m        # a = F/m
        vnew = atom.v + dt*a        # v(t+dt/2) = v(t-dt/2) + dt*a
        # Save the current velocity for later calc of T,Ek
        atom.v0 = 0.5*(vnew+atom.v) # v(t) = 0.5*(v(t-dt/2)+v(t+dt/2)
        atom.r += dt*vnew           # r(t+dt) = r(t) + dt*v(t+dt/2)
        atom.v = vnew
    return

def set_forces(atoms,F):
    for i in range(len(atoms)):
        fx,fy,fz = F[i]
        atoms[i].F = array((fx,fy,fz))
    return

def set_boltzmann_velocities(atoms,T):
    from random import gauss,randint
    Eavg = Rgas*T/2000 # kT/2 per degree of freedom (kJ/mol)

    vels = []
    for atom in atoms:
        m = atom.mass()
        vavg = sqrt(2*Eavg*fconst/m)
        
        stdev = 0.01 #I'm setting the std dev wrong here
        atom.v = array((pow(-1,randint(0,1))*gauss(vavg,stdev),
                        pow(-1,randint(0,1))*gauss(vavg,stdev),
                        pow(-1,randint(0,1))*gauss(vavg,stdev)))

    subtract_com_velocity(atoms)
    rescale_velocities(atoms,T)
    return

def subtract_com_velocity(atoms):
    vcom = get_vcom(atoms)
    for atom in atoms: atom.v -= vcom
    return

def rescale_velocities(atoms,T):
    Tact = get_temperature(atoms)
    scalef = sqrt(T/Tact)
    for atom in atoms: atom.v *= scalef
    return

def get_vcom(atoms):
    "Compute the Center of Mass Velocity"
    vcom = zeros(3,'d')
    totm = 0
    for atom in atoms:
        m = atom.mass()
        vcom += m*atom.v
        totm += m
    return vcom/totm

if __name__ == '__main__':
    from MINDO3 import get_energy_forces
    from Molecule import Molecule
    rdx = Molecule('RDX',filename='/home/rmuller/gallery/rdx.xyz')
    Dynamics(rdx,get_energy_forces,nsteps=3,Ti=4000)
