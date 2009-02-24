"""\

NAME
      leapfrog.py

SYNOPSIS
      "leap frog" numerical integration for performing dynamics
      
DESCRIPTION
      This module provides a "leap frog" numerical integrator for doing 
      molecular dynamics.  This integration technique is described in the
      Feynman lectures, Vol 1, Section 9-6.  
	  
AUTHOR
      Hatem H. Helal, hhh23@cam.ac.uk

REPORT BUGS
      Report bugs to hhh23@cam.ac.uk

COPYRIGHT

"""

from NumWrap import array
"""
a note on units
assume forces are stored as H/Bohr, mass is AMU, and our time step is given in ps
to get quantities like acceleration in Bohr/ps^2 we need to convert our force to
units of AMU*Bohr/ps^2 using:
1 H    = 4.35974417(75)X10^-18 J
1 amu  = 1.660538782(83)X10^-27 Kg
1 Bohr = 5.291772108(18)X10^-11 m
we find 

1 H/Bohr = 937583.984 amu*Bohr/ps^2 
"""
unitfactor = 937583.984
def leapfrog(mol,t,dt):
    ke = 0.0 
    if t==0.: #initialize velocity at half step
        for atom in mol.atoms:
            m   = atom.mass()
            x   = atom.pos()          #x(0)
            v0  = atom.velocity() #v at t=0
            a0  = atom.force()/m*unitfactor
            
            #velocity at t=dt/2 is set as v(dt/2) = v(0) + (dt/2)*a(0)
            vnew = v0 + a0*dt*0.5
            atom.set_velocity(vnew)
            
            xnew = x+dt*vnew #now calculate x(dt)
            atom.update_coords(xnew)
            
            vavg = 0.5*(vnew+v0)
            ke+= 0.5*m*sum(vavg**2)/unitfactor
    else:
        for atom in mol.atoms:
            #leap frog integration is summarized as:
        	# x(t + dt)   = x(t) + dt*v(t+dt/2)
            # v(t + dt/2) = v(t - dt/2) + dt*a(t)
            # 	a(t)  = f(t)/m
            #so we use the velocities at the 1/2 steps to get positions at full steps
            m = atom.mass()
            x = atom.pos()          #x(t)
            v = atom.velocity()     #v(t-dt/2)
            a = atom.force()/m*unitfactor      #a(t)
            
            vnew = v+dt*a    #calculate v(t+dt/2)
            xnew = x+dt*vnew #now calculate x(t+dt)
            atom.set_velocity(vnew) 
            atom.update_coords(xnew)
            
            vavg = 0.5*(vnew+v)
            ke+= 0.5*m*sum(vavg**2)/unitfactor

    return ke
    

