#!/usr/bin/env python
# Test the cint module

from PyQuante.PGBF import PGBF
from PyQuante.cints import coulomb_repulsion

p1 = PGBF(0.5,(0,0,0),(0,0,0))
p2 = PGBF(1.0,(0,0,1.),(1,0,0))
p3 = PGBF(0.1,(0,0,1.),(1,1,0))
p4 = PGBF(0.2,(0,1.,0),(0,0,1))
p5 = PGBF(0.5,(0,0,1.),(0,0,0))

print "Norms: ",p1.norm(),p2.norm(),p3.norm(),p4.norm()

print "Testing overlap",p1.overlap(p1),p1.overlap(p5)

print "Testing kinetic energy"
print p1.kinetic(p1),p2.kinetic(p2)

print "Testing cints coulomb repulsion"
print coulomb_repulsion(p1.origin(),p1.norm(),p1.powers(),p1.exp(),
                        p2.origin(),p2.norm(),p2.powers(),p2.exp(),
                        p3.origin(),p3.norm(),p3.powers(),p3.exp(),
                        p4.origin(),p4.norm(),p4.powers(),p4.exp())



