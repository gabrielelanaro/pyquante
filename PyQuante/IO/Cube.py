"""\
 Code for reading/writing Gaussian Cube files
 
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# This code doesn't work!!! By which I mean that it works even
# less than the normal code in PyQuante ;-).
#
# Still in the process of being debugged.

def get_bbox(atoms,**opts):
    dbuff = opts.get('dbuff',5)
    big = opts.get('big',10000)
    xmin = ymin = zmin = big
    xmax = ymax = zmax = -big
    for atom in atoms:
        x,y,z = atom.pos()
        xmin = min(xmin,x)
        ymin = min(ymin,y)
        zmin = min(zmin,z)
        xmax = max(xmax,x)
        ymax = max(ymax,y)
        zmax = max(zmax,z)
    xmin -= dbuff
    ymin -= dbuff
    zmin -= dbuff
    xmax -= dbuff
    ymax -= dbuff
    zmax -= dbuff
    return (xmin,xmax),(ymin,ymax),(zmin,zmax)

def mesh_orb(atoms,bfs,orbs,index):
    (xmin,xmax),(ymin,ymax),(zmin,zmax) = get_bbox(atoms)
    dx,dy,dz = xmax-xmin,ymax-ymin,zmax-zmin
    ppb = 2.0 # Points per bohr
    spacing = 1.0/ppb
    nx,ny,nz = int(dx*ppb),int(dy*ppb),int(dz*ppb)
    print "CUBE FILE"
    print "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z"
    print "%5i %11.6f %11.6f %11.6f" %  (len(atoms),xmin,ymin,zmin)
    print "%5i %11.6f %11.6f %11.6f" %  (nx,spacing,0,0)
    print "%5i %11.6f %11.6f %11.6f" %  (ny,0,spacing,0)
    print "%5i %11.6f %11.6f %11.6f" %  (nz,0,0,spacing)

    # The second record here is the nuclear charge, which differs from the
    #  atomic number when a ppot is used. Since I don't have that info, I'll
    #  just echo the atno
    for atom in atoms:
        atno = atom.atno
        x,y,z = atom.pos()
        print "%5i %11.6f %11.6f %11.6f %11.6f" %  (atno,atno,x,y,z)
    nbf = len(bfs)
    print " ",
    for i in xrange(nx):
        xg = xmin + i*spacing
        for j in xrange(ny):
            yg = ymin + j*spacing
            for k in xrange(nz):
                zg = zmin + k*spacing
                amp = 0
                for ibf in xrange(nbf):
                    amp += bfs[ibf].amp(xg,yg,zg)*orbs[ibf,index]
                if abs(amp) < 1e-12: amp = 0
                print " %11.5e" % amp,
                if k % 6 == 5: print "\n ",
            print "\n ",
