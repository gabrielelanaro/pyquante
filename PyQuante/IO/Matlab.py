"""\
 Code for reading/writing Matlab file formats
 
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from NumWrap import zeros

def mtx2file(a,filename='bs.dat'):
    "General format for matrices: one space-delimited row per line"
    n,m = a.shape
    file = open(filename,'w')
    file.write("%d %d\n" % (n,m))
    for i in xrange(n):
        for j in xrange(m):
            file.write("%f " % a[i,j])
        file.write("\n")
    file.close()
    return

def rdmtx(filename):            
    file = open(filename)
    line = file.readline()
    n,m = map(int,line.split())
    A = zeros((n,m),'d')
    for i in xrange(n):
        line = file.readline()
        vals = map(float,line.split())
        for j in xrange(m):
            A[i,j] = vals[j]
    file.close()
    return A

def print_halfmat(A,name=None):
    "Print the lower half of a square matrix"
    if name: print "%s Matrix" % name
    for i in xrange(A.shape[0]):
        for j in xrange (i+1):
            print '%10.4f ' % A[i,j],
        print ''
    return


    
