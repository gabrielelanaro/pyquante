#!/usr/bin/env python

from PyQuante.NumWrap import arange
from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

energy = -1.132710
name = 'H2/DFT'

def dissoc(**opts):
    rmin = opts.get('rmin',0.65)
    rmax = opts.get('rmax',1.01)
    step = opts.get('step',0.05)
    doplot = opts.get('doplot',False)
    rs = arange(rmin,rmax,step)
    E = []
    for r in rs:
        en = main(r=r,**opts)
        print "%.3f %.5f" % (r,en)
        E.append(en)

    if doplot:
        from pylab import plot
        plot(rs,E,'go-')
    return

def main(**opts):
    r = opts.get('r',0.70)
    h2 = Molecule('h2',
                  atomlist = [(1,(0,0,r/2.)),
                              (1,(0,0,-r/2.))],
                  units='Angs')
    en,orbe,orbs = dft(h2,**opts)
    return en

if __name__ == '__main__': main()
    
                                   
    
