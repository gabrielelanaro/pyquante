#!/usr/bin/env python

from PyQuante.NumWrap import arange
from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

energy = -7.813809
name = "LiH/DFT"

def dissoc(**opts):
    rmin = opts.get('rmin',0.65)
    rmax = opts.get('rmax',1.01)
    step = opts.get('step',0.05)
    doplot = opts.get('doplot',False)
    rs = arange(rmin,rmax,step)
    E = []
    for r in rs:
        lih = Molecule('lih',
                       atomlist = [(3,(0,0,0)),
                                   (1,(0,0,r))],
                       units='Angs')
        en,orbe,orbs = dft(lih,**opts)
        print "%.3f %.5f" % (r,en)
        E.append(en)

    if doplot:
        from pylab import plot
        plot(rs,E,'bo-')
    return

def main(**opts):
    r = opts.get('r',1.0)
    lih = Molecule('lih',
                   atomlist = [(3,(0,0,0)),
                               (1,(0,0,r))],
                   units='Angs')
    en,orbe,orbs = dft(lih,**opts)

    return en

def profmain():
    import profile,pstats
    en = 0
    profile.run('main()','prof')
    pstats.Stats('prof').strip_dirs().sort_stats('time').print_stats(15)
    return

if __name__ == '__main__':
    main()
    #profmain()
    
                                   
    
