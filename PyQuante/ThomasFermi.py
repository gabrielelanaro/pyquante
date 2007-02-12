#!/usr/bin/env python
"""\
 Thomas-Fermi theory for atoms

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from math import sqrt,pi,exp
from Gnuplot import Gnuplot,Data
from NumWrap import array,zeros
from Pistol.Util import frange

# Chi is the dimensionless version of the potential

# I think these are originally from Landau/Lifshitz, but I took them
#  from Flugge
chi_flugge = [ # x, Chi(x)
    (0.,1.0), (0.02,0.9720), (0.04,0.9470), (0.06,0.9238), (0.08,0.9022),
    (0.10,0.8817), (0.2,0.7931), (0.3,0.7206), (0.4,0.6595), (0.5,0.6070),
    (0.6,0.5612), (0.7,0.5208), (0.8,0.4849), (0.9,0.4529), (1.0,0.4240),
    (1.2,0.3742), (1.4,0.3329), (1.6,0.2981), (1.8,0.2685), (2.0,0.2430),
    (2.2,0.2210), (2.4,0.2017), (2.6,0.1848), (2.8,0.1699), (3.0,0.1566),
    (3.2,0.1448), (3.4,0.1343), (3.6,0.1247), (3.8,0.1162), (4.0,0.1084),
    (4.5,0.0919), (5.0,0.0788), (5.5,0.0682), (6.0,0.0594), (6.5,0.0522),
    (7.0,0.0461), (7.5,0.0410), (8.0,0.0366), (8.5,0.0328), (9.0,0.0296),
    (9.5,0.0268), (10.0,0.0243)]

def chi_sommer(x):
    # Sommerfield's fit to Chi
    a = pow(12.,2./3.)
    d = 0.772
    c = 3/d
    return pow(1+pow(x/a,d),-c)

# Tietz's approximation of Chi
def chi_tietz(x): return pow(1+0.53625*x,-2)
    
def energy(Z): return -0.7687*pow(Z,7./3.)
def energy_eV(Z): return -20.93*pow(Z,7./3.)

def rho(r,Z):
    alpha = 1.1295*pow(Z,1./3.)
    x = alpha*r
    #chi = chi_tietz(x)
    chi = chi_sommer(x)
    V = -Z*chi/r
    return pow(-2*V,1.5)/3./pi/pi

def chi_integrate(xend):
    # This doesn't quite work yet, but you can get the general idea:
    #  use a Verlet scheme to integrate to +inf
    x = 0.00000001
    f = 1.0
    dx = 0.0001
    df = -1.58807
    xs = []
    fs = []
    while x < xend:
        xs.append(x)
        fs.append(f)        
        x += dx
        if f < 0: break
        ddf = pow(f,1.5)/sqrt(x)
        f += dx*df + 0.5*dx*dx*ddf
        df += dx*ddf
    return xs,fs

def chi_rk2(tend):
    # Same thing as chi_integrate, but uses a second-order Runge-Kutta
    #  scheme
    t = 0.001
    X = 1.0
    h = 0.00001
    v = -1.58807
    ts = []
    Xs = []
    while t < 3:
        ts.append(t)
        Xs.append(X)
        X_ = X + 0.5*h*v
        v_ = v + 0.5*h*pow(X,1.5)/sqrt(t)
        X = X + h*v_
        v = v + h*pow(X_,1.5)/sqrt(t+0.5*h)
        t += h
        if X < 0: break
    return ts,Xs

def F((x,v),t):
    if t < 1e-5: return array((v,0))
    return array((v,pow(x,1.5)/sqrt(t)))

def chi_rk4(tend):
    # Same thing as chi_integrate, but uses a fourth-order Runge-Kutta
    #  scheme
    t = 0
    X = 1.0
    h = 1e-4 # 1e-5 is not small enough and still takes forever!
    v = -1.58807
    ts = []
    Xs = []
    while t < 1:
        ts.append(t)
        Xs.append(X)
        A = array((X,v))
        k1 = h*F(A,t)
        k2 = h*F(A+k1/2,t+h/2)
        k3 = h*F(A+k2/2,t+h/2)
        k4 = h*F(A+k3,t+h)
        A += (k1+2*k2+2*k3+k4)/6
        X,v = A
        t += h
        if X < 0: break
    return ts,Xs

def test():
    xs = []
    chi1 = []
    chi2 = []
    chi4 = []
    chi5 = []
    for x,chi in chi_flugge:
        xs.append(x)
        chi1.append(chi)
        chi2.append(chi_sommer(x))
        chi4.append(chi_tietz(x))
    #xs3,chi3 = chi_rk4(10)
    g = Gnuplot()
    #g.title('Comparison of Thomas-Fermi Chi values')
    d1 = Data(xs,chi1,title='Analytic',with='points')
    #d2 = Data(xs,chi2,title='Sommerfield',with='linespoints')
    #d3 = Data(xs3,chi3,title="Integrated",with='lines')
    d3 = Data(xs,chi4,title='Tietz',with='lines')
    g.ylabel('Chi(x)')
    g.xlabel('x')
    g.plot(d1,d3)
    g("set term post")
    g("set output 'tietz.eps'")
    g.replot()
    raw_input('Press any key to continue')
    
if __name__ == '__main__': test()
