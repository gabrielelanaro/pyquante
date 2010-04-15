"""\
 Slater.py: Coefficients for fitting Gaussian functions to Slater
   Functions.These functions are STO-6G fits to Slater exponents
   with exponents of 1. To fit to exponents of \zeta, one need only
   multiply each exponent by \zeta^2. For STO-1G, STO-2G and other
   levels of fit, see the paper.

 Reference: RF Stewart, JCP 52, 431 (1970)

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from CGBF import CGBF
from MINDO3_Parameters import zetas,zetap
from PyQuante.cints import fact
from math import pi,sqrt,exp
from NumWrap import array


gauss_powers = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]

# Gaussian functions for fitting to Slaters. These functions are
# STO-6G fits to slater exponents with exponents of 1. To fit
# to exponents of \zeta, you need only multiply each
# exponent by \zeta^2
# The rest of these functions can be obtained from Stewart,
#  JCP 52, 431 (1970)

gexps_1s = [2.310303149e01,4.235915534e00,1.185056519e00,
            4.070988982e-01,1.580884151e-01,6.510953954e-02]
gcoefs_1s = [9.163596280e-03,4.936149294e-02,1.685383049e-01,
             3.705627997e-01,4.164915298e-01,1.303340841e-01]

gexps_2s = [2.768496241e01,5.077140627e00,1.426786050e00,
            2.040335729e-01,9.260298399e-02,4.416183978e-02]
gcoefs_2s = [-4.151277819e-03,-2.067024148e-02,-5.150303337e-02,
             3.346271174e-01,5.621061301e-01,1.712994697e-01]

gexps_2p = [5.868285913e00,1.530329631e00,5.475665231e-01,
            2.288932733e-01,1.046655969e-01,4.948220127e-02]
gcoefs_2p = [7.924233646e-03,5.144104825e-02,1.898400060e-01,
             4.049863191e-01,4.012362861e-01,1.051855189e-01]

gexps_3s = [3.273031938e00,9.200611311e-01,3.593349765e-01,
            8.636686991e-02,4.797373812e-02,2.724741144e-02]
gcoefs_3s = [-6.775596947e-03,-5.639325779e-02,-1.587856086e-01,
             5.534527651e-01,5.015351020e-01,7.223633674e-02]

gexps_3p = [5.077973607e00,1.340786940e00,2.248434849e-01,
            1.131741848e-01,6.076408893e-02,3.315424265e-02]
gcoefs_3p = [-3.329929840e-03,-1.419488340e-02,1.639395770e-01,
             4.485358256e-01,3.908813050e-01,7.411456232e-02]
gexps_3d = [2.488296923,7.981487853e-1,3.311327490e-1,
            1.559114463e-1,7.877734732e-2,4.058484363e-2]
gcoefs_3d = [7.283828112e-3,5.386799363e-2,2.072139149e-1,
             4.266269092e-1,3.843100204e-1,8.902827546e-2]

gexps_4s = [3.232838646,3.605788802e-1,1.717902487e-1,
            5.277666487e-2,3.163400284e-2,1.874093091e-2]
gcoefs_4s = [1.374817488e-3,-8.666390043e-2,-3.130627309e-1,
             7.812787397e-1,4.389247988-1,2.487178756e-2]
gexps_4p = [2.389722618, 7.960947826e-1,3.415541380e-1,
            8.847434525e-2,4.958248334e-2,2.816929784e-2]
gcoefs_4p = [-1.665913575e-3,-1.657464971e-2,-5.958513378e-2,
             4.053115554e-1,5.433958189e-1,1.20970491e-1]

gfitfuncs = {'1S':(gexps_1s,gcoefs_1s),
             '2S':(gexps_2s,gcoefs_2s),
             '2P':(gexps_2p,gcoefs_2p),
             '3S':(gexps_3s,gcoefs_3s),
             '3P':(gexps_3p,gcoefs_3p),
             '3D':(gexps_3d,gcoefs_3d),
             '4S':(gexps_4s,gcoefs_4s),
             '4P':(gexps_4p,gcoefs_4p)}

# Here are the STO-6G values from Hehre, Stewart, Pople JCP 51,2657 (1969),
# and Hehre, Ditchfield, Stewart, Pople JCP 52, 2769 (1970)
# which are a little different, in that they use the same exponent for
# 2s,2p, and 3s,3p, which makes the fit a bit different.
gexps_old_2 = [1.03087e1,2.04036,6.34142e-1,
               2.43977e-1,1.05960e-1,4.85690e-2]
gcoefs_old_2s = [-1.32528e-2,-4.69917e-2,-3.37854e-2,
                 2.50242e-1,2.95117e-1,2.40706e-1]
gcoefs_old_2p = [3.75970e-3,3.76794e-2,1.73897e-1,
                 4.18036e-1,4.25860e-1,1.017008e-1]
gexps_old_3 = [3.0817,8.24896e-1,3.09345e-1,
               1.38468e-1,6.85210e-2,3.53133e-2]
gcoefs_old_3s = [-7.94313e-3,-7.10026e-2,-1.78503e-1,
                 1.51064e-1,7.35491e-1,2.76059e-1]
gcoefs_old_3p = [-7.13936e-3,-1.82928e-2,7.62162e-2,
                 4.14510e-1,4.88962e-1,1.05882e-1]

NQN = [ None, 1, 1, # principle quantum number N
        2, 2, 2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3]

gexps = { # indexed by N,s_or_p:
    (1,0) : gexps_1s,
    (2,0) : gexps_2s,
    (2,1) : gexps_2p,
    (3,0) : gexps_3s,
    (3,1) : gexps_3p
    }

gcoefs = {  # indexed by N,s_or_p:
    (1,0) : gcoefs_1s,
    (2,0) : gcoefs_2s,
    (2,1) : gcoefs_2p,
    (3,0) : gcoefs_3s,
    (3,1) : gcoefs_3p
    }

gexps_old = { # indexed by N,s_or_p:
    (1,0) : gexps_1s,
    (2,0) : gexps_old_2,
    (2,1) : gexps_old_2,
    (3,0) : gexps_old_3,
    (3,1) : gexps_old_3
    }

gcoefs_old = {  # indexed by N,s_or_p:
    (1,0) : gcoefs_1s,
    (2,0) : gcoefs_old_2s,
    (2,1) : gcoefs_old_2p,
    (3,0) : gcoefs_3s,
    (3,1) : gcoefs_3p
    }

s_or_p = [0,1,1,1] # whether the func is s or p type, based on the L QN

# This is a translation of the mopac gover function, which uses an
# STO-6G expansion of slater functions and then computes the Gaussian
# overlap. The advantage of this function over the others is that
# derivatives are easier with Gaussians.
def get_overlap(atnoi,i,xyzi,atnoj,j,xyzj):
    bohr2ang = 0.52918
    xyzi = (xyzi[0]/bohr2ang,xyzi[1]/bohr2ang,xyzi[2]/bohr2ang)
    gi = CGBF(xyzi,gauss_powers[i])
    
    zi = gexps[(NQN[atnoi],s_or_p[i])]
    ci = gcoefs[(NQN[atnoi],s_or_p[i])]
    if i:
        zetai = zetap[atnoi]
    else:
        zetai = zetas[atnoi]
    
    xyzj = (xyzj[0]/bohr2ang,xyzj[1]/bohr2ang,xyzj[2]/bohr2ang)
    gj = CGBF(xyzj,gauss_powers[j])
    zj = gexps[(NQN[atnoj],s_or_p[j])]
    cj = gcoefs[(NQN[atnoj],s_or_p[j])]
    if j:
        zetaj = zetap[atnoj]
    else:
        zetaj = zetas[atnoj]

    # Multiply the functions by \zeta^2
    for a in xrange(6):
        gi.add_primitive(zi[a]*zetai*zetai,ci[a])
        gj.add_primitive(zj[a]*zetaj*zetaj,cj[a])
    gi.normalize()
    gj.normalize()
    return gi.overlap(gj)


class SlaterFunction:
    """Cartesian Slater basis functions"""
    nlm2powers = {
        (1,0,0) : (0,0,0,0),   # x,y,z,r
        (2,0,0) : (0,0,0,1),
        (3,0,0) : (0,0,0,2),
        (2,1,0) : (1,0,0,0),
        (2,1,1) : (0,1,0,0),
        (2,1,-1) : (0,0,1,0),
        (3,1,0) : (1,0,0,1),
        (3,1,1) : (0,1,0,1),
        (3,1,-1) : (0,0,1,1)
        }

    def __init__(self,N,L,M,origin,zeta):
        self.N = N # 1 2 3, etc.
        self.L = L
        self.M = M # [x,y,z] -> -1 -> z, 0->x
        self.origin = origin
        self.zeta = zeta
        self.normalize()
        return

    def normalize(self):
        N = self.N
        zeta = self.zeta
        s_norm = sqrt(pow(zeta,2*N+1)*pow(2,2*N-1)/pi/fact(2*N))
        if self.L == 0:
            self.norm = s_norm
        elif self.L == 1:
            self.norm = sqrt(3)*s_norm
        else:
            raise Exception("d,f,... Slater functions not implemented yet")
        return

    def amp(self,xyz):
        x = xyz[0]-self.origin[0]
        y = xyz[1]-self.origin[1]
        z = xyz[2]-self.origin[2]
        r = sqrt(x*x+y*y+z*z)
        a,b,c,d = self.nlm2powers[(self.N,self.L,self.M)]
        return pow(x,a)*pow(y,b)*pow(z,c)*pow(r,d)*\
               self.norm*exp(-self.zeta*r) 
    
    def laplacian(self,xyz):
        return self.laplacian_ratio(xyz)*self.amp(xyz)

    def laplacian_ratio(self,xyz):
        x2 = pow(xyz[0]-self.origin[0],2)
        y2 = pow(xyz[1]-self.origin[1],2)
        z2 = pow(xyz[2]-self.origin[2],2)
        r2 = x2+y2+z2
        r = sqrt(r2)
        a,b,c,d = self.nlm2powers[(self.N,self.L,self.M)]
        zeta = self.zeta
        return d*(2*a+2*b+2*c+2*d+1)/r2 - 2*zeta*(a+b+c+d+1)/r \
               + zeta*zeta + a*(a-1)/x2 + b*(b-1)/y2 + c*(c-1)/z2

    def grad(self,xyz):
        amp = self.amp(xyz)
        gx,gy,gz = self.grad_ratio(xyz)
        return array((amp*gx,amp*gy,amp*gz))

    def grad_ratio(self,xyz):
        x = xyz[0]-self.origin[0]
        y = xyz[1]-self.origin[1]
        z = xyz[2]-self.origin[2]
        r2 = x*x+y*y+z*z
        r = sqrt(r2)
        a,b,c,d = self.nlm2powers[(self.N,self.L,self.M)]
        zeta = self.zeta
        return (a/x + d*x/r2 - x*zeta/r,
                b/y + d*y/r2 - y*zeta/r,
                c/z + d*z/r2 - z*zeta/r)

    def gaussian(self):
        a,b,c,d = self.nlm2powers[(self.N,self.L,self.M)]
        cgbf = CGBF(self.origin,(a,b,c))
        expos,coefs = gexps[(self.N,self.L)],gcoefs[(self.N,self.L)]
        zet2 = self.zeta*self.zeta
        for i in xrange(6):
            cgbf.add_primitive(zet2*expos[i],coefs[i])
        cgbf.normalize()
        return cgbf        

if __name__ == '__main__':
    print get_overlap(8,0,(0,0,0),1,0,(1.,0,0)),0.2521088
    print get_overlap(8,1,(0,0,0),1,0,(1.,0,0)),0.3994043


