"""\
 pyints.py  Python implementations of work functions for Gaussian
  integrals in the PyQuante package.

 The equations herein are based upon
 'Gaussian Expansion Methods for Molecular Orbitals.' H. Taketa,
 S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
 [THO paper].

 The accompanying module cints contains the same functions written in C
 to be faster.
   
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from math import sqrt,pow,log,exp,pi,floor

def fact(i):
    "Normal factorial"
    val = 1
    while (i>1):
        val = i*val
        i = i-1
    return val

def fact2(i):
    "Double factorial (!!) function = 1*3*5*...*i"
    val = 1
    while (i>0):
        val = i*val
        i = i-2
    return val

def dist2(A,B):
    return pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2)

def dist(A,B): return sqrt(dist2(A,B))

def binomial_prefactor(s,ia,ib,xpa,xpb):
    "From Augspurger and Dykstra"
    sum = 0
    for t in xrange(s+1):
        if s-ia <= t <= ib:
            sum = sum + binomial(ia,s-t)*binomial(ib,t)* \
                  pow(xpa,ia-s+t)*pow(xpb,ib-t)
    return sum

def binomial(a,b):
    "Binomial coefficient"
    return fact(a)/fact(b)/fact(a-b)

def Fgamma(m,x):
    "Incomplete gamma function"
    SMALL=0.00000001
    x = max(abs(x),SMALL)
    val = gamm_inc(m+0.5,x)
    return 0.5*pow(x,-m-0.5)*val;

def gammln(x):
    "Numerical recipes, section 6.1"
    cof = [76.18009172947146,-86.50532032941677,
           24.01409824083091,-1.231739572450155,
           0.1208650973866179e-2,-0.5395239384953e-5]
    y=x
    tmp=x+5.5
    tmp = tmp - (x+0.5)*log(tmp)
    ser=1.000000000190015 # don't you just love these numbers?!
    for j in xrange(6):
        y = y+1
        ser = ser+cof[j]/y
    return -tmp+log(2.5066282746310005*ser/x);

def gamm_inc(a,x):
    "Incomple gamma function \gamma; computed from NumRec routine gammp."
    gammap,gln = gammp(a,x)
    return exp(gln)*gammap
    
def gammp(a,x):
    "Returns the incomplete gamma function P(a;x). NumRec sect 6.2."
    assert (x > 0 and a >= 0), "Invalid arguments in routine gammp"

    if x < (a+1.0): #Use the series representation
        gamser,gln = _gser(a,x)
        return gamser,gln
    #Use the continued fraction representation
    gammcf,gln = _gcf(a,x)
    return 1.0-gammcf ,gln  #and take its complement.

def _gser(a,x):
    "Series representation of Gamma. NumRec sect 6.1."
    ITMAX=100
    EPS=3.e-7

    gln=gammln(a)
    assert(x>=0),'x < 0 in gser'
    if x == 0 : return 0,gln

    ap = a
    delt = sum = 1./a
    for i in xrange(ITMAX):
        ap=ap+1.
        delt=delt*x/ap
        sum=sum+delt
        if abs(delt) < abs(sum)*EPS: break
    else:
        print 'a too large, ITMAX too small in gser'
    gamser=sum*exp(-x+a*log(x)-gln)
    return gamser,gln

def _gcf(a,x):
    "Continued fraction representation of Gamma. NumRec sect 6.1"
    ITMAX=100
    EPS=3.e-7
    FPMIN=1.e-30

    gln=gammln(a)
    b=x+1.-a
    c=1./FPMIN
    d=1./b
    h=d
    for i in xrange(1,ITMAX+1):
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if abs(d) < FPMIN: d=FPMIN
        c=b+an/c
        if abs(c) < FPMIN: c=FPMIN
        d=1./d
        delt=d*c
        h=h*delt
        if abs(delt-1.) < EPS: break
    else:
        print 'a too large, ITMAX too small in gcf'
    gammcf=exp(-x+a*log(x)-gln)*h
    return gammcf,gln
    
    
def erf(z): 
    if z==0.0: return 0.0
    else: return gammp(0.5,z**2)[0]

def ijkl2intindex(i,j,k,l):
    "Indexing into the get2ints long array"
    if i<j: i,j = j,i
    if k<l: k,l = l,k
    ij = i*(i+1)/2+j
    kl = k*(k+1)/2+l
    if ij < kl: ij,kl = kl,ij
    return ij*(ij+1)/2+kl

def fB(i,l1,l2,P,A,B,r,g):
    return binomial_prefactor(i,l1,l2,P-A,P-B)*B0(i,r,g)

def B0(i,r,g): return fact_ratio2(i,r)*pow(4*g,r-i)

def fact_ratio2(a,b): return fact(a)/fact(b)/fact(a-2*b)

def kinetic(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B):
    term0 = alpha2*(2*(l2+m2+n2)+3)*\
            overlap(alpha1,(l1,m1,n1),A,\
                           alpha2,(l2,m2,n2),B)
    term1 = -2*pow(alpha2,2)*\
            (overlap(alpha1,(l1,m1,n1),A,
                            alpha2,(l2+2,m2,n2),B)
             + overlap(alpha1,(l1,m1,n1),A,
                              alpha2,(l2,m2+2,n2),B)
             + overlap(alpha1,(l1,m1,n1),A,
                              alpha2,(l2,m2,n2+2),B))
    term2 = -0.5*(l2*(l2-1)*overlap(alpha1,(l1,m1,n1),A,
                                           alpha2,(l2-2,m2,n2),B) +
                  m2*(m2-1)*overlap(alpha1,(l1,m1,n1),A,
                                           alpha2,(l2,m2-2,n2),B) +
                  n2*(n2-1)*overlap(alpha1,(l1,m1,n1),A,
                                           alpha2,(l2,m2,n2-2),B))
    return term0+term1+term2
    

def overlap(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B):
    "Taken from THO eq. 2.12"
    rab2 = dist2(A,B)
    gamma = alpha1+alpha2
    P = gaussian_product_center(alpha1,A,alpha2,B)

    pre = pow(pi/gamma,1.5)*exp(-alpha1*alpha2*rab2/gamma)

    wx = overlap_1D(l1,l2,P[0]-A[0],P[0]-B[0],gamma)
    wy = overlap_1D(m1,m2,P[1]-A[1],P[1]-B[1],gamma)
    wz = overlap_1D(n1,n2,P[2]-A[2],P[2]-B[2],gamma)
    return pre*wx*wy*wz

def multipole_ints((kx,ky,kz),alpha1,(l1,m1,n1),A1,alpha2,(l2,m2,n2),A2):
    total = 0
    for ix in xrange(kx+1):
        for iy in xrange(ky+1):
            for iz in xrange(kz+1):
                total += binomial(kx,ix)*binomial(ky,iy)*binomial(kz,iz)*\
                         A1[0]**(kx-ix)*A1[1]**(ky-iy)*A1[2]**(kz-iz)*\
                         overlap(alpha1,(l1+ix,m1+iy,n1+iz),A1,
                                 alpha2,(l2,m2,n2),A2)
    return total    

def overlap_1D(l1,l2,PAx,PBx,gamma):
    "Taken from THO eq. 2.12"
    sum = 0
    for i in xrange(1+int(floor(0.5*(l1+l2)))):
        sum = sum + binomial_prefactor(2*i,l1,l2,PAx,PBx)* \
              fact2(2*i-1)/pow(2*gamma,i)
    return sum

    
def gaussian_product_center(alpha1,A,alpha2,B):
    gamma = alpha1+alpha2
    return (alpha1*A[0]+alpha2*B[0])/gamma,\
           (alpha1*A[1]+alpha2*B[1])/gamma,\
           (alpha1*A[2]+alpha2*B[2])/gamma

def nuclear_attraction((x1,y1,z1),norm1,(l1,m1,n1),alpha1,
                       (x2,y2,z2),norm2,(l2,m2,n2),alpha2,
                       (x3,y3,z3)):

        gamma = alpha1+alpha2

        xp,yp,zp = gaussian_product_center(alpha1,(x1,y1,z1),alpha2,(x2,y2,z2))
        rab2 = dist2((x1,y1,z1),(x2,y2,z2))
        rcp2 = dist2((x3,y3,z3),(xp,yp,zp))

        Ax = A_array(l1,l2,xp-x1,xp-x2,xp-x3,gamma)
        Ay = A_array(m1,m2,yp-y1,yp-y2,yp-y3,gamma)
        Az = A_array(n1,n2,zp-z1,zp-z2,zp-z3,gamma)

        sum = 0.
        for I in xrange(l1+l2+1):
            for J in xrange(m1+m2+1):
                for K in xrange(n1+n2+1):
                    sum = sum + Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K,rcp2*gamma)

        return -norm1*norm2*\
               2*pi/gamma*exp(-alpha1*alpha2*rab2/gamma)*sum
    
def A_term(i,r,u,l1,l2,PAx,PBx,CPx,gamma):
    "THO eq. 2.18"
    return pow(-1,i)*binomial_prefactor(i,l1,l2,PAx,PBx)*\
           pow(-1,u)*fact(i)*pow(CPx,i-2*r-2*u)*\
           pow(0.25/gamma,r+u)/fact(r)/fact(u)/fact(i-2*r-2*u)

def A_array(l1,l2,PA,PB,CP,g):
    "THO eq. 2.18 and 3.1"
    Imax = l1+l2+1
    A = [0]*Imax
    for i in xrange(Imax):
        for r in xrange(int(floor(i/2)+1)):
            for u in xrange(int(floor((i-2*r)/2)+1)):
                I = i-2*r-u
                A[I] = A[I] + A_term(i,r,u,l1,l2,PA,PB,CP,g)
    return A

def grad_A_array(l1,l2,PA,PB,CP,g):
    "several pages of algebra were necessary to find this result..."
    Imax = l1+l2+1
    A = [0]*Imax
    for i in xrange(Imax):
        for r in xrange(int(floor(i/2)+1)):
            for u in xrange(int(floor((i-2*r)/2.)+1)):
                I = i-2*r-u
                A[I] = A[I] + (i-2.*r+1.)/(i-2.*r-2.*u+1.)*CP*A_term(i,r,u,l1,l2,PA,PB,CP,g)
    return A
               

def grad_nuc_att((x1,y1,z1),(l1,m1,n1),alpha1,
                 (x2,y2,z2),(l2,m2,n2),alpha2,
                 (x3,y3,z3)):
                 
    gamma = alpha1+alpha2

    xp,yp,zp = gaussian_product_center(alpha1,(x1,y1,z1),alpha2,(x2,y2,z2))
    rab2 = dist2((x1,y1,z1),(x2,y2,z2))
    rcp2 = dist2((x3,y3,z3),(xp,yp,zp))

    Ax = A_array(l1,l2,xp-x1,xp-x2,xp-x3,gamma)
    Ay = A_array(m1,m2,yp-y1,yp-y2,yp-y3,gamma)
    Az = A_array(n1,n2,zp-z1,zp-z2,zp-z3,gamma)

    grad_Ax = grad_A_array(l1,l2,xp-x1,xp-x2,xp-x3,gamma)
    grad_Ay = grad_A_array(m1,m2,yp-y1,yp-y2,yp-y3,gamma)
    grad_Az = grad_A_array(n1,n2,zp-z1,zp-z2,zp-z3,gamma)
    
    sum_x = 0.
    sum_y = 0.
    sum_z = 0.
    for I in xrange(l1+l2+1):
        for J in xrange(m1+m2+1):
            for K in xrange(n1+n2+1):
                sum_x += grad_Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K+1,rcp2*gamma)
                sum_y += Ax[I]*grad_Ay[J]*Az[K]*Fgamma(I+J+K+1,rcp2*gamma)
                sum_z += Ax[I]*Ay[J]*grad_Az[K]*Fgamma(I+J+K+1,rcp2*gamma)

    gradV_x = -4*pi*exp(-alpha1*alpha2*rab2/gamma)*sum_x
    gradV_y = -4*pi*exp(-alpha1*alpha2*rab2/gamma)*sum_y
    gradV_z = -4*pi*exp(-alpha1*alpha2*rab2/gamma)*sum_z
        
    return gradV_x,gradV_y,gradV_z
    
    
def contr_coulomb(aexps,acoefs,anorms,xyza,powa,
                  bexps,bcoefs,bnorms,xyzb,powb,
                  cexps,ccoefs,cnorms,xyzc,powc,
                  dexps,dcoefs,dnorms,xyzd,powd):

    Jij = 0.
    for i in xrange(len(aexps)):
        for j in xrange(len(bexps)):
            for k in xrange(len(cexps)):
                for l in xrange(len(dexps)):
                    incr = coulomb_repulsion(xyza,anorms[i],powa,aexps[i],
                                             xyzb,bnorms[j],powb,bexps[j],
                                             xyzc,cnorms[k],powc,cexps[k],
                                             xyzd,dnorms[l],powd,dexps[l])
                    Jij = Jij + acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]*incr
    return Jij
    

def coulomb_repulsion((xa,ya,za),norma,(la,ma,na),alphaa,
                      (xb,yb,zb),normb,(lb,mb,nb),alphab,
                      (xc,yc,zc),normc,(lc,mc,nc),alphac,
                      (xd,yd,zd),normd,(ld,md,nd),alphad):

    rab2 = dist2((xa,ya,za),(xb,yb,zb))
    rcd2 = dist2((xc,yc,zc),(xd,yd,zd))
    xp,yp,zp = gaussian_product_center(alphaa,(xa,ya,za),alphab,(xb,yb,zb))
    xq,yq,zq = gaussian_product_center(alphac,(xc,yc,zc),alphad,(xd,yd,zd))
    rpq2 = dist2((xp,yp,zp),(xq,yq,zq))
    gamma1 = alphaa+alphab
    gamma2 = alphac+alphad
    delta = 0.25*(1/gamma1+1/gamma2)

    Bx = B_array(la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta)
    By = B_array(ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta)
    Bz = B_array(na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta)

    sum = 0.
    for I in xrange(la+lb+lc+ld+1):
        for J in xrange(ma+mb+mc+md+1):
            for K in xrange(na+nb+nc+nd+1):
                sum = sum + Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,0.25*rpq2/delta)

    return 2*pow(pi,2.5)/(gamma1*gamma2*sqrt(gamma1+gamma2)) \
           *exp(-alphaa*alphab*rab2/gamma1) \
           *exp(-alphac*alphad*rcd2/gamma2)*sum*norma*normb*normc*normd

def B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,Px,Ax,Bx,Qx,Cx,Dx,gamma1,gamma2,delta):
    "THO eq. 2.22"
    return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1) \
           *pow(-1,i2)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2) \
           *pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u) \
           *pow(Qx-Px,i1+i2-2*(r1+r2)-2*u) \
           /pow(delta,i1+i2-2*(r1+r2)-u)

def B_array(l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta):
    Imax = l1+l2+l3+l4+1
    B = [0]*Imax
    for i1 in xrange(l1+l2+1):
        for i2 in xrange(l3+l4+1):
            for r1 in xrange(i1/2+1):
                for r2 in xrange(i2/2+1):
                    for u in xrange((i1+i2)/2-r1-r2+1):
                        I = i1+i2-2*(r1+r2)-u
                        B[I] = B[I] + B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
                                             p,a,b,q,c,d,g1,g2,delta)
    return B


def three_center_1D(xi,ai,alphai,xj,aj,alphaj,xk,ak,alphak):
    gamma = alphai+alphaj+alphak
    dx = exp(-alphai*alphaj*(xi-xj)**2/gamma) * \
         exp(-alphai*alphak*(xi-xk)**2/gamma) * \
         exp(-alphaj*alphak*(xj-xk)**2/gamma)
    px = (alphai*xi+alphaj*xj+alphak*xk)/gamma
    
    xpi = px-xi
    xpj = px-xj
    xpk = px-xk
    int = 0
    for q in xrange(ai+1):
        for r in xrange(aj+1):
            for s in xrange(ak+1):
                n,rem = divmod(q+r+s,2)
                if (q+r+s) %2 == 0:
                    i_qrs = binomial(ai,q)*binomial(aj,r)*binomial(ak,s)*\
                            pow(xpi,ai-q)*pow(xpj,aj-r)*pow(xpk,ak-s)*\
                            fact2(2*n-1)/pow(2*gamma,n)*sqrt(pi/gamma)
                    int += i_qrs
    return dx*int

if __name__ == '__main__':
    print multipole_ints((0,0,1),1.0,(0,0,0),(0.,0.,0.),1.0,(0,0,0),(0.,0.,1.))
    
