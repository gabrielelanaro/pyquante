# mobihf v0.1, Quantum Chemistry on mobile!, (c) V.Ganesh, GPL
# Based on quantumj (http://quantumj.dev.java.net)
#       and PyQuante (http://pyquante.sf.net)
from math import sin,cos,acos,pi,sqrt,floor,pow,exp,log
from string import *

radii = {"H":0.23, "C":0.77, "N":0.75, "O":0.73}
atmnum = {"H":1, "C":6, "N":7, "O":8}
    
basis_sto = {
    'H': [('S',[
           (3.425251, 0.154329),
           (0.623914, 0.535328),
           (0.168855, 0.444635)])],   
    'C': [('S',[
           (71.616837, 0.154329),
           (13.045096, 0.535328),
           (3.530512, 0.444635)]),
         ('S',[
           (2.941249, -0.099967),
           (0.683483, 0.399513),
           (0.222290, 0.700115)]),
         ('P',[
           (2.941249, 0.155916),
           (0.683483, 0.607684),
           (0.222290, 0.391957)])],
    'N': [('S',[
           (99.106169, 0.154329),
           (18.052312, 0.535328),
           (4.885660, 0.444635)]),
         ('S',[
           (3.780456, -0.099967),
           (0.878497, 0.399513),
           (0.285714, 0.700115)]),
         ('P',[
           (3.780456, 0.155916),
           (0.878497, 0.607684),
           (0.285714, 0.391957)])],
    'O': [('S',[
           (130.709321, 0.154329),
           (23.808866, 0.535328),
           (6.443608, 0.444635)]),
         ('S',[
           (5.033151, -0.099967),
           (1.169596, 0.399513),
           (0.380389, 0.700115)]),
         ('P',[
           (5.033151, 0.155916),
           (1.169596, 0.607684),
           (0.380389, 0.391957)])]
    }

curbas = basis_sto

powers = {
    'S': [(0, 0, 0)],
    'P': [(1, 0, 0), (0, 1, 0), (0, 0, 1)]    
    }

class molecule:
    def __init__(self):
        self.noOfAtoms = 0
        self.title = ""
        self.atoms = []

class Integral:
    def __init__(self):
        self.MAX_ITERATION = 100
        self.EPS = 3.0e-7
        self.SMALL = 0.00000001
        self.FPMIN = 1.0e-30
        self.cof = [76.18009172947146, -86.50532032941677, \
        24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5]
    
    def factorial(self, n):
        value = 1        
        while(n > 1): value = value * n; n-=1
        return value
    
    def factorial2(self, n):
        value = 1        
        while(n > 0): value = value * n; n-=2
        return value
    
    def factorialRatioSquared(self, a, b):
        return (self.factorial(a) / self.factorial(b) / self.factorial(a-2*b))
    
    def binomial(self, i, j):
        return (self.factorial(i) / self.factorial(j) / self.factorial(i-j))

    def mpow(self, x, y):
        if (round(y,10)==0.0): return 1.0
        if (round(x,10)==0.0): return 0.0
        return pow(x, y)
        
    def binomialPrefactor(self, s, ia, ib, xpa, xpb):
        sum = 0.0        
        for t in range(0, s+1):
            if (((s-ia) <= t) and (t <= ib)):
                sum += self.binomial(ia, s-t) * self.binomial(ib, t) \
                        * self.mpow(xpa, ia-s+t) * self.mpow(xpb, ib-t)
        return sum
    
    def dist2(self, a, b):
        x=a[0]-b[0]; y=a[1]-b[1]; z=a[2]-b[2]
        return (x*x+y*y+z*z)

    def dist(self, a, b):
        return (sqrt(self.dist2(a, b)))
    
    def ovrlp1D(self, l1, l2, pax, pbx, gamma):
        sum = 0.0
        k = 1 + int(floor(0.5*(l1+l2)));        
        for i in range(0, k):
            sum += (self.binomialPrefactor(2*i, l1, l2, pax, pbx) \
                    * self.factorial2(2*i-1)) / self.mpow(2*gamma, i)
        return sum

    def gaussianProductCenter(self, alp1, a, alp2, b):
        gamma = alp1 + alp2
        return ((alp1*a[0]+alp2*b[0])/gamma, \
                (alp1*a[1]+alp2*b[1])/gamma, \
                (alp1*a[2]+alp2*b[2])/gamma)
        
    def ovrlp(self, alp1, pwr, a, alp2, pwr2, b):
        rab2 = self.dist2(a, b)
        gamma = alp1 + alp2
        prd = self.gaussianProductCenter(alp1, a, alp2, b)

        wx = self.ovrlp1D(pwr[0], pwr2[0], prd[0]-a[0], prd[0]-b[0], gamma)
        wy = self.ovrlp1D(pwr[1], pwr2[1], prd[1]-a[1], prd[1]-b[1], gamma)
        wz = self.ovrlp1D(pwr[2], pwr2[2], prd[2]-a[2], prd[2]-b[2], gamma)
        
        return (pow(pi/gamma, 1.5)*exp((-alp1*alp2*rab2)/gamma)*wx*wy*wz)

    def conATrm(self, i, r, u, l1, l2, pax, pbx, cpx, gamma):
        return (self.mpow(-1.0, i) * self.binomialPrefactor(i, l1, l2, pax, pbx) \
                 * self.mpow(-1.0, u) * self.factorial(i) \
                 * self.mpow(cpx, i-2*r-2*u) \
                 * self.mpow(0.25/gamma, r+u) \
                 / self.factorial(r) \
                 / self.factorial(u) \
                 / self.factorial(i-2*r-2*u))
        
    def conAArr(self, l1, l2, pa, pb, cp, gamma):
        iMax = l1 + l2 + 1
        a = []
        for i in range(0, iMax): a.append(0.0)      
        
        for i in range(0, iMax):
            for r in range(0, int(floor(i/2.0)+1.0)):
                for u in range(0, int(floor((i-2.0*r)/2.0)+1.0)):
                    idx = i-2 * r-u                    
                    a[idx] += self.conATrm(i, r, u, l1, l2, pa, pb, cp, gamma)
        return a

    def gammln(self, x):
        y = x; tmp = x+5.5        
        tmp -= (x+0.5)*log(tmp)        
        ser = 1.000000000190015
        for j in range(0, 6):
            y+=1; ser += self.cof[j]/y
            
        return (-tmp+log((2.5066282746310005*ser)/x))
        
    def gammaIncomplete(self, a, x):
        gammap = 0.0; gln = self.gammln(a)
        
        if(x < (a + 1.0)):
            if(x != 0.0):
                ap = a; delta = sum = 1.0/a                
                for i in range(0, self.MAX_ITERATION):
                    ap+=1; delta *= x/ap; sum += delta
                    if(abs(delta) < abs(sum)*self.EPS): break                
                gammap = sum * exp(-x + a*log(x)-gln)
            else: gammap = 0.0
        else:
            b = (x+1.0)-a; c = 1.0/self.FPMIN
            d = 1.0/b; h=d            
            for i in range(1, self.MAX_ITERATION+1):
                an = -i*(i-a); b += 2.0; d = an * d + b
                if(abs(d) < self.FPMIN): d = self.FPMIN
                c = b+an/c
                if(abs(c) < self.FPMIN): c = self.FPMIN                
                d = 1.0/d; delta = d*c; h *= delta;                
                if(abs(delta-1.0) < self.EPS): break            
            gammap = 1.0-(exp(-x+a*log(x)-gln)*h)
        
        return (exp(gln)*gammap)
        
    def fGamma(self, m, x):
        x = max(abs(x), self.SMALL)        
        return (0.5 * self.mpow(x, -m-0.5) * self.gammaIncomplete(m+0.5, x))

    def nuc(self, a, nrm1, pwr1, alp1, b, nrm2, pwr2, alp2, c):
        prd = self.gaussianProductCenter(alp1, a, alp2, b)
        rab2 = self.dist2(a, b)
        rcp2 = self.dist2(c, prd)        
        gamma = alp1 + alp2

        ax = self.conAArr(pwr1[0], pwr2[0], prd[0]-a[0], prd[0]-b[0], prd[0]-c[0], gamma)
        ay = self.conAArr(pwr1[1], pwr2[1], prd[1]-a[1], prd[1]-b[1], prd[1]-c[1], gamma)
        az = self.conAArr(pwr1[2], pwr2[2], prd[2]-a[2], prd[2]-b[2], prd[2]-c[2], gamma)

        sum = 0.0
        for i in range(0, len(ax)):
            for j in range(0, len(ay)):
                for k in range(0, len(az)):
                    sum += ax[i]*ay[j]*az[k]*self.fGamma(i+j+k, rcp2*gamma)
                    
        return (-nrm1*nrm2*2.0*pi/gamma*exp(-alp1*alp2*rab2/gamma)*sum)

    def kinetic(self, alp1, pwr1, a, alp2, pwr2, b):
        l1 = pwr1[0]; m1 = pwr1[1]; n1 = pwr1[2]
        l2 = pwr2[0]; m2 = pwr2[1]; n2 = pwr2[2]        
        term = alp2*(2*(l2+m2+n2)+3)*self.ovrlp(alp1, pwr1, a, alp2, pwr2, b)        
        term += -2.0 * self.mpow(alp2, 2.0) \
                * (self.ovrlp(alp1, pwr1, a, alp2, (l2+2, m2, n2), b) \
                   + self.ovrlp(alp1, pwr1, a, alp2, (l2, m2+2, n2), b) \
                   + self.ovrlp(alp1, pwr1, a, alp2, (l2, m2, n2+2), b))        
        term += -0.5 * ((l2 * (l2-1)) * self.ovrlp(alp1, pwr1, a, alp2, (l2-2, m2, n2), b) \
                        + (m2 * (m2-1)) * self.ovrlp(alp1, pwr1, a, alp2, (l2, m2-2, n2), b) \
                        + (n2 * (n2-1)) * self.ovrlp(alp1, pwr1, a, alp2, (l2, m2, n2-2), b))        
        return term

    def funB(self, i, l1, l2, p, a, b, r, g):
        return (self.binomialPrefactor(i, l1, l2, p-a, p-b) * self.funB0(i, r, g))

    def funB0(self, i, r, g):
        return (self.factorialRatioSquared(i, r) * self.mpow(4*g, r-i))
    
    def conBTrm(self, i1, i2, r1, r2, u, l1, l2, l3, l4, px, ax, bx, \
                qx, cx, dx, gamma1, gamma2, delta):         
        return (self.funB(i1, l1, l2, px, ax, bx, r1, gamma1) \
                * pow(-1, i2) \
                * self.funB(i2, l3, l4, qx, cx, dx, r2, gamma2) \
                * pow(-1, u) \
                * self.factorialRatioSquared(i1+i2-2*(r1+r2), u) \
                * self.mpow(qx-px, i1+i2-2*(r1+r2)-2*u) \
                / self.mpow(delta, i1+i2-2*(r1+r2)-u))
        
    def conBArr(self, l1, l2, l3, l4, p, a, b, q, c, d, g1, g2, delta):
        iMax = l1+l2+l3+l4+1;
        bArr = []
        for i in range(0, iMax): bArr.append(0.0)
        
        for i1 in range(0, (l1+l2+1)):
            for i2 in range(0, (l3+l4+1)):
                for r1 in range(0, (i1/2+1)):
                    for r2 in range(0, (i2/2+1)):
                        for u in range(0, ((i1+i2)/2-r1-r2+1)):
                            idx = i1+i2-2*(r1+r2)-u
                            bArr[idx] += self.conBTrm(i1, i2, r1, r2, u, \
                                          l1, l2, l3, l4, p, a, b, q, c, d, g1, g2, delta)        
        return bArr
        
    def coulombRepulsion(self, a, aNorm, aPower, aAlpha, \
           b, bNorm, bPower, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dPower, dAlpha):               
        rab2 = self.dist2(a, b)
        rcd2 = self.dist2(c, d)        
        p = self.gaussianProductCenter(aAlpha, a, bAlpha, b);
        q = self.gaussianProductCenter(cAlpha, c, dAlpha, d);        
        rpq2 = self.dist2(p, q)        
        gamma1 = aAlpha + bAlpha;
        gamma2 = cAlpha + dAlpha;
        delta  = 0.25 * (1/gamma1 + 1/gamma2);
        
        bx = self.conBArr(aPower[0], bPower[0], cPower[0], dPower[0], \
                   p[0], a[0], b[0], q[0], c[0], d[0], gamma1, gamma2, delta)
        by = self.conBArr(aPower[1], bPower[1], cPower[1], dPower[1], \
                   p[1], a[1], b[1], q[1], c[1], d[1], gamma1, gamma2, delta)
        bz = self.conBArr(aPower[2], bPower[2], cPower[2], dPower[2], \
                   p[2], a[2], b[2], q[2], c[2], d[2], gamma1, gamma2, delta)
        
        sum = 0.0
        for i in range(0, len(bx)):
            for j in range(0, len(by)):
                for k in range(0, len(bz)):
                    sum += bx[i] * by[j] * bz[k] * self.fGamma(i+j+k, 0.25*rpq2/delta)
        
        return (2 * pow(pi, 2.5) \
                  / (gamma1 * gamma2 * sqrt(gamma1+gamma2)) \
                  * exp(-aAlpha*bAlpha*rab2/gamma1) \
                  * exp(-cAlpha*dAlpha*rcd2/gamma2) \
                  * sum * aNorm * bNorm * cNorm * dNorm)
        
    def coulomb(self, a, b, c, d):
         jij = 0.0         
         aExps = a.exps; aCoefs = a.cofs; aNorms = a.primnrm
         aOrigin = a.org; aPower = a.pwr

         bExps = b.exps; bCoefs = b.cofs; bNorms = b.primnrm
         bOrigin = b.org; bPower = b.pwr

         cExps = c.exps; cCoefs = c.cofs; cNorms = c.primnrm
         cOrigin = c.org; cPower = c.pwr         

         dExps = d.exps; dCoefs = d.cofs; dNorms = d.primnrm
         dOrigin = d.org; dPower = d.pwr
                   
         for i in range(0, len(aExps)):
             iaCoef = aCoefs[i]; iaExp = aExps[i]; iaNorm = aNorms[i]
             for j in range(0, len(bExps)):
                 jbCoef = bCoefs[j]; jbExp = bExps[j]; jbNorm = bNorms[j]
                 for k in range(0, len(cExps)):
                     kcCoef = cCoefs[k]; kcExp = cExps[k]; kcNorm = cNorms[k]             
                     for l in range(0, len(dExps)):
                        repulsionTerm = self.coulombRepulsion( \
                                         aOrigin, iaNorm, aPower, iaExp, \
                                         bOrigin, jbNorm, bPower, jbExp, \
                                         cOrigin, kcNorm, cPower, kcExp, \
                                         dOrigin, dNorms[l], dPower, dExps[l])                        
                        jij += iaCoef*jbCoef*kcCoef*dCoefs[l]*repulsionTerm
         
         return (a.norfac*b.norfac*c.norfac*d.norfac*jij)
         
    def ijkl2intindex(self, i, j, k, l):
        temp = 0
        
        if (i<j):
            temp=i; i=j; j=temp
        if (k<l):
            temp=k; k=l; l=temp
        
        ij = i*(i+1)/2+j
        kl = k*(k+1)/2+l
        
        if (ij < kl):
            temp=ij; ij=kl; kl=temp
        
        return (ij * (ij+1) / 2+kl)

integral = None

try:
    import cmobihf
    from cmobihf import *
    
    class CIntegral:
        def __init__(self): pass
        
        def factorial(self, n):
            return fact(n)
    
        def factorial2(self, n):
            return fact2(n)
    
        def factorialRatioSquared(self, a, b):
            return (self.factorial(a) / self.factorial(b) / self.factorial(a-2*b))
    
        def binomial(self, i, j):
            return (self.factorial(i) / self.factorial(j) / self.factorial(i-j))

        def mpow(self, x, y):
            if (round(y,10)==0.0): return 1.0
            if (round(x,10)==0.0): return 0.0
            return pow(x, y)
        
        def binomialPrefactor(self, s, ia, ib, xpa, xpb):
            return binomial_prefactor(s, ia, ib, xpa, xpb)
    
        def dist2(self, a, b):
            return dist2((a[0],a[1],a[2]),(b[0],b[1],b[2]))
        
        def dist(self, a, b):
            return dist((a[0],a[1],a[2]),(b[0],b[1],b[2]))
        
        def ovrlp(self, alp1, pwr, a, alp2, pwr2, b):
            return overlap(alp1, (pwr[0],pwr[1],pwr[2]), (a[0],a[1],a[2]), \
                           alp2, (pwr2[0],pwr2[1],pwr2[2]), (b[0],b[1],b[2]))
        
        def fGamma(self, m, x):
            return FGamma(m, x)
        
        def nuc(self, a, nrm1, pwr1, alp1, b, nrm2, pwr2, alp2, c):
            return nuclear_attraction((a[0],a[1],a[2]), nrm1, (pwr1[0],pwr1[1],pwr1[2]), alp1, \
                                      (b[0],b[1],b[2]), nrm2, (pwr2[0],pwr2[1],pwr2[2]), alp2, \
                                      (c[0],c[1],c[2]))
    
        def kinetic(self, alp1, pwr1, a, alp2, pwr2, b):
            return kinetic(alp1, (pwr1[0],pwr1[1],pwr1[2]), (a[0],a[1],a[2]), \
                           alp2, (pwr2[0],pwr2[1],pwr2[2]), (b[0],b[1],b[2]))
        
        def coulombRepulsion(self, a, aNorm, aPower, aAlpha, \
                             b, bNorm, bPower, bAlpha, c, cNorm, cPower, cAlpha, d, dNorm, dPower, dAlpha):               
            return coulomb_repulsion((a[0],a[1],a[2]), aNorm, (aPower[0],aPower[1],aPower[2]), aAlpha, \
                                     (b[0],b[1],b[2]), bNorm, (bPower[0],bPower[1],bPower[2]), bAlpha, \
                                     (c[0],c[1],c[2]), cNorm, (cPower[0],cPower[1],cPower[2]), cAlpha, \
                                     (d[0],d[1],d[2]), dNorm, (dPower[0],dPower[1],dPower[2]), dAlpha)
        
        def coulomb(self, a, b, c, d):
            return contr_coulomb(a.exps, a.cofs, a.primnrm, (a.org[0],a.org[1],a.org[2]), (a.pwr[0],a.pwr[1],a.pwr[2]), \
                                 b.exps, b.cofs, b.primnrm, (b.org[0],b.org[1],b.org[2]), (b.pwr[0],b.pwr[1],b.pwr[2]), \
                                 c.exps, c.cofs, c.primnrm, (c.org[0],c.org[1],c.org[2]), (c.pwr[0],c.pwr[1],c.pwr[2]), \
                                 d.exps, d.cofs, d.primnrm, (d.org[0],d.org[1],d.org[2]), (d.pwr[0],d.pwr[1],d.pwr[2]))
        
        def ijkl2intindex(self, i, j, k, l):
            return ijkl2intindex(i, j, k, l)

    integral = CIntegral()
except:
    integral = Integral()
        
class pg:
    def __init__(self, org, pwr, exp, cof):
        self.org=org; self.pwr=pwr
        self.exp=exp; self.cof=cof
        self.norm()

    def ovrlp(self, og):
        return (self.norfac * og.norfac * \
                 integral.ovrlp(self.exp, self.pwr, self.org, og.exp, og.pwr, og.org))

    def norm(self):
        l=self.pwr[0]; m=self.pwr[1]; n=self.pwr[2]
        self.norfac = sqrt(pow(2, 2*(l+m+n)+1.5) \
                           * pow(self.exp, l+m+n+1.5) \
                           / integral.factorial2(2*l-1) \
                           / integral.factorial2(2*m-1) \
                           / integral.factorial2(2*n-1) / pow(pi, 1.5))

    def nuc(self, og, cen):
        return (integral.nuc(self.org, self.norfac, self.pwr, \
                 self.exp, og.org, og.norfac, og.pwr, og.exp, cen))

    def kinetic(self, og):
        return (self.norfac*og.norfac*integral.kinetic( \
                self.exp, self.pwr, self.org, og.exp, og.pwr, og.org))
        
class cg:
    def __init__(self, org, pwr):
        self.org=org; self.pwr=pwr
        self.pgs = []; self.exps = []; self.cofs = []
        self.norfac = 1.0
        self.primnrm = []

    def addpg(self, exp, cof):        
        self.pgs.append(pg(self.org, self.pwr, exp, cof))
        self.exps.append(exp); self.cofs.append(cof)

    def norm(self):
        self.norfac = 1.0 / sqrt(self.ovrlp(self))        
        for i in range(0, len(self.pgs)):
            self.primnrm.append(self.pgs[i].norfac)

    def ovrlp(self, og):
        sij = 0.0
        for i in range(0, len(self.pgs)):
            ipg = self.pgs[i]
            for j in range(0, len(og.pgs)):
                jpg = og.pgs[j]
                sij += ipg.cof*jpg.cof*ipg.ovrlp(jpg)
        return self.norfac*og.norfac*sij

    def kinetic(self, og):
        tij=0.0
        for i in range(0, len(self.pgs)):
            ipg = self.pgs[i]
            for j in range(0, len(og.pgs)):
                jpg = og.pgs[j]
                tij += ipg.cof*jpg.cof*ipg.kinetic(jpg)
        return self.norfac*og.norfac*tij

    def nuc(self, og, cen):
        vij=0.0
        for i in range(0, len(self.pgs)):
            ipg = self.pgs[i]
            for j in range(0, len(og.pgs)):
                jpg = og.pgs[j]
                vij += ipg.cof*jpg.cof*ipg.nuc(jpg, cen)
        return self.norfac*og.norfac*vij
    
class atmbas:
    def __init__(self, sym, num):
        self.sym=sym; self.num=num;
        self.orbs=[]
                 
class basis:
    def __init__(self, mol, bas):
        self.bfs=[]

        for atm in mol.atoms:
            atmbas = curbas[atm[0]]
            for ab in atmbas:
                pwrlst = powers[ab[0]]
                for pwr in pwrlst:
                    c = cg(atm[1], pwr)
                    for orb in ab[1]:
                        c.addpg(orb[0], orb[1])
                    c.norm()
                    self.bfs.append(c)
                
class oneEI:
    def __init__(self, bfs, mol):
        self.H = []; self.S = []        
        nbf = len(bfs)
        for i in range(0, nbf):
            self.H.append([]); self.S.append([])
            for j in range(0, nbf):
                self.H[i].append(0.0);
                self.S[i].append(0.0);
                
        for i in range(0, nbf):
            ibf=bfs[i]
            for j in range(0, nbf):
                jbf=bfs[j]
                self.S[i][j] = ibf.ovrlp(jbf)
                self.H[i][j] = ibf.kinetic(jbf)
                
                for k in range(0, len(mol.atoms)):
                    self.H[i][j] += atmnum[mol.atoms[k][0]] * ibf.nuc(jbf, mol.atoms[k][1])

class twoEI:
    def __init__(self, bfs):
        nbf = len(bfs)
        self.nint = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
        self.ints = []
        for i in range(0, self.nint): self.ints.append(0.0)

        for i in range(0, nbf):
            bfi = bfs[i]            
            for j in range(0, (i+1)):
                bfj = bfs[j]
                ij = i*(i+1)/2+j                
                for k in range(0, nbf):
                    bfk = bfs[k]                    
                    for l in range(0, (k+1)):
                        bfl = bfs[l]                        
                        kl = k*(k+1)/2+l;
                        if (ij >= kl):
                            ijkl = integral.ijkl2intindex(i, j, k, l)
                            self.ints[ijkl] = integral.coulomb(bfi, bfj, bfk, bfl)

class hf:
    def __init__(self, mol, e1, e2):
        self.e1=e1; self.e2=e2
        ele=0
        for i in range(0, len(mol.atoms)):
            ele += atmnum[mol.atoms[i][0]]
        self.nocc=ele/2

        self.eNuc=0.0
        for i in range(0, len(mol.atoms)):
            ati = mol.atoms[i]
            for j in range(0, i):
                atj = mol.atoms[j]                
                self.eNuc += atmnum[ati[0]] * atmnum[atj[0]] / integral.dist(ati[1], atj[1])
                
    def doRotate(self, a, i, j, k, l, sin, tau):
        g = a[i][j]; h = a[k][l]        
        a[i][j] = g - sin * (h + g*tau)
        a[k][l] = h + sin * (g - h*tau)

    def sortEval(self, eval, evec):
        n = len(evec)
        for i in range(0, n):
            k=i; p=eval[i]
            for j in range(i+1, n):
                if (eval[j] <= p):
                    k=j; p=eval[j]
            if (k!=i):
                eval[k] = eval[i]; eval[i] = p                
                for j in range(0, n):
                    p = evec[j][i]; evec[j][i] = evec[j][k]; evec[j][k] = p
        return eval, evec
                    
    def diagonalize(self, m):
        maxItr=50
        evec=[]; a=[]; n=len(m)
        for i in range(0, n):
            evec.append([]); a.append([])
            for j in range(0, n):
                evec[i].append(0.0)
                a[i].append(m[i][j])
            evec[i][i] = 1.0
        eval=[]; b=[]; z=[]        
        for i in range(0, n):
            eval.append(a[i][i]); b.append(a[i][i]); z.append(0.0)
        for sweeps in range(0, maxItr):
            sum = 0.0
            for i in range(0, n-1):
                for j in range(i+1, n):
                    sum += abs(a[i][j])
            if (sum == 0.0): break
            
            if (sweeps < 3): zeroTol = 0.2*sum/(n*n)
            else: zeroTol = 0.0
            
            for ip in range(0, n-1):
                for iq in range(ip+1, n):
                    g = 100.0 * abs(a[ip][iq])
                    
                    if ((sweeps > 4) and (float(abs(eval[ip])+g)==float(abs(eval[ip]))) \
                        and (float(abs(eval[iq])+g)==float(abs(eval[iq])))): a[ip][iq] = 0.0
                    elif (abs(a[ip][iq]) > zeroTol):
                        h = eval[iq]-eval[ip]                        
                        if (float((abs(h)+g))==float(abs(h))): t = a[ip][iq] / h
                        else:
                            theta = 0.5 * h / a[ip][iq]
                            t = 1.0 / (abs(theta)+sqrt(1.0 + theta*theta))                            
                            if (theta < 0.0): t = -t
                            
                        cos = 1.0 / sqrt(1.0 + t*t)
                        sin = t * cos
                        tau = sin / (1.0 + cos)
                        h   = t * a[ip][iq]                        
                        z[ip] -= h; z[iq] += h
                        eval[ip] -= h; eval[iq] += h                        
                        a[ip][iq] = 0.0
                                                
                        for j in range(0, ip): self.doRotate(a, j, ip, j, iq, sin, tau)
                        for j in range(ip+1, iq): self.doRotate(a, ip, j, j, iq, sin, tau)
                        for j in range(iq+1, n): self.doRotate(a, ip, j, iq, j, sin, tau)
                        for j in range(0, n): self.doRotate(evec, j, ip, j, iq, sin, tau)

            for ip in range(0, n):
                b[ip] += z[ip]; eval[ip] = b[ip]; z[ip] = 0.0
        eval, evec = self.sortEval(eval, evec)
        return eval, self.trans(evec)

    def trans(self, m):
        nm=[]; rc=len(m); cc=len(m[0])
        for i in range(0, cc):
            nm.append([])
            for j in range(0, rc):
                nm[i].append(m[j][i])
        return nm
    
    def mul(self, m1, m2):
        nm=[]; rc=len(m1); brc=len(m2); bcc=len(m2[0])
        for i in range(0, rc):
            nm.append([])
            for j in range(0, bcc):
                cij=0.0
                for k in range(0, brc):
                    cij += m1[i][k]*m2[k][j]
                nm[i].append(cij)
        return nm

    def add(self, m1, m2):
        nm=[]; rc=len(m1);
        for i in range(0, rc):
            nm.append([])
            for j in range(0, rc):
                nm[i].append(m1[i][j]+m2[i][j])
        return nm

    def tr(self, m):
        r=0.0
        for i in range(len(m)): r+=m[i][i]
        return r
    
    def simiTransT(self, m1, m2):
        return self.mul(self.mul(self.trans(m1), m2), m1)
    
    def simiTrans(self, m1, m2):
        return self.mul(self.mul(m1, m2), self.trans(m1))
    
    def symOrth(self, m):
        eval, evec = self.diagonalize(m)
        sHalf=[]; rc=len(m)
        for i in range(0, rc):
            sHalf.append([])
            for j in range(0, rc):
                sHalf[i].append(0.0)
            sHalf[i][i] = 1.0 / sqrt(eval[i])
            
        return (self.simiTransT(evec, sHalf))
        
    def compOrbES(self, h, s):
        x = self.symOrth(s)
        a = self.simiTrans(x, h)
        eval, evec = self.diagonalize(a)
        return eval, self.mul(evec, x)

    def makD(self, mos):
        d = []
        for i in range(0, self.nocc):
            d.append([])
            for j in range(0, len(mos)):
                d[i].append(mos[i][j])
        return self.mul(self.trans(d), d)

    def dot(self, v1, v2):
        res = 0.0
        for i in range(0, len(v1[0])): res += v1[0][i]*v2[0][i]
        return res
    
    def makG(self, mos, D):
        nbf=len(D); G=[]; oneD=[]; oneD.append([]); tmpV=[]; tmpV.append([])
        for i in range(0, nbf):
            G.append([])
            for j in range(0, nbf):
                G[i].append(0.0)
                oneD[0].append(D[i][j])
                tmpV[0].append(0.0)
                
        for i in range(0, nbf):
            for j in range(0, i+1):
                for k in range(0, nbf*nbf): tmpV[0][i]=0.0
                kl = 0                
                for k in range(0, nbf):
                    for l in range(0, nbf):
                        idxJ   = integral.ijkl2intindex(i, j, k, l)
                        idxK1  = integral.ijkl2intindex(i, k, j, l)
                        idxK2  = integral.ijkl2intindex(i, l, k, j)
                        tmpV[0][kl] = 2.0*self.e2.ints[idxJ]-0.5*self.e2.ints[idxK1]-0.5*self.e2.ints[idxK2]
                        kl+=1
                G[i][j] = G[j][i] = self.dot(tmpV, oneD)
                
        return G
    
    def scf(self):                
        maxItr=30; eneTol=1.0e-4
        orbE, mos = self.compOrbES(self.e1.H, self.e1.S)
        oldEnergy=0.0
        for scfIteration in range(0, maxItr):
            D = self.makD(mos);  
            G = self.makG(mos, D)
            F = self.add(self.e1.H, G)
            orbE, mos = self.compOrbES(F, self.e1.S)
            eOne = self.tr(self.mul(D, self.e1.H))
            eTwo = self.tr(self.mul(D, F))
            energy = eOne + eTwo + self.eNuc
            print scfIteration, energy
            if (abs(energy - oldEnergy) < eneTol): break            
            oldEnergy = energy
                            
class mobihf:
    def __init__(self):
        self.filename = "c:\\system\\apps\\python\\my\\h2.inp"
        self.filename = "e:\\others\\input\\wat.inp"
        #self.filename = "h2.inp"

    def run(self):
        self.do_run()

    def do_run(self):
        self.read_inp()
        self.chk_inp()
        self.bfs = basis(self.mol, self.basis).bfs
        self.ei1 = oneEI(self.bfs, self.mol)
        self.ei2 = twoEI(self.bfs)
        hf(self.mol, self.ei1, self.ei2).scf()

    def read_inp(self):
        global radii
        xyzf = open(self.filename, "r")
        lines = xyzf.readlines()
        xyzf.close()
        
        self.mol = molecule()
        self.mol.noOfAtoms = atoi(lines[0])
        self.mol.title, self.level, self.basis = split(strip(lines[1]))
        self.level = upper(self.level); self.basis = upper(self.basis)
        for i in range(2, len(lines)):
            words = split(lines[i])
            self.mol.atoms.append((upper(words[0]), map(atof, words[1:]), []))
         
        for i in range(0, self.mol.noOfAtoms):
            atm1 = self.mol.atoms[i]
            for j in range(0, i):
                atm2 = self.mol.atoms[j] 
                x = atm1[1][0]-atm2[1][0]
                y = atm1[1][1]-atm2[1][1]
                z = atm1[1][2]-atm2[1][2]
                dist = sqrt(x*x+y*y+z*z)
                radsum = radii[atm1[0]] + radii[atm2[0]]
                if (((radsum - 0.4) < dist) and (dist < (radsum + 0.4))):
                    atm1[2].append(j)
                    atm2[2].append(i)
                     
        lines = None

        print self.mol.title, self.level, self.basis        

    def chk_inp(self):
        if (self.level != "HF"):
            print "Unsupported level"
        if (self.basis != "STO-3G"):
            print "Unsupported basis"
            
if __name__ == '__main__':
    mobihf().run()
