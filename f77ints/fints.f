c     Module fints.f
c
c     A F77 implementation of the functions in the cints/pyints
c     routines of PyQuante.
c     This program is part of the PyQuante quantum chemistry suite.
c     Copyright 2002 California Institute of Technology.  All Rights
c     Reserved.  You may contact the author, Richard P. Muller at
c     rpm@wag.caltech.edu.
c
c     This program is free software; you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the Free Software Foundation provided that the above copyright notice
c     and this complete notice appear in all copies.
c
c     This program is distributed in the hope that it will be useful, and in
c     no event shall California Institute of Technology be liable to any
c     party for direct, indirect, special, incidental or consequential
c     damages, including lost profits, arising out of the use of this
c     software and its documentation, even if the California Institute of
c     Technology has been advised of the possibility of such damage. The
c     California Institute of Technology specifically disclaims any
c     warranties, including the implied warranties or merchantability and
c     fitness for a particular purpose. The software and documentation
c     provided hereunder is on an AS IS basis, and the California Institute
c     of Technology has no obligations to provide maintenance, support,
c     updates, enhancements or modifications.
c
c     You should have received a copy of the GNU General Public License
c     along with this program; if not, write to the Free Software
c     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307

c     Functions in this module:
c     coulomb_repulsion: Compute the coulomb repulsion between four
c       primitive gaussian functions.
c     overlap: Compute the overlap between two primitive gaussians
c     kinetic: Compute the kinetic energy matrix element between two
c      primitive gaussian functions

c--------------------------------------------------------------------------

      real*8 function overlap(alpha1,l1,m1,n1,x1,y1,z1,
     $     alpha2,l2,m2,n2,x2,y2,z2)
c     Note: this function does multiply by the normalization constant

      implicit none

c     Arguments and local variables
      real*8 alpha1,x1,y1,z1,alpha2,x2,y2,z2,
     $     rab2,gamma,xp,yp,zp,factor,wx,wy,wz,pi
      integer l1,m1,n1,l2,m2,n2

c     Functions
      real*8 product_center_1D,dist2,overlap1D

c     Constants
      pi = 3.1415927

      rab2 = dist2(x1,y1,z1,x2,y2,z2)

      gamma = alpha1+alpha2
      xp = product_center_1D(alpha1,x1,alpha2,x2)
      yp = product_center_1D(alpha1,y1,alpha2,y2)
      zp = product_center_1D(alpha1,z1,alpha2,z2)

      factor = ((pi/gamma)**1.5)*exp(-alpha1*alpha2*rab2/gamma)

      wx = overlap1D(l1,l2,xp-x1,xp-x2,gamma)
      wy = overlap1D(m1,m2,yp-y1,yp-y2,gamma)
      wz = overlap1D(n1,n2,zp-z1,zp-z2,gamma)
      overlap = factor*wx*wy*wz;
      return
      end

c--------------------------------------------------------------------------

      real*8 function kinetic(alpha1,l1,m1,n1,x1,y1,z1,
     $     alpha2,l2,m2,n2,x2,y2,z2)

      implicit none

c     Arguments and local variables
      real*8 alpha1,x1,y1,z1,alpha2,x2,y2,z2,
     $     term0,term1,term2x,term2y,term2z
      integer l1,m1,n1,l2,m2,n2

c     Constants
      real*8 zero,one,two,three,half

c     Functions
      real*8 overlap

      zero = 0.0d0
      one = 1.0d0
      two = 2*one
      three = 3*one
      half = one/two

      term0 = alpha2*(two*(l2+m2+n2)+three)*
     $     overlap(alpha1,l1,m1,n1,x1,y1,z1,
     $     alpha2,l2,m2,n2,x2,y2,z2)

      term1 = -two*alpha2*alpha2*(
     $     overlap(alpha1,l1,m1,n1,x1,y1,z1,alpha2,l2+2,m2,n2,x2,y2,z2)+
     $     overlap(alpha1,l1,m1,n1,x1,y1,z1,alpha2,l2,m2+2,n2,x2,y2,z2)+
     $     overlap(alpha1,l1,m1,n1,x1,y1,z1,alpha2,l2,m2,n2+2,x2,y2,z2))

      if (l2 .ge. 2) then
         term2x = -half*l2*(l2-1)*overlap(alpha1,l1,m1,n1,x1,y1,z1,
     $        alpha2,l2-2,m2,n2,x2,y2,z2)
      else
         term2x = zero
      endif

      if (m2 .ge. 2) then
         term2y = -half*m2*(m2-1)*overlap(alpha1,l1,m1,n1,x1,y1,z1,
     $        alpha2,l2,m2-2,n2,x2,y2,z2)
      else
         term2y = zero
      endif

      if (n2 .ge. 2) then
         term2z = -half*n2*(n2-1)*overlap(alpha1,l1,m1,n1,x1,y1,z1,
     $        alpha2,l2,m2,n2-2,x2,y2,z2)
      else
         term2z = zero
      endif

      kinetic = term0+term1+term2x+term2y+term2z

      return
      end

c--------------------------------------------------------------------------

      real*8 function overlap1D(l1,l2,xap,xbp,gamma)
      implicit none
      integer l1,l2,i
      real*8 xap,xbp,gamma,sum
      real*8 binomial_prefactor
      integer fact2
      sum = 0.0d0
      do i=0, int(0.5*(l1+l2))
         sum = sum + binomial_prefactor(2*i,l1,l2,xap,xbp)*
     $        fact2(2*i-1)/((2*gamma)**i)
      enddo
      overlap1D = sum
      return 
      end

c--------------------------------------------------------------------------

      real*8 function coulomb_repulsion(
     $     xa, ya, za, norma, la, ma, na, alphaa,
     $     xb, yb, zb, normb, lb, mb, nb, alphab,
     $     xc, yc, zc, normc, lc, mc, nc, alphac,
     $     xd, yd, zd, normd, ld, md, nd, alphad)

      implicit none

c     Parameters
      integer MAXLEN
      parameter (MAXLEN=25)

c     Arguments
      real*8 xa,ya,za,norma,alphaa,xb,yb,zb,normb,alphab,
     $     xc,yc,zc,normc,alphac,xd,yd,zd,normd,alphad
      integer la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd

c     Local variables
      real*8 rab2, rcd2,rpq2,xp,yp,zp,xq,yq,zq,gamma1,gamma2,delta,sum,
     $     xval
      real*8 Bx(0:MAXLEN-1), By(0:MAXLEN-1), Bz(0:MAXLEN-1)
      integer i,j,k

c     Functions
      real*8 dist2,product_center_1D,Fgamma

c     Constants
      real*8 pi

      pi = 3.1415927

      rab2 = dist2(xa,ya,za,xb,yb,zb)
      rcd2 = dist2(xc,yc,zc,xd,yd,zd)
      xp = product_center_1D(alphaa,xa,alphab,xb)
      yp = product_center_1D(alphaa,ya,alphab,yb)
      zp = product_center_1D(alphaa,za,alphab,zb)
      xq = product_center_1D(alphac,xc,alphad,xd)
      yq = product_center_1D(alphac,yc,alphad,yd)
      zq = product_center_1D(alphac,zc,alphad,zd)
      rpq2 = dist2(xp,yp,zp,xq,yq,zq)
      gamma1 = alphaa+alphab
      gamma2 = alphac+alphad
      delta = (1./gamma1+1./gamma2)/4.

      call B_array(la+lb+lc+ld+1,Bx,
     $     la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta)
      call B_array(ma+mb+mc+md+1,By,
     $     ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta)
      call B_array(na+nb+nc+nd+1,Bz,
     $     na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta)

      sum = 0.d0
      do i = 0,la+lb+lc+ld
         do j = 0, ma+mb+mc+md
            do k = 0, na+nb+nc+nd
               xval = float(i+j+k) ! cast to double
               sum = sum +
     $              Bx(i)*By(j)*Bz(k)
     $              *Fgamma(xval,0.25*rpq2/delta)
            enddo
         enddo
      enddo

      coulomb_repulsion = 2.*pi**2.5
     $     /(gamma1*gamma2*sqrt(gamma1+gamma2))
     $     *exp(-alphaa*alphab*rab2/gamma1)
     $     *exp(-alphac*alphad*rcd2/gamma2)
     $     *sum*norma*normb*normc*normd

      return
      end

c--------------------------------------------------------------------------

      real*8 function dist2(xa,ya,za,xb,yb,zb)

      implicit none
      real*8 xa,ya,za,xb,yb,zb,dx,dy,dz

      dx = xa-xb
      dy = ya-yb
      dz = za-zb
      dist2 = dx*dx+dy*dy+dz*dz

      return
      end

c--------------------------------------------------------------------------

      real*8 function dist(xa,ya,za,xb,yb,zb)

      implicit none
      real*8 xa,ya,za,xb,yb,zb

      real*8 dist2

      dist = sqrt(dist2(xa,ya,za,xb,yb,zb))

      return
      end

c--------------------------------------------------------------------------

      real*8 function product_center_1D(alphaa,xa,alphab,xb)

      implicit none

      real*8 alphaa,xa,alphab,xb

      product_center_1D = (alphaa*xa+alphab*xb)/(alphaa+alphab)

      return 
      end

c--------------------------------------------------------------------------

      subroutine B_array(n,barry,l1,l2,l3,l4,
     $     p,a,b,q,c,d,g1,g2,delta)

      implicit none

      integer n,l1,l2,l3,l4,i,i1,i2,j1,j2,k,index
      real*8 barry(0:n-1),p,a,b,q,c,d,g1,g2,delta,zero

      real*8 b_term ! External function

      zero = 0.

      do i = 0, n-1
         barry(i) = zero
      enddo

      do i1=0, l1+l2
         do i2=0, l3+l4
            do j1=0, i1/2
               do j2=0, i2/2
                  do k=0, (i1+i2)/2-j1-j2
                     index = i1+i2-2*(j1+j2)-k
                     barry(index) = barry(index)
     $                    + b_term(i1,i2,j1,j2,k,l1,l2,l3,l4,
     $                             p,a,b,q,c,d,g1,g2,delta)
                  enddo
               enddo
            enddo
         enddo
      enddo      

      return 
      end

c--------------------------------------------------------------------------

      real*8 function b_term(i1,i2,j1,j2,k,l1,l2,l3,l4,
     $     p,a,b,q,c,d,g1,g2,delta)

      integer i1,i2,j1,j2,k,l1,l2,l3,l4
      real*8 p,a,b,q,c,d,g1,g2,delta
      real*8 fB
      integer fact_ratio2

      b_term = fB(i1,l1,l2,p,a,b,j1,g1)
     $     *(-1**i2)*fB(i2,l3,l4,q,c,d,j2,g2)
     $     *(-1**k)*fact_ratio2(i1+i2-2*(j1+j2),k)
     $     *(q-p)**(i1+i2-2*(j1+j2)-2*k)
     $     /(delta)**(i1+i2-2*(j1+j2)-k)
      return
      end

c--------------------------------------------------------------------------

      real*8 function fB(i, l1, l2, px, ax, bx, ir, g)
      implicit none

      integer i,l1,l2,ir
      real*8 px,ax,bx,g

      real*8 binomial_prefactor,Bfunc

      fB = binomial_prefactor(i,l1,l2,px-ax,px-bx)*Bfunc(i,ir,g)

      return
      end

c--------------------------------------------------------------------------

      real*8 function Bfunc(i,ir,g)
      implicit none
      integer i,ir
      real*8 g
      integer fact_ratio2
      
      Bfunc = fact_ratio2(i,ir)*((4*g)**(ir-i))

      return
      end

c--------------------------------------------------------------------------

      real*8 function binomial_prefactor(s,ia,ib,xpa,xpb)
      implicit none
      integer s,ia,ib,j
      real*8 xpa,xpb,sum
      integer binomial

      sum = 0.
      do j=0,s
         if ((s-ia .le. j) .and. (j .le. ib))
     $        sum = sum + binomial(ia,s-j)*binomial(ib,j)
     $        *(xpa**(ia-s+j))*(xpb**(ib-j))
      enddo
      binomial_prefactor = sum
      return
      end

c--------------------------------------------------------------------------

      integer function binomial(i,j)
      implicit none
      integer i,j
      integer fact

      binomial = fact(i)/fact(j)/fact(i-j)      

      return
      end

c--------------------------------------------------------------------------

      integer function fact_ratio2(i,j)
      implicit none
      integer i,j
      integer fact

      fact_ratio2 = fact(i)/fact(j)/fact(i-2*j)

      return
      end

c--------------------------------------------------------------------------

      integer function fact2(n)
      implicit none
      integer i,n

      fact2=1
      if (n .le. 1) return
      do i = 1,n,2
         fact2 = fact2*i
      enddo
      return
      end

c--------------------------------------------------------------------------

      integer function fact(n)
      implicit none
      integer i,n

      fact=1
      if (n .le. 1) return
      do i=1,n
         fact = fact*i
      enddo
      return
      end

c--------------------------------------------------------------------------

      real*8 function Fgamma(a,x)
      implicit none
      real*8 a,x,val
      real*8 gammp

      if (x .lt. 0.00001d0) x = 0.000001d0

      val = gammp(a+0.5,x)
      Fgamma = 0.5*(x**(-a-0.5))*val

      return
      end

c--------------------------------------------------------------------------

c     The following are from Numerical Recipes
c     gammp is hacked a bit to return exp(gln)*gamma
      FUNCTION gammp(a,x)
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln

      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=exp(gln)*gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=exp(gln)*(1.-gammcf)
      endif

      return
      END

c--------------------------------------------------------------------------

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

c--------------------------------------------------------------------------

      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      real*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END

c--------------------------------------------------------------------------

      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

c--------------------------------------------------------------------------
      
