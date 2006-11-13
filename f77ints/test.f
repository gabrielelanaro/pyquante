      program test

      implicit none
      real*8 alphaa,xa,ya,za,xnorma,
     $     alphab,xb,yb,zb,xnormb,
     $     alphac,xc,yc,zc,xnormc,
     $     alphad,xd,yd,zd,xnormd,
     $     alphae,xe,ye,ze,xnorme,
     $     answ,one
      integer la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd,le,me,ne,i

      real*8 coulomb_repulsion,overlap,kinetic
      integer fact,fact2,binomial

      one = 1.0d0

      print*, "Testing null loop"
      do i=0,0
         print*, "In loop ",i
      enddo

      print*, "Testing factorial"
      print*, "fact(-1) = ",fact(-1)
      print*, "fact(0) = ",fact(0)
      print*, "fact(1) = ",fact(1)
      print*, "fact(2) = ",fact(2)
      print*, "fact2(-1) = ",fact2(-1)
      print*, "fact2(0) = ",fact2(0)
      print*, "fact2(1) = ",fact2(1)
      print*, "fact2(3) = ",fact2(3)

      print*, "Testing binomial"
      print*, binomial(2,1),binomial(2,2),binomial(3,1),binomial(3,2),
     $     binomial(3,3),binomial(3,0)



c     This program tests calling the cints library from fortran.
cp1 = PGBF(0.5,(0,0,0),(0,0,0))
      alphaa = 0.5d0
      xa = 0.d0
      ya = 0.d0
      za = 0.d0
      la = 0
      ma = 0
      na = 0
      xnorma = one/sqrt(overlap(alphaa,la,ma,na,xa,ya,za,
     $     alphaa,la,ma,na,xa,ya,za))

cp2 = PGBF(1.0,(0,0,1),(1,0,0))
      alphab = 1.0
      xb = 0.
      yb = 0.
      zb = 1.
      lb = 1
      mb = 0
      nb = 0
      xnormb = one/sqrt(overlap(alphab,lb,mb,nb,xb,yb,zb,
     $     alphab,lb,mb,nb,xb,yb,zb))


cp3 = PGBF(0.1,(0,0,1),(1,1,0))
      alphac = 0.1
      xc = 0.
      yc = 0.
      zc = 1.
      lc = 1
      mc = 1
      nc = 0
      xnormc = one/sqrt(overlap(alphac,lc,mc,nc,xc,yc,zc,
     $     alphac,lc,mc,nc,xc,yc,zc))

cp4 = PGBF(0.2,(0,1,0),(0,0,1))
      alphad = 0.2
      xd = 0.
      yd = 1.
      zd = 0.
      ld = 0
      md = 0
      nd = 1
      xnormd = one/sqrt(overlap(alphad,ld,md,nd,xd,yd,zd,
     $     alphad,ld,md,nd,xd,yd,zd))

      alphae = 0.5d0
      xe = 0.d0
      ye = 0.d0
      ze = 1.d0
      le = 0
      me = 0
      ne = 0
      xnorme = one/sqrt(overlap(alphae,le,me,ne,xe,ye,ze,
     $     alphae,le,me,ne,xe,ye,ze))

      print*, "Norms:"
      print*, xnorma,xnormb,xnormc,xnormd

      print*, "Testing overlap function...should be 1:"
      answ = xnorma*xnorma*overlap(alphaa,la,ma,na,xa,ya,za,
     $     alphaa,la,ma,na,xa,ya,za)
      print*, answ

      print*, "Testing overlap function...should be 0.778800783"
      print*, xnorma*xnorme*overlap(alphaa,la,ma,na,xa,ya,za,
     $     alphae,le,me,ne,xe,ye,ze)

      print*, "Testing kinetic energy matrix element,"
      print*, "...should be 0.75"
      answ = xnorma*xnorma*kinetic(alphaa,la,ma,na,xa,ya,za,
     $     alphaa,la,ma,na,xa,ya,za)
      print*, answ
       	
      print*, "Testing kinetic energy matrix element,"
      print*, "...should be 2.5"
      answ = xnormb*xnormb*kinetic(alphab,lb,mb,nb,xb,yb,zb,
     $     alphab,lb,mb,nb,xb,yb,zb)
      print*, answ
       	
      print*, "Testing coulomb repulsion from Fortran"
      print*, "...correct answer should be 0.00285248"
      answ = coulomb_repulsion(xa, ya, za, xnorma,la, ma, na, alphaa,
     $     xb, yb, zb, xnormb,lb, mb, nb, alphab,
     $     xc, yc, zc, xnormc,lc, mc, nc, alphac,
     $     xd, yd, zd, xnormd,ld, md, nd, alphad)
      print*, answ

      stop
      end
