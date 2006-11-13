      program test_cint

      implicit double precision(a-h,o-z)

c     This program tests calling the cints library from fortran.
cp1 = PGBF(0.5,(0,0,0))
      xa = 0.d0
      ya = 0.d0
      za = 0.d0
      la = 0
      ma = 0
      na = 0
      alphaa = 0.5d0
      xnorma = 0.423777208124d0

cp2 = PGBF(1.0,(0,0,1))
      xb = 0.
      yb = 0.
      zb = 1.
      lb = 0
      mb = 0
      nb = 0
      alphab = 1.0
      xnormb = 0.712705470355 


cp3 = PGBF(0.1,(0,0,1))
      xc = 0.
      yc = 0.
      zc = 1.
      lc = 0
      mc = 0
      nc = 0
      alphac = 0.1
      xnormc = 0.126738946335 

cp4 = PGBF(0.2,(0,1,0))
      xd = 0.
      yd = 1.
      zd = 0.
      ld = 0
      md = 0
      nd = 0
      alphad = 0.2
      xnormd = 0.213148651293


      print*, "Testing cints coulomb repulsion from Fortran"
      call coulomb_repulsion_f(answ,xa, ya, za, norma,la, ma, na, alphaa,
     $     xb, yb, zb, normb,lb, mb, nb, alphab,
     $     xc, yc, zc, normc,lc, mc, nc, alphac,
     $     xd, yd, zd, normd,ld, md, nd, alphad)
      print*, answ
		
      stop
      end
