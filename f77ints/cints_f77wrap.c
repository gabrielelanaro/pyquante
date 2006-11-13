#include <stdio.h>
#include <stdlib.h>
#include "cints.h"

void coulomb_repulsion_f__(double *answ,
			   double *xa, double *ya, double *za, double *norma,
			   int *la, int *ma, int *na, double *alphaa,
			   double *xb, double *yb, double *zb, double *normb,
			   int *lb, int *mb, int *nb, double *alphab,
			   double *xc, double *yc, double *zc, double *normc,
			   int *lc, int *mc, int *nc, double *alphac,
			   double *xd, double *yd, double *zd, double *normd,
			   int *ld, int *md, int *nd, double *alphad){

  *answ = coulomb_repulsion(*xa,  *ya,  *za,  *norma,
			    *la,  *ma,  *na,  *alphaa,
			    *xb,  *yb,  *zb,  *normb,
			    *lb,  *mb,  *nb,  *alphab,
			    *xc,  *yc,  *zc,  *normc,
			    *lc,  *mc,  *nc,  *alphac,
			    *xd,  *yd,  *zd,  *normd,
			    *ld,  *md,  *nd,  *alphad);
}

