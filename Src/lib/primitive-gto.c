/**
 * @file primitive-gto.c
 * @author Gabriele Lanaro <gabriele.lanaro@gmail.com>
 *
 * This source code is part of the PyQuante Quantum Chemistry suite.
 *  
 * Written by Gabriele Lanaro, 2009-2010
 * Copyright (c) 2009-2010, Gabriele Lanaro
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 */

#include "primitive-gto.h"
#include <stdlib.h>
#include <math.h>
#include "utils/math.h"
//#include "utils/crys.h"

/**
 * Allocates the PrimitiveGTO
 */
PrimitiveGTO *primitive_gto_new(void)
{
  return malloc(sizeof(PrimitiveGTO));
}

/**
 * Initialize the PrimitiveGTO
 */
void primitive_gto_init(PrimitiveGTO *gto, double alpha,
			double x0, double y0, double z0,
			int l, int m, int n,
			double coef)
{
    /*
     * Initialization
     */

    gto->x0 = x0;
    gto->y0 = y0;
    gto->z0 = z0;

    gto->l = l;
    gto->m = m;
    gto->n = n;

    gto->alpha = alpha;
    gto->coef = coef;

    gto->norm = 1.; /* Un-normalized */
    primitive_gto_normalize(gto); /* Normalizing in the constructo */
}

PrimitiveGTO *primitive_gto_copy(PrimitiveGTO *pgto)
{
  PrimitiveGTO *pgto_new;
  
  pgto_new = primitive_gto_new();
  primitive_gto_init(pgto_new, pgto->alpha,
		     pgto->x0, pgto->y0, pgto->z0,
		     pgto->l, pgto->m, pgto->n,
		     pgto->coef);
  return pgto_new;

}

void primitive_gto_free(PrimitiveGTO *pgto)
{
    free(pgto);
}

void primitive_gto_normalize(PrimitiveGTO *pgto)
{
    /* Simple formula */
    pgto->norm = sqrt(pow(2, 2 * (pgto->l + pgto->m + pgto->n) + 1.5) *
                      pow(pgto->alpha, pgto->l + pgto->m + pgto->n + 1.5) /
                      fact2(2 * pgto->l - 1) / fact2(2 * pgto->m - 1) /
                      fact2(2 * pgto->n - 1) / pow(M_PI, 1.5));
}

/**
 * Compute the overlap with another gto
 *
 * the overlap function should be copied here
 */
double primitive_gto_overlap(PrimitiveGTO *pgto1, PrimitiveGTO *pgto2)
{
    double N;
    N = pgto1->norm * pgto2->norm;

    return N * overlap(pgto1->alpha,
                       pgto1->l, pgto1->m, pgto1->n,
                       pgto1->x0, pgto1->y0, pgto1->z0,
                       pgto2->alpha,
                       pgto2->l, pgto2->m, pgto2->n,
                       pgto2->x0, pgto2->y0, pgto2->z0
                      );
}


void primitive_gto_recenter(PrimitiveGTO *pgto, double x, double y, double z)
{
    pgto->x0 = x;
    pgto->y0 = y;
    pgto->z0 = z;

    primitive_gto_normalize(pgto);
}

void primitive_gto_set_powers(PrimitiveGTO * pgto, int l, int m, int n)
{
    pgto->l = l;
    pgto->m = m;
    pgto->n = n;
}

double primitive_gto_amp(PrimitiveGTO *pgto, double x, double y, double z)
{
    return  pgto->norm *
            pow(x - pgto->x0, pgto->l) *
            pow(y - pgto->y0, pgto->m) *
            pow(z - pgto->z0, pgto->n) *
            exp(- pgto->alpha * dist2(x, y, z, pgto->x0, pgto->y0, pgto->z0));
}

/*
 * Removing it because this function is accomplished by other libraries
 */

/*
double primitive_gto_coulomb(PrimitiveGTO *pgto1, PrimitiveGTO *pgto2, PrimitiveGTO *pgto3, PrimitiveGTO *pgto4)
{
    return coulomb_repulsion(pgto1->x0, pgto1->y0, pgto1->z0, pgto1->norm,
                             pgto1->l, pgto1->m, pgto1->n, pgto1->alpha,
                             pgto2->x0, pgto2->y0, pgto2->z0, pgto2->norm,
                             pgto2->l, pgto2->m, pgto2->n, pgto2->alpha,
                             pgto3->x0, pgto3->y0, pgto3->z0, pgto3->norm,
                             pgto3->l, pgto3->m, pgto3->n, pgto3->alpha,
                             pgto4->x0, pgto4->y0, pgto4->z0, pgto4->norm,
                             pgto4->l, pgto4->m, pgto4->n, pgto4->alpha
                            );
			    }
*/

int primitive_gto_angular_momentum(PrimitiveGTO *pgto)
{
    return pgto->l + pgto->m + pgto->n;
}

void primitive_gto_R(PrimitiveGTO *pgto, double output[3])
{
    output[0] = pgto->x0;
    output[1] = pgto->y0;
    output[2] = pgto->z0;
}

/**********************************************************************
 * Other internal functions used but pending for removal
 **********************************************************************/

double overlap_1D(int l1, int l2, double PAx,
                  double PBx, double gamma)
{
    /*Taken from THO eq. 2.12*/
    int i;
    double sum;
    sum = 0.;
    for (i = 0; i < (1 + floor(0.5 *(l1 + l2))); i++)
        sum += binomial_prefactor(2 * i, l1, l2, PAx, PBx) *
               fact2(2 * i - 1) / pow(2 * gamma, i);
    return sum;
}

double product_center_1D(double alphaa, double xa,
                         double alphab, double xb)
{
    return (alphaa * xa + alphab * xb) / (alphaa + alphab);
}

double overlap(double alpha1, int l1, int m1, int n1,
               double xa, double ya, double za,
               double alpha2, int l2, int m2, int n2,
               double xb, double yb, double zb)
{
    /*Taken from THO eq. 2.12*/
    double rab2, gamma, xp, yp, zp, pre, wx, wy, wz;

    rab2 = dist2(xa, ya, za, xb, yb, zb);
    gamma = alpha1 + alpha2;
    xp = product_center_1D(alpha1, xa, alpha2, xb);
    yp = product_center_1D(alpha1, ya, alpha2, yb);
    zp = product_center_1D(alpha1, za, alpha2, zb);

    pre = pow(M_PI / gamma, 1.5) * exp(-alpha1 * alpha2 * rab2 / gamma);

    wx = overlap_1D(l1, l2, xp - xa, xp - xb, gamma);
    wy = overlap_1D(m1, m2, yp - ya, yp - yb, gamma);
    wz = overlap_1D(n1, n2, zp - za, zp - zb, gamma);
    return pre * wx * wy * wz;
}




