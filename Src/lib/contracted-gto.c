/**
 * @file contracted-gto.c
 * 
 * This source code is part of the PyQuante Quantum Chemistry suite.
 * 
 * Most function basically sums/products over the primitive values.
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

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "contracted-gto.h"
#include "primitive-gto.h"

/**
 * Allocated also the first element of primitives to permit further
 * realloc'ing
 */
ContractedGTO * contracted_gto_new(void)
{
    ContractedGTO *cgto;

    cgto = malloc(sizeof(ContractedGTO));
    cgto->primitives = malloc(sizeof(PrimitiveGTO *));

    return cgto;
}

void contracted_gto_init(ContractedGTO *cgto,  double x, double y, double z, int l, int m, int n, int atid)
{
    cgto->l = l;
    cgto->m = m;
    cgto->n = n;
    cgto->x0 = x;
    cgto->y0 = y;
    cgto->z0 = z;

    cgto->atid = atid;

    cgto->nprims = 0;
}

void contracted_gto_from_primitives(ContractedGTO *cgto, PrimitiveGTO ** pgtos, int ngtos)
{
  int i;

    for ( i = 0; i < ngtos ; i++) 
      {
	contracted_gto_add_primitive(cgto, pgtos[i], pgtos[i]->coef);
      }

    cgto->nprims = ngtos;

    /* Initializing */
    cgto->norm = 1.0;
    contracted_gto_normalize(cgto);

}

void contracted_gto_add_primitive(ContractedGTO *cgto, PrimitiveGTO *pgto, double coef)
{
    PrimitiveGTO *pgto_copy;
    
    cgto->primitives = realloc(cgto->primitives, sizeof(PrimitiveGTO *) * (cgto->nprims + 1));

    /* Allocating the new pgto */
    pgto_copy = primitive_gto_new();
    
    /* Copying the gto */
    memcpy(pgto_copy, pgto, sizeof(PrimitiveGTO));
    
    /* Adding the coef and packing */
    pgto_copy->coef = coef;
    cgto->primitives[cgto->nprims] = pgto_copy;
    cgto->nprims += 1;
}

/**
 * Since all the primitives were allocated by the cgto, we have also
 * to deallocate them.
 */
void contracted_gto_free(ContractedGTO *cgto)
{
  int i;
    for ( i = 0; i < cgto->nprims ; i++) 
      {
        primitive_gto_free(cgto->primitives[i]);
      }
    
    free(cgto->primitives);
    free(cgto);
}

void contracted_gto_normalize(ContractedGTO *cgto)
{

    double S;
    /* Normalizing primitives another time ?*/

    S = contracted_gto_overlap(cgto, cgto);
    cgto->norm = cgto->norm / sqrt(S);
}

double contracted_gto_overlap(ContractedGTO *cgto1, ContractedGTO *cgto2)
{
    double S;
    int i, j;

    S = 0.0 ;
    for (i = 0; i < cgto1->nprims ; i++)
    {
        for (j = 0; j < cgto2->nprims ; j++)
        {
            S += cgto1->primitives[i]->coef * cgto2->primitives[j]->coef *
                 primitive_gto_overlap(cgto1->primitives[i],
                                       cgto2->primitives[j]);
        }
    }
    return S * (cgto1->norm) * (cgto2->norm);

}

double contracted_gto_amp(ContractedGTO *cgto, double x, double y, double z)
{
    double amp = 0;
    int i;
    double N = cgto->nprims;
    for (i = 0; i < N; ++i)
    {
        amp += cgto->primitives[i]->coef *
               primitive_gto_amp(cgto->primitives[i], x, y, z);
    }
    return cgto->norm * amp;
}

/*
double contracted_gto_coulomb(ContractedGTO *cgto1, ContractedGTO *cgto2, ContractedGTO *cgto3, ContractedGTO *cgto4)
{
    int i, j, k, l;
    double J = 0., incr = 0.;

    for (i = 0; i < cgto1->nprims; i++)
        for (j = 0; j < cgto2->nprims; j++)
            for (k = 0; k < cgto3->nprims; k++)
                for (l = 0; l < cgto4->nprims; l++)
                {
                    incr = primitive_gto_coulomb(cgto1->primitives[i],
                                                 cgto2->primitives[j],
                                                 cgto3->primitives[k],
                                                 cgto4->primitives[l]);
                    J += incr *
                         cgto1->primitives[i]->coef *
                         cgto2->primitives[j]->coef *
                         cgto3->primitives[k]->coef *
                         cgto4->primitives[l]->coef;
                }
    return J;
}
*/

int contracted_gto_angular_momentum(ContractedGTO *cgto)
{
    return primitive_gto_angular_momentum(cgto->primitives[0]);
}

void contracted_gto_R(ContractedGTO *cgto, double output[3])
{
    primitive_gto_R(cgto->primitives[0], output);
}

/* Auxiliary functions, to Test */

void contracted_gto_recenter(ContractedGTO *cgto, double x, double y, double z)
{
    int i;
    for (i = 0; i < cgto->nprims ; i++)
    {
        primitive_gto_recenter(cgto->primitives[i], x, y, z);
    }
    contracted_gto_normalize(cgto);
}

void contracted_gto_set_powers(ContractedGTO *cgto, int l, int m, int n)
{
    int i;
    cgto->l = l;
    cgto->m = m;
    cgto->n = n;

    for (i = 0; i < cgto->nprims ; i++)
    {
        primitive_gto_set_powers(cgto->primitives[i], l, m, n);
    }
    contracted_gto_normalize(cgto);
}


/**
 * Libint related, part of contracted_gto_renorm_prefactor
 */
double contracted_gto_libint_renorm(ContractedGTO *cgto)
{
    double renorm;
    int oldl, oldm, oldn;

    renorm = cgto->norm;

    oldl = cgto->primitives[0]->l;
    oldm = cgto->primitives[0]->m;
    oldn = cgto->primitives[0]->n;

    contracted_gto_set_powers(cgto ,
                              cgto->l + cgto->m + cgto->n, 0, 0);
    renorm /= cgto->norm;

    /* Put back the l,y,z */
    contracted_gto_set_powers(cgto, oldl, oldm, oldn);

    return renorm;
}

/**
 * Libint related
 */
double contracted_gto_renorm_prefactor(ContractedGTO *cgto1, ContractedGTO *cgto2,
                                       ContractedGTO *cgto3, ContractedGTO *cgto4)
{
    return contracted_gto_libint_renorm(cgto1) *
           contracted_gto_libint_renorm(cgto2) *
           contracted_gto_libint_renorm(cgto3) *
           contracted_gto_libint_renorm(cgto4);
}

