/**
 * @file clibint.c
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




#include "shell.h"
#include <utils/swap.h>
#include <libint.h>
#include "clibint.h"

#include <stdlib.h>

#include <math.h>
#include <utils/math.h>

/**
 * The swapping is done in reverse order to return to the original
 * state.
 */
int swapped_ijkl(int i, int j, int k, int l,
                 int ni, int nj, int nk, int nl, int swapflag)
{
    if (swapflag & ABCDSWAP)
    {
        swap(&i, &k, sizeof(int));
        swap(&ni, &nk, sizeof(int));
        swap(&j, &l, sizeof(int));
        swap(&nj, &nl, sizeof(int));
    }
    if (swapflag & CDSWAP)
    {
        swap(&k, &l, sizeof(int));
        swap(&nk, &nl, sizeof(int));
    }
    if (swapflag & ABSWAP)
    {
        swap(&i, &j, sizeof(int));
        swap(&ni, &nj, sizeof(int));
    }

    return ((i * nj + j) * nk + k) * nl + l;
}

void shell_libint_renorm(Shell *shell, double *output)
{
    int i;
    for (i = 0; i < shell->nfuncs; i++)
    {
        output[i] = contracted_gto_libint_renorm(shell->functions[i]) ;
    }

}

/**
 * This function must accomplish some tasks, this is not cohesive, may
 * be refactored.
 *
 * It allocate and initialize the necessary stuff, fills the intobj
 * with the necessary informations, calculate the ERIs and put them in
 * the output buffer in the correct order.
 */
void shell_compute_eri(Shell *shell1, Shell *shell2, Shell *shell3, Shell *shell4, double *output)
{

    Libint_t *intobj; /* Integral object to pass to libint */
    int max_am, max_num_prim_comb; /* max_am = max angular momentum
				   max_num_prim_comb = max number of contractions^4 */

    double AB[3], CD[3];		/* Geometric quantities */
    prim_data pdata;  		/* Support variable that represents
				   primitive data */

    ContractedGTO *cgto_i, *cgto_j, *cgto_k, *cgto_l;
    int ijkl, ijklout;
    int size;
    int primitive_number, buffer_offset;

    double *integrals;
    double integral;
    int i, j, k, l;
    int p, q, r, s;
    int abswap = 0, cdswap = 0, abcdswap = 0;

    /*
     * Initializing intobj
     */
    intobj = malloc(sizeof(Libint_t));

    max_am = max4(shell1->ang_mom, shell2->ang_mom,
                  shell3->ang_mom, shell4->ang_mom);
    max_num_prim_comb = pow(max4(shell_max_num_prim(shell1),
                                 shell_max_num_prim(shell2),
                                 shell_max_num_prim(shell3),
                                 shell_max_num_prim(shell4)), 4);

    init_libint(intobj, max_am, max_num_prim_comb);

    /*
     * Preparing the data
     */
    /* Swapping shells if necessary */
    if (shell1->ang_mom < shell2->ang_mom)
    {
        swap(&shell1, &shell2, sizeof(Shell *));
        abswap = 1;
    }
    if (shell3->ang_mom < shell4->ang_mom)
    {
        swap(&shell3, &shell4, sizeof(Shell *));
        cdswap = 1;
    }
    if ((shell1->ang_mom + shell2->ang_mom) > (shell3->ang_mom + shell4->ang_mom))
    {
        swap(&shell1, &shell3, sizeof(Shell *));
        swap(&shell2, &shell4, sizeof(Shell *));
        abcdswap = 1;
    }

    /* Computing AB CD */
    vec_subtract(shell1->r, shell2->r, AB);
    vec_subtract(shell3->r, shell4->r, CD);
    for (i = 0; i < 3 ; i++)
    {
        intobj->AB[i] = AB[i];
        intobj->CD[i] = CD[i];
    }

    /* Size of the return array, combination of shells */
    size = shell1->nfuncs * shell2->nfuncs *
           shell3->nfuncs * shell4->nfuncs;

    buffer_offset = 0;

    primitive_number = 0; /* Total number of primitive */
    i = j = k = l = 0;

    /* only the first function as a sample */
    cgto_i = shell1->functions[i];
    cgto_j = shell2->functions[j];
    cgto_k = shell3->functions[k];
    cgto_l = shell4->functions[l];





    for (p = 0; p < cgto_i->nprims ; p++)
    {
        for (q = 0; q < cgto_j->nprims ; q++)
        {
            for (r = 0; r < cgto_k->nprims; r++)
            {
                for (s = 0; s < cgto_l->nprims; s++)
                {
                    /* Compute primitive data for Libint */
                    pdata = compute_primitive_data(cgto_i->primitives[p],
                                                   cgto_j->primitives[q],
                                                   cgto_k->primitives[r],
                                                   cgto_l->primitives[s]);

                    intobj->PrimQuartet[primitive_number++] = pdata;
                }
            }
        }
    }

    /* Compute the integrals */
    if (shell1->ang_mom == 0 && shell2->ang_mom == 0
            && shell3->ang_mom == 0 && shell4->ang_mom == 0)
    {
        /* Compute the F[0] stuff */
        integral = 0;
        for (i = 0; i < primitive_number ; i++)
        {
            integral += intobj->PrimQuartet[i].F[0];
        }
        output[0] = integral;
    }
    else
    {

        integrals = build_eri
                    [shell1->ang_mom]
                    [shell2->ang_mom]
                    [shell3->ang_mom]
                    [shell4->ang_mom](intobj, primitive_number);

        double shell1_renorm[shell1->nfuncs];
        double shell2_renorm[shell2->nfuncs];
        double shell3_renorm[shell3->nfuncs];
        double shell4_renorm[shell4->nfuncs];

        /* Computing the renormalization prefactors */
        shell_libint_renorm(shell1, shell1_renorm);
        shell_libint_renorm(shell2, shell2_renorm);
        shell_libint_renorm(shell3, shell3_renorm);
        shell_libint_renorm(shell4, shell4_renorm);

        /* Retrieve the primitive integrals and form the contracted-integral value*/
        for (i = 0; i < shell1->nfuncs ; i++)
        {
            for (j = 0; j < shell2->nfuncs ; j++)
            {
                for (k = 0; k < shell3->nfuncs ; k++)
                {
                    for (l = 0; l < shell4->nfuncs ; l++)
                    {

                        ijklout = swapped_ijkl(i, j, k, l,
                                               shell1->nfuncs,
                                               shell2->nfuncs,
                                               shell3->nfuncs,
                                               shell4->nfuncs,
                                               abswap | (2 * cdswap) | (4 * abcdswap));

                        ijkl = ((i * shell2->nfuncs + j) * shell3->nfuncs + k) * shell4->nfuncs + l;
                        output[ijklout] = integrals[ijkl]  * (shell1_renorm[i] * shell2_renorm[j] *
                                                              shell3_renorm[k] * shell4_renorm[l]);
                    }
                }
            }
        }
    }

    free_libint(intobj);
    free(intobj);
}

/**
 * Computes the data to libint, it is somewhat simple.
 */
prim_data compute_primitive_data(PrimitiveGTO *a, PrimitiveGTO *b,
                                 PrimitiveGTO *c, PrimitiveGTO *d)
{
    prim_data pdata;
    double zeta, eta, rho, S12, S34;
    double Ca, Cb, Cc, Cd;
    double A[3], B[3], C[3], D[3];
    double P[3], Q[3], W[3];
    int m;
    int i;

    primitive_gto_R(a, A);
    primitive_gto_R(b, B);
    primitive_gto_R(c, C);
    primitive_gto_R(d, D);

    zeta = a->alpha + b->alpha;
    eta = c->alpha + d->alpha;
    rho = zeta * eta / (zeta + eta);

    for (i = 0; i < 3 ; i++)
    {
        P[i] = (A[i] * a->alpha + B[i] * b->alpha) / zeta;
        Q[i] = (C[i] * c->alpha + D[i] * d->alpha) / eta;
        W[i] = (P[i] * zeta + Q[i] * eta) / (zeta + eta);
    }

    /* Coefficinet including normalization */
    Ca = a->coef * a->norm;
    Cb = b->coef * b->norm;
    Cc = c->coef * c->norm;
    Cd = d->coef * d->norm;

    /* Overlap ? Normalization costant*/
    S12 = pow(M_PI / zeta, 1.5) * exp(- a->alpha * b->alpha / zeta * vec_dist2(A, B));
    S34 = pow(M_PI / eta, 1.5) * exp(- c->alpha * d->alpha / eta * vec_dist2(C, D));

    pdata.twozeta_a = 2 * a->alpha;
    pdata.twozeta_b = 2 * b->alpha;
    pdata.twozeta_c = 2 * c->alpha;
    pdata.twozeta_d = 2 * d->alpha;

    pdata.oo2z = 1.0 / (2 * zeta);
    pdata.oo2n = 1.0 / (2 * eta);
    pdata.oo2zn = 1.0 / (2 * (zeta + eta));
    pdata.poz = rho / zeta;
    pdata.pon = rho / eta;
    pdata.oo2p = 1.0 / (2 * rho);

    vec_subtract(P, A,  pdata.U[0]);
    vec_subtract(Q, C,  pdata.U[2]);
    vec_subtract(W, P,  pdata.U[4]);
    vec_subtract(W, Q,  pdata.U[5]);

    m = primitive_gto_angular_momentum(a) +
        primitive_gto_angular_momentum(b) +
        primitive_gto_angular_momentum(c) +
        primitive_gto_angular_momentum(d);

    for (i = 0; i <= m ; i++)
    {
        pdata.F[i] = 2 * Fgamma(i, rho * vec_dist2(P, Q)) * sqrt(rho / M_PI) * S12 * S34 *
                     Ca * Cb * Cc * Cd ;
    }
    return pdata;
}
