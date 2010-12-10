/**
 * @file shell.c
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
 */

#include "shell.h"
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include <utils/math.h>

/**
 * Instantiate the Shell struct
 */
Shell *shell_new(void)
{
  return malloc(sizeof(Shell));
}

/**
 * Initializate the struct
 * 
 * - shell->functions
 * - shell->basis_index
 *
 * and initializing the nfunc and ang_mom attribute properly.
 *
 * The stuff is not copied because ContractedGTO are not strictly
 * related to shells.
 */
void shell_init(Shell *shell, int ang_mom)
{
    int size = (ang_mom + 1) * (ang_mom + 2) / 2; /* Size of the shell from the ang_mom */

    shell->functions = malloc(size * sizeof(ContractedGTO *));
    shell->basis_index = malloc(size * sizeof(int));
    shell->ang_mom = ang_mom;
    shell->nfuncs = size;
}


/**
 * Frees the Shell struct deallocing it and the other stuff
 * initialized
 */
void shell_free(Shell *shell)
{
    free(shell->functions);
    free(shell->basis_index);
    free(shell);
}

/**
 * Places the cgto in the shell set of functions, in the correct order.
 * It also store the index in which the cgto reside in the basis set
 */
void shell_append(Shell *shell, ContractedGTO *cgto, int index)
{
    int nx, ny, nz;
    int i, j, ind = 0;

    contracted_gto_R(cgto, shell->r);

    /* Canonical libint order */
    for (i = 0; i <= shell->ang_mom ; i++)
    {
        nx = shell->ang_mom - i;
        for (j = 0; j <= i ; j++)
        {
            ny = i - j;
            nz = j;
            if (nx == cgto->l && ny == cgto->m && nz == cgto->n)
            {
                shell->functions[ind] = cgto;
                shell->basis_index[ind] = index;
                return;
            }
            ind++;
        }
    }
}


int shell_max_num_prim(Shell *shell)
{
    int i;
    int ret;
    ContractedGTO *f;

    ret = 0;

    for (i = 0; i < shell->nfuncs; i++)
    {
        f = shell->functions[i];
        if (f->nprims > ret)
            ret = f->nprims;
    }
    return ret;
}

