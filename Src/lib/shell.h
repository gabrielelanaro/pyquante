/** 
 * @file shell.h
 * @author Gabriele Lanaro <gabriele.lanaro@gmail.com>
 *
 * This source code is part of the PyQuante Quantum Chemistry suite.
 * 
 * Defines the struct Shell and the related functions.
 *
 * Implementation of the shell_compute_eri function is done in the
 * clibint.h. libint library is used to compute integrals in a fast
 * way
 *
 * [http://www.ccmst.gatech.edu/evaleev/libint/] 
 *  
 * Written by Gabriele Lanaro, 2009-2010
 * Copyright (c) 2009-2010, Gabriele Lanaro
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * 
 */

#ifndef _SHELL_H_
#define _SHELL_H_

#include "contracted-gto.h"


/**
 * Represents the Cartesian Shell for a given cartesian number, it is
 * basically a set of ContractedGTOs.
 * 
 * The shell with ang_mom 1 is this set:
 * ContractedGTO -> l,m,n = 1,0,0
 * ContractedGTO -> l,m,n = 0,1,0
 * Contractedgto -> l,m,n = 0,0,1
 *
 * In addition they must stay on the same center and they have to
 * share the same set of primitives.
 */
typedef struct  {
  /**
   * Set of ContractedGTOs
   */
  ContractedGTO **functions;
  /**
   * Set of index from the BasisSet
   * 
   * The BasisSet is defined in python, this attribute is mainly an
   * hack to get the correct order of integrals
   */
  int *basis_index;
  /**
   * Angular momentum of the shell, S->0, P->1, D->2 etc...
   */
  int ang_mom;
  /**
   * Number of functions in the shell.
   * Treated as the length of functions.
   */
  int nfuncs;
  /**
   * Center of the shell
   * TODO: pending from removal, this is redundant and not so useful.
   */
  double r[3];
} Shell ;


/**
 * Instantiate the shell
 */
Shell *shell_new(void);

/**
 * Initializate the shell malloc'ing the necessary space.
 * 
 * @param ang_mom : Angular momentum of the shell
 */
void shell_init(Shell *shell, int ang_mom);

/**
 * Deallocate the shell
 */
void shell_free(Shell *shell);

/**
 * Append a cgto to the shell 
 *
 * @param index : this is the index from the original BasisSet. This
 * is used by python code, it can be safely set to 0.
 */
void shell_append(Shell *shell, ContractedGTO *cgto, int index);

/**
 * Compute the max number of primitives in the shell. This function is
 * redundant since in the shell the cgtos are from the same family.
 *
 * @todo pending for removal!
 */
int shell_max_num_prim(Shell *shell);

#endif /* _SHELL_H_ */
