/* This is a template file that is read in by make_hf_driver.py and
 * populated with data. You can write it pretty much like normal C
 * code, except when you see a percent sign. If you're printing out
 * data from the python driver, you need to use the Python "dictionary
 * style" of expanding strings. If you're putting in a print statement
 * to be executed by the C program, you need to double the percent
 * signs, because Python is reading and then rewriting the text.
 *
 * This program is part of PyQuante, and is distributed under a
 * modified BSD license. Copyright (c) 2006, Richard P. Muller.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "cints.h"
#include "qc.h"
#include "linalg.h"

int main(){
  /* Simple example program for using cints.c routines */

  int i=0;
  int nat = 3;
  int nbf = 7;
  int nocc=5,info;
  int iter,maxiter=5;
  double h[49],S[49],T[49],V[49],Ints[406];
  double D[49],J[49],K[49];
  double orbe[7],orbs[49],temp[49];
  double temp1[7], temp2[7];
  double H[49];
  double energy,eone,ej,ek;
  double enuke=8.841067;

  /*  Data over atoms  */
  int atno[] = {
    8,
    1,
    1,
  };
  double x[] = {
    0.0,
    1.88971616463,
    0.0,
  };
  double y[] = {
    0.0,
    0.0,
    1.88971616463,
  };
  double z[] = {
    0.0,
    0.0,
    0.0,
  };

  /*  Data over contracted bfns  */
  double xcenter[] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    1.88971616463,
    0.0,
  };  
  double ycenter[] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    1.88971616463,
  };  
  double zcenter[] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
  };  
  int lpower[] = {
    0,
    0,
    1,
    0,
    0,
    0,
    0,
  };
  int mpower[] = {
    0,
    0,
    0,
    1,
    0,
    0,
    0,
  };
  int npower[] = {
    0,
    0,
    0,
    0,
    1,
    0,
    0,
  };
  double normc[] = {
    0.99999971621,
    1.0000001911,
    1.00000019021,
    1.00000019021,
    1.00000019021,
    0.999999946968,
    0.999999946968,
  };
  int istart[] = {
    0,
    3,
    6,
    9,
    12,
    15,
    18,
  };
  int nprim[] = {
    3,
    3,
    3,
    3,
    3,
    3,
    3,
  };

  /*  Data over primitive bfns  */
  double alpha[] = {
    130.709321,
    23.808866,
    6.443608,
    5.033151,
    1.169596,
    0.380389,
    5.033151,
    1.169596,
    0.380389,
    5.033151,
    1.169596,
    0.380389,
    5.033151,
    1.169596,
    0.380389,
    3.425251,
    0.623914,
    0.168855,
    3.425251,
    0.623914,
    0.168855,
  };
  double coef[] = {
    0.154329,
    0.535328,
    0.444635,
    -0.099967,
    0.399513,
    0.700115,
    0.155916,
    0.607684,
    0.391957,
    0.155916,
    0.607684,
    0.391957,
    0.155916,
    0.607684,
    0.391957,
    0.154329,
    0.535328,
    0.444635,
    0.154329,
    0.535328,
    0.444635,
  };
  double normp[] = {
    27.5511677588,
    7.68181997711,
    2.88241776816,
    2.39491476866,
    0.801561774379,
    0.345208161164,
    10.7458317829,
    1.73374383911,
    0.425818989418,
    10.7458317829,
    1.73374383911,
    0.425818989418,
    10.7458317829,
    1.73374383911,
    0.425818989418,
    1.79444186758,
    0.500326654719,
    0.187735124967,
    1.79444186758,
    0.500326654719,
    0.187735124967,
  };

  getS(nbf,nprim,normc,istart,xcenter,ycenter,zcenter, 
       lpower,mpower,npower,normp,coef,alpha,S);
  getT(nbf,nprim,normc,istart,xcenter,ycenter,zcenter, 
       lpower,mpower,npower,normp,coef,alpha,T);
  getV(nbf,nprim,normc,istart,xcenter,ycenter,zcenter, 
       lpower,mpower,npower,normp,coef,alpha,
       nat,atno,x,y,z,V);

  getInts(nbf,nprim,normc,istart,xcenter,ycenter,zcenter, 
	  lpower,mpower,npower,normp,coef,alpha,Ints);

  /* Form h */
  for (i=0; i<nbf*nbf; i++) h[i] = T[i] + V[i];

  /* Form the density matrix from the guess orbs */
  getD(nbf,nocc,orbs,D);

  /* diagonalize h to get guess orbs */
  for (i=0; i<nbf*nbf; i++) {
    orbs[i] = h[i];
    H[i] = h[i];
    temp[i] = S[i];
  }
  info = gjacobi(H,temp,nbf,orbe,orbs,temp1,temp2);

  for (iter=0; iter<maxiter; iter++){
    /* Form the density matrix from the guess orbs */
    getD(nbf,nocc,orbs,D);

    /* Form the J and K matrices */
    getJ(nbf,D,Ints,J);
    getK(nbf,D,Ints,K);

    /* Form the Fock matrix in the orbs space */
    for (i=0; i<nbf*nbf; i++) {
      H[i] = h[i] + 2*J[i] - K[i];
      temp[i] = S[i];
    };
    /* Diagonalize F */
    info = gjacobi(H,temp,nbf,orbe,orbs,temp1,temp2);

    /* Compute the HF energy */
    eone = 0.;
    for (i=0; i<nbf*nbf; i++) eone += h[i]*D[i];
    ej = 0.;
    for (i=0; i<nbf*nbf; i++) ej += D[i]*J[i];
    ek = 0.;
    for (i=0; i<nbf*nbf; i++) ek += D[i]*K[i];
    energy = enuke+2*eone+2*ej-ek;
    printf("%d %f %f %f %f %f\n",iter,energy,enuke,eone,ej,ek);

  }
  return 0;
}
