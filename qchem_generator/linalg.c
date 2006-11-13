/* linalg.c - Simple linear algebra routines in C 
 * This program is part of the PyQuante Program
 * Copyright (c) 2006, Richard P. Muller. All Rights Reserved.
 *
 * Note: the matrices here are stored as single strips of continuous
 *  memory, in "column major" format.
 *
 * One of the goals is to keep the code trivially simple,
 * so that I can implement the same code in f77 and C
 * (and maybe C++ and Python, someday)
 * 
 * This means passing in all workspace! No mallocs here!
 *
 * Obviously, there are close connections to LAPACK here.
 * One major difference is that the matrices here are all
 * expected to be exactly the proper length.
 */

#include <math.h>
#include <stdio.h> /* printf */
#include <stdlib.h>

void copy(int N, double *A, int Astride, double *B, int Bstride){
  /* Copy elements of A into B
   * N       Length of the vectors
   * A       First vector
   * Astride Stride between elements of A vector
   * B       Second vector
   * Bstride Stride between elements of B vector
   */
  int i;
  for (i=0; i<N; i++) B[i*Bstride] = A[i*Astride];
}

double dot(int N, double *A, int Astride, double *B, int Bstride){
  /* Compute the dot product of two vectors 
   * N       Length of the vectors
   * A       First vector
   * Astride Stride between elements of A vector
   * B       Second vector
   * Bstride Stride between elements of B vector
   */
  int i;
  double val=0;
  for (i=0; i<N; i++) val += A[i*Astride]*B[i*Bstride];
  return val;
}

void matmult(int N, int M, int L, double *A, double *B, double *C){
  /* Matrix multiply C = A*B
   * A is a NxL matrix
   * B is a LxM matrix
   * C is a NxM matrix
   */
  int i,j;
  for (i=0; i<N; i++){
    for (j=0; j<M; j++)
      C[i+j*N] = dot(L,&A[i],N,&B[j*L],1);
  }
}

void matmult_tn(int N, int M, int L, double *A, double *B, double *C){
  /* Matrix multiply C = At*B
   * A is a LxN matrix
   * B is a LxM matrix
   * C is a NxM matrix
   */
  int i,j;
  for (i=0; i<N; i++){
    for (j=0; j<M; j++){
      C[i+j*N] = dot(L,&A[i*L],1,&B[j*L],1);
    }
  }
}

/* also need matmult_nt, matmult_tt */

void printmat(int N, int M, double *A, char *title){
  /* Print NxM matrix A */
  int i,j;
  printf("%s\n",title);
  for (i=0; i<N; i++){
    for (j=0; j<M; j++)
      printf("%8.2f ",A[i+j*N]);
    printf("\n");
  }
}

void zero(int N, double *A){
  /* Zero a vector of length N */
  int i;
  for (i=0; i<N; i++) A[i] = 0.;
}

void eye(int N, double *A){
  /* Initialize a NxN matrix A to the identity */
  int i;
  zero(N*N,A);
  for (i=0; i<N; i++) A[i+i*N] = 1.;
}

void transpose(int n, double *a){
  int i,j;
  double tmp;
  for (i=0; i<n; i++){
    for (j=0; j<=i; j++){
      tmp = a[i+j*n];
      a[i+j*n] = a[j+i*n];
      a[j+i*n] = tmp;
    }
  }
}

void evsort(int n, double *d, double *v){
  /*  Sort the eigenvalues in ascending order, and sort the
   *  eigenvalues accordingly.
   *  This uses the awful "buble sort" routine.
   */
  int max_sweeps,isweep,nswap,i,j;
  double tmp;

  if (n > 30)
    max_sweeps = n;
  else
    max_sweeps = 30;

  for (isweep=0; isweep<max_sweeps; isweep++){
    nswap = 0;
    for (i=0; i<n-1; i++){
      if (d[i]>d[i+1]){
	nswap++;
	tmp = d[i];
	d[i] = d[i+1];
	d[i+1] = tmp;
	for (j=0; j<n; j++){
	  tmp = v[j+i*n];
	  v[j+i*n] = v[j+(i+1)*n];
	  v[j+(i+1)*n] = tmp;
	}
      }
    }
    if (nswap == 0) break;
  }
  if (isweep == max_sweeps-1) printf("Warning: max_sweeps exceeded in evsort\n");
  return;
}


int jacobi(double *a, int n, double d[], double *v, double *b, double *z){
  int j,iq,ip,i,nrot=0;
  double tresh,theta,tau,t,sm,s,h,g,c;

  eye(n,v);
  zero(n,z);
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip+ip*n];
  }
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) sm += fabs(a[iq+ip*n]);
    }
    if (sm == 0.0) {
      transpose(n,v);
      evsort(n,d,v);
      return nrot;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	g=100.0*fabs(a[iq+ip*n]);
	if (i > 4 && (fabs(d[ip])+g) == fabs(d[ip])
	    && (fabs(d[iq])+g) == fabs(d[iq]))
	  a[iq+ip*n]=0.0;
	else if (fabs(a[iq+ip*n]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((fabs(h)+g) == fabs(h))
	    t=(a[iq+ip*n])/h;
	  else {
	    theta=0.5*h/(a[iq+ip*n]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[iq+ip*n];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[iq+ip*n]=0.0;
	  for (j=0;j<ip;j++) {
	    g=a[ip+j*n];
	    h=a[iq+j*n];
	    a[ip+j*n]=g-s*(h+g*tau); 
	    a[iq+j*n]=h+s*(g-h*tau);
	  }
	  for (j=ip+1;j<iq;j++) {
	    g=a[j+ip*n];
	    h=a[iq+j*n];
	    a[j+ip*n]=g-s*(h+g*tau); 
	    a[iq+j*n]=h+s*(g-h*tau);
	  }
	  for (j=iq+1;j<n;j++) {
	    g=a[j+ip*n];
	    h=a[j+iq*n];
	    a[j+ip*n]=g-s*(h+g*tau); 
	    a[j+iq*n]=h+s*(g-h*tau);
	  }
	  for (j=0;j<n;j++) {
	    g=v[ip+j*n];
	    h=v[iq+j*n];
	    v[ip+j*n]=g-s*(h+g*tau); 
	    v[iq+j*n]=h+s*(g-h*tau);
	  }
	  nrot++;
	}
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  printf("Too many iterations in routine jacobi");
  exit(1);
}

int gjacobi(double *a, double *s, int n, double d[], 
	    double *v, double *b, double *z){
  /* Use a Jacobi method and a canonical orthogonalization to solve
   * the generalized eigenproblem Ac = ScE
   */

  int info,i,j;

  /* Form the eigen val/vec of S */
  info = jacobi(s,n,d,v,b,z);
  /* For the orthogonalizing matrix X; overwrite S */
  for (j=0; j<n; j++){
    for (i=0; i<n; i++) s[i+j*n] = v[i+j*n]/sqrt(d[j]);
  }

  /* Similarity transform A -> XtAX; use V as temp and overwrite A */
  matmult(n,n,n,a,s,v);
  matmult_tn(n,n,n,s,v,a);

  /* Form the eigen val/vec of XtAX */
  info = jacobi(a,n,d,v,b,z);
  /* overwrite v with xv */
  copy(n*n,v,1,a,1);
  matmult(n,n,n,s,a,v);
  return info;
}  

