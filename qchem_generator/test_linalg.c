#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>

int main(){
  double *B, *C;
  int N = 3; double A[] = { 1, 2, 3, 2, 4, 5, 3, 5, 6 };
  //int N = 2; double A[] = {1,2,2,3};
  double val = 0;
  double D[N],TMP1[N],TMP2[N];
  int info,i;

  B = (double *)malloc(N*N*sizeof(double));
  C = (double *)malloc(N*N*sizeof(double));

  eye(N,B);
  eye(N,C);

  printmat(N,N,A,"Matrix A");
  printmat(N,N,B,"Matrix B");

  matmult(N,N,N,A,A,C);
  printmat(N,N,C,"Matrix A*A");

  matmult(N,N,N,A,B,C);
  printmat(N,N,C,"Matrix A*B");

  matmult(N,N,N,B,B,C);
  printmat(N,N,C,"Matrix B*B");

  matmult_tn(N,N,N,A,B,C);
  printmat(N,N,C,"Matrix At*B");

  matmult_tn(N,N,N,A,A,C);
  printmat(N,N,C,"Matrix At*A");

  val = dot(N*N,A,1,A,1);
  printf("Dot product results = %f\n",val);

  info = jacobi(A,N,D,B,TMP1,TMP2);
  for (i=0; i<N; i++) printf("D[%d] = %f\n",i,D[i]);
  printmat(N,N,B,"Evecs(A)");
  
  free(B);
  free(C);

  return 0;
}
