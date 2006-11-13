void copy(int N, double *A, int Astride, double *B, int Bstride);
double dot(int N, double *A, int Astride, double *B, int Bstride);
void matmult(int n, int m, int l, double *a, double *b, double *c);
void matmult_tn(int n, int m, int l, double *a, double *b, double *c);
void printmat(int N, int M, double *A, char *title);
void zero(int N, double *A);
void eye(int N, double *A);
void transpose(int n, double *a);

int jacobi(double *a, int n, double d[], double *v, double *b, double *z);
int gjacobi(double *a, double *s, int n, double d[], double *v, 
	    double *b, double *z);

