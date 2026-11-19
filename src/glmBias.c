#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/RS.h>
#ifndef FCONE
# define FCONE
#endif

#define LWORK 10000

int * jpvt=NULL;
double * tau=NULL;
int lwork;
double * work=NULL;
int info;
int rank;
double * hatwq=NULL;
double * Q =NULL;
int np;
//double *coef; /* p x 1 vector */


/* BLAS/LAPACK constants */
char * const L="L";
char * const T="T";
char * const U="U";
char * const N="N";
const double dOne=1.0;
const int iOne=1;
const double dZero=0.0;

void getBias(int* n, int* p, double* rtwx, double* wrt, double* coef)
{
	int i;
	
	/* QR decomposition with column pivoting using level 3 BLAS */
	F77_CALL(dgeqp3)(n, p, rtwx, n, jpvt, tau, work, &lwork, &info); 
	/* test error */
	if(info != 0) error("'dgeqp3' QR decomposition error");
	
	/* Find Q matrix */
	F77_CALL(dcopy)(&np, rtwx, &iOne, Q, &iOne);
	F77_CALL(dorgqr)(n, p, &rank, Q, n, tau, work, &lwork, &info);
	
	/* Find diagonal of hat matrix, multiply by inverse sqrt of weight and -0.5 */
	for(i=(*n)-1; i>=0; i--)
		hatwq[i] = -0.5 * pow(F77_CALL(dnrm2)(p, Q+i, n), 2.0) / wrt[i];
	
	/* Find Q'(RHS) */
	F77_CALL(dgemv)(T, n, p, &dOne, Q, n, hatwq, &iOne, &dZero, coef, &iOne FCONE);
	
	/* Solve triangular system for coefficients */
	F77_CALL(dtrsv)(U, N, N, p, rtwx, n, coef, &iOne FCONE FCONE FCONE);
	
	/* Reverse the pivoting */
	for(i=*p-1; i>=0; i--)
		hatwq[jpvt[i]-1]=coef[i];  /* re-using hatwq space */
	F77_CALL(dcopy)(p, hatwq, &iOne, coef, &iOne);
	
}

void initQRdecomp(int* n, int* p)
{
	lwork = LWORK;
	rank=*p; np=(*n) * (*p); 

	if(p) R_Free(jpvt);
	jpvt=R_Calloc(*p, int); 
	
	if(tau) R_Free(tau);
	tau=R_Calloc(*p, double); 
	
	if(work) R_Free(work);
	work=R_Calloc(lwork, double);
	
	if(hatwq) R_Free(hatwq);
	hatwq=R_Calloc(*n, double);
	
	if(Q) R_Free(Q);
	Q=R_Calloc(np, double);
}
void finalQRdecomp()
{
	R_Free(jpvt); jpvt=NULL;
	R_Free(tau);  tau=NULL;
	R_Free(work); work=NULL;
	R_Free(hatwq); hatwq=NULL;
	R_Free(Q); Q=NULL;
}

SEXP getGlmBias(SEXP rtwx, SEXP wrt, SEXP ngood, SEXP rk)
{
	SEXP out;
	out = PROTECT(allocVector(REALSXP, *(INTEGER(rk))));
	getBias(INTEGER(ngood), INTEGER(rk), REAL(rtwx), REAL(wrt), REAL(out));
	UNPROTECT(1);
	return out;
}
