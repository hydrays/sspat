/* eigensub.h
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define square(a) ((a)*(a))
#define FOR(i,n) for(i=0; i<n; i++)
#define FPN(file) fputc('\n', file)
#define F0 stdout
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))


double sum (double x[], int n);
int xtoy (double x[], double y[], int n);
int abyx (double a, double x[], int n);
double distance (double x[], double y[], int n);
int matinv (double x[], int n, int m, double space[]);
void SetSeed (int seed);
double rndu (void);

int matby (double a[], double b[], double c[], int n,int m,int k);
int matout (FILE * file, double x[], int n, int m);
int eigen (int job, double A[], int n, double rr[], double ri[],
          double vr[], double vi[], double w[]);

typedef struct { double re, im; } complex;
#define csize(a) (fabs(a.re)+fabs(a.im))

complex compl (double re,double im);
complex conj (complex a);
complex cplus (complex a, complex b);
complex cminus (complex a, complex b);
complex cby (complex a, complex b);
complex cdiv (complex a,complex b);
complex cexp (complex a);
complex cfactor (complex x, double a);
int cxtoy (complex x[], complex y[], int n);
int cmatby (complex a[], complex b[], complex c[], int n,int m,int k);
int cmatout (FILE * fout, complex x[], int n, int m);
int cmatinv( complex x[], int n, int m, double space[]);

