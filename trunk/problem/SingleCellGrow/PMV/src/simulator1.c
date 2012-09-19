#include <R.h>
#include <Rmath.h>


void getrate(int *x, double *c, double *a){
  a[0] = c[0]*x[0];
  a[1] = c[1]*x[1];
  a[2] = c[2]*x[0];
  a[3] = c[3]*x[1];
  a[4] = c[4]*x[0];
}

void simulator1(double *Tin, double *l1in, double *l2in, double *d1in,
		double *d2in, double *vin, double *win, int *x){
  double T = Tin[0];
  double l1 = l1in[0];
  double l2 = l2in[0];
  double d1 = d1in[0];
  double d2 = d2in[0];
  double v = vin[0];
  double w = win[0];  

  int nr = 5;
  int ns = 2;
  double c[5] = {l1, l2, d1, d2, v};
  double a[5];
  double ca[5];
  int nu[5][2] = {{1,0}, {0,1}, {-1,0}, {0,-1}, {-1,1}};
  
  int rindex, j;
  double u;
  double t_, tau;

  GetRNGstate();

  /* determing initial cell */
  u = unif_rand();
  if ( u < w ) {
    x[0] = 1;
  }
  else {
    x[1] = 1;
  }

  /* Begin main loop */
  t_ = 0.0;
  while ( TRUE ){
    getrate(x, c, a);
    ca[0] = a[0];

    for (j=1; j<nr; j++)
      ca[j] = ca[j-1] + a[j];

    if ( ca[4] == 0 ) break;

    tau = exp_rand();
    tau = tau/ca[4];
    t_ = t_ + tau;
    
    if ( t_ > T ) {
      break;
    }
    
    u = unif_rand();
    u = u*ca[4];
    rindex = -1;
    for ( j=0; j<5; j++ ){
      if ( u < ca[j] ) {
        rindex = j;
	break;
      }
    }
    x[0] = x[0] + nu[rindex][0];
    x[1] = x[1] + nu[rindex][1];
  }

  PutRNGstate();
}

