#include <R.h>
#include <Rmath.h>


void getrate(int *x, double *c, double *a){
  a[0] = c[0]*x[0];
  a[1] = c[1]*x[1];
  a[2] = c[2]*x[0];
  a[3] = c[3]*x[1];
  a[4] = c[4]*x[0];
}

void simulator1p(double *Tin, double *l1in, double *l2in, double *d1in,
	   double *d2in, double *vin, int *x){
  double T = Tin[0];
  double l1 = l1in[0];
  double l2 = l2in[0];
  double d1 = d1in[0];
  double d2 = d2in[0];
  double v = vin[0];
  //double w = win[0];  
  
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
  /*u = unif_rand();
  if ( u < w ) {
    x[0] = 1;
  }
  else {
    x[1] = 1;
  }
  */

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

void simulator1s(double *Tin, double *l1in, double *l2in, double *d1in,
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
    x[1] = 0;
  }
  else {
    x[1] = 1;
    x[0] = 0;
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

void simulator1n(double *Tin, double *l1in, double *l2in, double *d1in,
		 double *d2in, double *vin, double *win, int *Nin, 
		 int *xout){
  double T = Tin[0];
  double l1 = l1in[0];
  double l2 = l2in[0];
  double d1 = d1in[0];
  double d2 = d2in[0];
  double v = vin[0];
  double w = win[0];  
  int N = Nin[0];

  int i;
  int y[2];

  for (i=0; i<N; i++){
    simulator1s(&T, &l1, &l2, &d1, &d2, &v, &w, &y[0]);
    xout[i] = y[0] + y[1];
      }
}

void gridopslow(double *Tin, double *l2truein, double *d2truein, 
		double *l2in, double *d2in, 
		int *msin, int *nsin, double *xout){
  double T = Tin[0];
  double l2true = l2truein[0];
  double d2true = d2truein[0];
  double l2min = l2in[0];
  double l2max = l2in[1];
  double l2step = l2in[2];
  double d2min = d2in[0];
  double d2max = d2in[1];
  double d2step = d2in[2];

  int ms = msin[0];
  int ns = nsin[0];

  int i, j;
  int sp1[ms];
  int sp2[ns];
  double l2, d2;
  double l1, d1;
  double v, w;

  int nl2, nd2;

  double res;

  xout[0] = -1.0;
  xout[1] = -1.0;

  l1 = 0.0;
  d1 = 0.0;
  v = 0.0;
  w = 0.0;
  
  nl2 = (int)((l2max - l2min)/l2step) + 1;
  nd2 = (int)((d2max - d2min)/d2step) + 1;

  simulator1n(&T, &l1, &l2true, &d1, &d2true, &v, &w, &ms, sp1);
  printf("[%f %f] \n", l2true, d2true);
  for(i=0; i<ms; i++) printf("%4d ", sp1[i]);

  printf("\n");

  for (i=0; i<=nl2; i++){
    l2 = l2min + i*l2step;
    for (j=0; j<=nd2; j++){
      d2 = d2min + j*d2step;
      simulator1n(&T, &l1, &l2, &d1, &d2, &v, &w, &ns, sp2);      
      res = kstest(sp1, ms, sp2, ns);
      printf("--[%f %f : %f] \n", l2, d2, res);
    }
  }
}

