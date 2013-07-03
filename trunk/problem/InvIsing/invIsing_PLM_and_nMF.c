#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <atlas/clapack.h>

#define p 2 /*system size*/
#define connprob 1.0 /*sparsity factor*/
#define n 10000	/*number of samples*/
#define JUMPS 10000 /*sampling period*/
#define WARMUPS 5000000 /*MC warmups*/
#define B 1.0 /*inverse temperature*/
#define TOL 0.000001 
#define SEED 2
#define PI 3.14159265
#define HES 1


long long int samples[n];



inline double unirnd() {
  return (double)(random()+1.0)/(RAND_MAX+1.0);
}


double randn() {
  return sqrt(-2*log(unirnd()))*cos(2*PI*unirnd());
}

main()
{
  int i,j,k,l,r,info,ipiv[p-1],index[p],X[p];
  char states[2]={-1,1};
  double C[p][p]={0},m[p]={0},delta=0,J[p][p],Jest[p][p],Jnaive[p][p];
  double h[p];

  FILE *fpJ;	
  FILE *fp;	

  //srandom(time(NULL));	
  /* srandom(SEED);	 */
  /* for(i=0;i<p;i++) { */
  /*   for(j=i;j<p;j++) { */
  /*     if (i==j) J[i][j]=0; */
  /*     else {	 */
  /* 	if (unirnd()<connprob)  { */
  /* 	  J[i][j]=J[j][i]=1.0/sqrt(p*connprob)*randn(); */
  /* 	} */
  /* 	else J[i][j]=J[j][i]=0; */
				
  /*     } */
  /*   } */

  /* } */

  J[0][0] = 0.0;
  J[1][1] = 0.0;
  J[2][2] = 0.0;
  J[0][1] = 0.0;
  J[1][0] = 0.0;
  //J[0][2] = J[2][0] = 1.0;
  //J[1][2] = J[2][1] = 0.0;
  h[0] = 1.0;
  h[1] = 0.0;
  //h[2] = 0.0;

  for(i=0;i<p;i++) {	
    X[i]=states[random() % 2];	
  }

  for(k=0;k<WARMUPS;k++) {
    for(j=0;j<p;j++) {
      double deltaE;
      deltaE=0;
      for(i=0;i<p;i++) {
	deltaE+=J[i][j]*X[i];	
      }
      deltaE = (deltaE+h[j])*2*B*X[j];
      if (deltaE < 0 || unirnd() < exp(-deltaE)) {
	X[j]=-X[j];
      }
    }
    if ((k % 1000) == 0) printf("%f\n",(100.0*k)/WARMUPS);
  }
  for(k=0;k<n;k++) {
    for(j=0;j<JUMPS;j++) {
      int o = random() % p;
      double deltaE;
      deltaE=0;
      for(i=0;i<p;i++) {
	deltaE+=J[i][o]*X[i];	
      }
      deltaE = (deltaE+h[o])*2*B*X[o];
      //deltaE*=2*B*X[o];
      if (deltaE < 0 || unirnd() < exp(-deltaE)) {
	X[o]=-X[o];
      }
    }
    samples[k]=0;
    for (j=p-1;j>=0;j--) {
      samples[k]<<=1;
      samples[k] |=(X[j]==1);
    }
    if ((k % 1000) == 0) printf("%f\n",(100.0*k)/n);
		
  }	
  printf("\n Done sampling \n");
  
  for(r=0;r<p;r++) {	
    double thetas[p-1];
    int done=0,count=0;
		
    for(j=0;j<p-1;j++) {
      thetas[j]=0.0;				
    }

    while(!done) {
      double f[p-1]={0},JAC[p-1][p-1]={0},gradnorm=0,gradmax=0;
      int rnd=random() % HES;

      for(l=0;l<n;l++) {
	double z=0,A1,A2;			
	for(j=0;j<p-1;j++) {
	  z+=thetas[j]*(((samples[l] >> j+(j>=r))&1)*2-1);				
	}
	z=exp(2*B*(((samples[l] >> r)&1)*2-1)*z);
	A1=-1.0/n*2.0*B*(((samples[l] >> r)&1)*2-1)/(z+1.0);
				
	if ((l % HES) == rnd) {
	  A2=1.0/n*(4.0*B*B*z)/((z+1.0)*(z+1.0));
	  for(j=0;j<p-1;j++) {
	    f[j]+=(((samples[l] >> j+(j>=r))&1)*2-1)*A1;	
	    double A3=(((samples[l] >> j+(j>=r))&1)*2-1)*A2;			
	    for(k=0;k<j+1;k++) {
	      JAC[j][k]+=HES*(((samples[l] >> k+(k>=r))&1)*2-1)*A3;
	    }
	  }
	}
	else {
	  for(j=0;j<p-1;j++) {
	    f[j]+=(((samples[l] >> j+(j>=r))&1)*2-1)*A1;	
	  }
	}
      }
      for(j=0;j<p-1;j++) {
	for(k=j+1;k<p-1;k++) {
	  JAC[j][k]=JAC[k][j];
	}
	gradnorm+=f[j]*f[j];
	if (gradmax<sqrt(f[j]*f[j])) gradmax=sqrt(f[j]*f[j]);
      }
      gradnorm=sqrt(gradnorm);
      printf("%.10f %.10f\n",gradmax,gradnorm);
      info = clapack_dgesv(CblasRowMajor, p-1, 1, &JAC[0][0], p-1, ipiv, &f[0], p-1);
      if (info != 0) fprintf(stderr, "failure with error %d\n", info);


      done=1;
      for(j=0;j<p-1;j++) {
	thetas[j]-=f[j];
	if (f[j]>TOL || f[j]<-TOL) done=0;
      }
			
      count++;
    }
    printf("%d %d\n",count,r);
    for(j=0;j<p-1;j++) {
      Jest[j+(j>=r)][r]=thetas[j];				
    }
    /*put thetas as columns in jest*/
  }
  for(j=0;j<p;j++) {
    Jest[j][j]=0;				
  }




  /*calc delta*/
  for(j=0;j<p;j++) {
    for(i=0;i<j;i++) {
      delta+=((Jest[j][i]+Jest[i][j])/2.0-J[j][i])*((Jest[j][i]+Jest[i][j])/2.0-J[j][i]);
    }
  }
  delta=sqrt(delta/((p-1)/(2.0*connprob)));
  printf("%f\n",delta);



  /*naive part*/
  for(k=0;k<n;k++) {

    for(j=0;j<p;j++) {
      m[j]+=(((samples[k] >> j)&1)*2-1);
      for(i=0;i<p;i++) {
	C[j][i]+=(((samples[k] >> j)&1)*2-1)*(((samples[k] >> i)&1)*2-1);
      }
    }
  }
  for(j=0;j<p;j++) {
    m[j]/=n;
  }
  for(j=0;j<p;j++) {
    for(i=0;i<p;i++) {
      C[j][i]/=n;
      C[j][i]-=m[j]*m[i];
    }
  }

  /* fpJ=fopen("Bconnprob.dat","w"); */
  /* fprintf(fpJ,"%f %f ",B,connprob); */
  /* fclose(fpJ); */

  /* fpJ=fopen("parametersest.dat","w"); */
  /* for(j=0;j<p;j++) { */
  /*   for(i=0;i<p;i++) { */
  /*     fprintf(fpJ,"%.12lf ",Jest[j][i]);		 */
  /*   } */
  /*   fprintf(fpJ,"\n");	 */
				
  /* } */
  /* fclose(fpJ);	 */

  fpJ=fopen("parameters.dat","w");

  for(j=0;j<p;j++) {
    fprintf(fpJ,"%.16lf ",h[j]);
    for(i=0;i<p;i++) {
      fprintf(fpJ,"%.16lf ",J[j][i]);
    }
    fprintf(fpJ,"\n");	
  }
  fclose(fpJ);

  fp=fopen("sample.dat","w");
  /* Output samples to file */
  for (l=0; l<n; l++) {
    for (k=0; k<p; k++) {
      fprintf(fpJ,"%4d", (int)((samples[l] >> k)&1)*2-1);
    }
    fprintf(fp,"\n");	
  }
  fclose(fp);

  /* fp=fopen("m.dat","w"); */
  /* fpJ=fopen("C.dat","w"); */

  /* for(j=0;j<p;j++) { */
  /*   fprintf(fp,"%.16lf\n",m[j]); */
  /*   for(i=0;i<p;i++) { */
  /*     fprintf(fpJ,"%.16lf ",C[j][i]); */
  /*   } */
  /*   fprintf(fpJ,"\n");	 */
  /* } */
  /* fclose(fpJ); */
  /* fclose(fp); */
  /* system("octave egenpnaive.m"); */
	
}	
	

