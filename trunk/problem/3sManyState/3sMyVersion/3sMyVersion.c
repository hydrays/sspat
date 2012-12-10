/* 3s.c
   This is an ML program to estimate theta0, theta1, tau0, and tau1 for the 3-species problem, 
   using numerical integration to calculate the likelihood.
	
   cl -O2 3s.c tools.c
   cl -o 3s -mcpu=G5 -O4 -funroll-loops -fomit-frame-pointer -finline-functions 3s.c tools.c
   cc -o 3s -m64 -march=opteron -mtune=opteron -ansi -O3 -funroll-loops -fomit-frame-pointer -finline-functions 3s.c tools.c -lm
	  
   3s
   3s <ctlfile>
		
   npoints can be 4, 8, 16, 32, 64, 128, ..., 1024.  The default is 32.
*/

#include "paml.h"
#include "eigensub.h"

#define NS            3
#define NGENE         1  /* required by ReadSeq, but not really in use */
#define LSPNAME       50
#define NCODE         4
/*zzz*/
#define DIM           22

extern double Small_Diff;
extern int noisy, NFunCall;

/* zzz: variables from pt*/
double mc1[]=
  {-3,0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, -1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
   0, 0, -1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, -1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

double mc2[]=
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, -1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, -1,0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, -1,0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, -3,0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,0, 0, 0, 0, 0, 0, 0, 0, 1,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,0, 0, 0, 0, 0, 0, 0, 1,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};


double mw1[]=
  {-3,1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0,-2, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0,-2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0,-1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0,-2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0,-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1 , 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

double mw2[]=
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 1, 1,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   1, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 1, 0, 0, 1,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 1, 0, 1, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 1, 0, 1, 1,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 1, 0, 0, 1, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 1, 0, 0, 1, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1};
//double c1=400,c2=400,w=800;
//double Q[DIM*DIM]={0};

/* zzz */

struct CommonInfo {
  char *z[3], *spname[3], outf[128], seqf[128], ratef[128], ctlf[128], fix_locusrate;
  int model, ncode, cleandata, seed, npoints, ncatBeta, UseMedianBeta, getSE;
  int ndata, ngene, seqtype, ns, ls, posG[1+1], lgene[1], *pose, npatt, readpattern;
  int *Nij, nGtree;
  /* zzz define new variables*/
  int *state;
  double *fpatt, kappa, alpha, rho, rgene[1], pi[4], piG[1][4];
  double *lnLmax, *locusrate;
  double *pDclass, *tau1beta, *bp0124[18], wwprior[18][32*32][8]; //zzz  igrid could be 16 or 32
}  com;

int LASTROUND, multiplier=100;
double para[8];  /* theta0 theta1 tau0 tau1 qbeta */
int debug = 0;

enum {M0, M1DiscreteBeta, M2SIM3s} MODELS;
enum {G1a, G1b, G1c,G2a, G2b, G2c,G3a, G3b, G3c,G4a, G4b, G4c,G5a, G5b, G5c,G6a, G6b, G6c} GTREES;   //zzz
char *ModelStr[3] = {"M0", "DiscreteBeta", "SIM3s"};
char *GtreeStr[18] = {"G1a", "G1b", "G1c","G2a", "G2b", "G2c","G3a", "G3b", "G3c","G4a", "G4b", "G4c","G5a", "G5b", "G5c","G6a", "G6b", "G6c"};//zzz

#define REALSEQUENCE
#include "treesub.c"

int GetOptions (char *ctlf);
int ReadSeqData(FILE*fout, char seqfile[], char ratefile[], int cleandata);
int Initialize3s(double space[]);
double lfun(double x[], int np);
void p0124Fromb0b1 (double p[5], double b[2]);
double Models0123(FILE *fout, FILE *frub, FILE *frst, double x[], double space[]);
int Simulation(FILE *fout, FILE *frub, double space[]);


FILE *fout, *frub, *frst;

int main (int argc, char* argv[])
{
  char VerStr[32] = "Version 2.1, August 2011";
  int i, j;
  double x[7]={1,1,1,1,1,1,1}, a,space[1000]; //zzz
	
  strcpy(com.ctlf, "3s.ctl");
  if(argc>1) strcpy(com.ctlf, argv[1]);
  starttimer();
  GetOptions(com.ctlf);
	
  if(com.seed<=0) com.seed = abs((2*(int)time(NULL)+1));
  SetSeed(com.seed, 0);
	
  fout=gfopen(com.outf, "w");
  frst=gfopen("rst", "w");
  frub=(FILE*)gfopen("rub","w");
  printf("3s (%s)\n", VerStr);
  fprintf(fout, "3s (%s)\n", VerStr);
	
  ReadSeqData(fout, com.seqf, com.ratef, com.cleandata);
  for (i=0; i<5; i++){
    for (j=0; j<5; j++){
      printf("  %d  ", com.Nij[i*5+j]);
    }
    printf("\n");
  }
  getchar();
	
  /* Debugging
     printf("\n------"); 
     printf("nloci=%4d\n",com.ndata);
     for(i=0;i<com.ndata*5;i++){
     printf("%3d  ",com.Nij[i]);
     if(i%5==4) printf("\n");
     }
  */
	
  a=Models0123(fout, frub, frst, x, space);
  printf("2*Dln=%f\n",a);
  free(com.Nij);
	
  fclose(frst);
  fclose(frub);
  fclose(fout);
  return 0;
}


int GetOptions (char *ctlf)
{
  com.ncode = 4;
  com.npoints = 4;
  com.ncatBeta = 5;
  com.UseMedianBeta = 0;   /* 1 to use the median */
  com.cleandata = 1;
  com.fix_locusrate = 0;
  com.seed = -1;
  com.ndata = 5;
  strcpy(com.outf, "out");
  strcpy(com.seqf, "seqfile.txt");
  com.getSE = 1;
  Small_Diff = 0.5e-6;
  return(0);
}

double Models0123 (FILE *fout, FILE *frub, FILE *frst, double x[], double space[])
{
  int np, i,j, s, K=com.npoints,nG=com.nGtree;
  char timestr[96];
  double *var, lnL, lnL0=0, e=1e-8;
  double xb[7][2];
  int my_a;

  printf("Inside Model0123: begin\n");
  getchar();

  com.pDclass = (double*)malloc((com.ncatBeta*com.ndata+com.ncatBeta)*sizeof(double));
  if(com.pDclass==NULL) error2("oom Models01231");
  com.tau1beta = com.pDclass + com.ncatBeta*com.ndata;
  s = (com.fix_locusrate ? K*K*2 : K*K*5);  /* 2 for b0 & b1; 5 for p0124 */
  com.bp0124[0] = (double*)malloc(18*s*sizeof(double));
  for(i=1;i<18;i++) com.bp0124[i] = com.bp0124[i-1] + s;
  /*for(j=0;j<18;j++){
    com.wwprior[j] = (double**)malloc(K*K*sizeof(double*)); 
    for(i=0;i<K*K;i++){
    com.wwprior[j][i] = malloc(8 * sizeof(double));
    }*/
  //for(i=1;i<nG;i++){
  //com.bp0124[i] = com.bp0124[i-1] + s;
  // com.wwprior[i] = com.wwprior[i-1] + K*K*8;
  //}//zzz
  /*com.bp0124[1] = com.bp0124[0] + s;
    com.bp0124[2] = com.bp0124[1] + s;
    com.wwprior[0] = (double*)malloc(3*K*K*sizeof(double));
    com.wwprior[1] = com.wwprior[0] + K*K;
    com.wwprior[2] = com.wwprior[1] + K*K;*/
  if(com.bp0124[0]==NULL || com.wwprior[0]==NULL) 
    error2("oom Models01232");

  Initialize3s(space);
  printf("\n here is 2.2\n");
	
  // zzz change the com.model=0 to com.model=2 to test the model 2 only.
  for(com.model=0; com.model<1; com.model+=2) {
    if(com.model==M0){   np=6; com.nGtree=18; }
    else if(com.model==M1DiscreteBeta)   { np=5;  com.nGtree=2;  }
    else if(com.model==M2SIM3s)  { np=7;  com.nGtree=18;  }
		
    LASTROUND = 0;
    printf("\n\n*** Model %d (%s) ***\n", com.model, ModelStr[com.model]);
    GetInitials (np, x, xb);
    x[0] = 200.0;
    x[1] = 30.0;
    x[2] = 4.0;
    x[3] = 0.5;
    x[4] = 2.0;
    x[5] = 200.0;
    var = space + np;
		
    lnL = lfun(x, np);

    /* for(x[0]=0.01; x[0] < 100; x[0] = x[0] + 0.01){ */
    /*   lnL = lfun(x, np); */
    /*   FOR(i,np) printf(" %9.6f", x[i]); */
    /*   printf("  %f\n", lnL); */
    /* } */

    /* printf("np = %d \n", np); */
    /* FOR(i,np) printf(" %9.6f", x[i]); FPN(F0); */
    /* FOR(i,np) printf(" %9.5f", xb[i][0]);  FPN(F0); */
    /* FOR(i,np) printf(" %9.5f", xb[i][1]);  FPN(F0); */
    /* printf("\nlnL0 = %12.6f\n", -lnL); */
	
    NFunCall=0;

    printf("ming2: start\n");
    getchar();
    my_a = ming2(frub, &lnL, lfun, NULL, x, xb, space, e, np);
    //getchar();
    printf("ming2 returned %d\n", my_a);
    printf("HDebug@ming2: end\n");

    if(com.model==M0) lnL0 = lnL;
    LASTROUND = 1;
    x[3] *= x[2];     /* xtau1 -> tau1 */
    printf("\nlnL  = %12.6f %+12.6f\nMLEs\n    ", -lnL, lnL0-lnL);
    if(com.model==M0)      printf("theta4    theta5      tau0      tau1     theta1      theta2(x100)\n");
    if(com.model==M1DiscreteBeta)  printf("theta0    theta1      tau0      tau1 (x100) qbeta\n");
    if(com.model==M2SIM3s) printf("theta4    theta5      tau0      tau1     theta1      theta2(x100) M1&2\n"); 
    for(i=0; i<np; i++)    printf(" %9.6f", x[i]);
	
    if(fout && com.getSE) {
      Hessian (np, x, lnL, space, var, lfun, var+np*np);
      matinv(var, np, np, var+np*np);
      matout2(F0, var, np, np, 10, 5);
      printf("\nSEs:\n");
      for(i=0;i<np;i++) printf(" %9.6f",(var[i*np+i]>0.?sqrt(var[i*np+i]):-1));
      printf("\n");
      fprintf(fout, "SEs:\n");
      for(i=0;i<np;i++) fprintf(fout, " %9.6f", (var[i*np+i]>0. ? sqrt(var[i*np+i]) : -1));
      fprintf(fout, "\n\nCorrelation matrix\n");
      for(i=0;i<np;i++) for(j=0;j<i;j++) 
			  var[j*np+i] = var[i*np+j] =
			    (var[i*np+i]>0 && var[j*np+j]>0 ? var[i*np+j]/sqrt(var[i*np+i]*var[j*np+j]) : -9);
      for(i=0;i<np;i++) var[i*np+i]=1;
      matout2(fout, var, np, np, 10, 5);
      fflush(fout);
    }

    if(com.model==M2SIM3s) {
      for(i=0; i<np; i++) fprintf(frst, "%.6f", x[i]);  
      printf("\nThe other local peak is at\n");
      fprintf(fout, "\nThe other local peak is at\n");
      x[4] = x[4]/(8*x[6]);  /*  theta/(8M)  */
      x[5] = x[5]/(8*x[6]);  /*  theta/(8M)  */
      x[6] = 1/(64*x[6]);    /*  1/(64M)  */
      lnL = lfun(x, np);
      for(i=0; i<np; i++) printf(" %9.6f", x[i]);  
      printf(" %12.6f\n", -lnL);
      for(i=0; i<np; i++) fprintf(fout, " %9.6f", x[i]);  
      fprintf(fout, " %12.6f\n", -lnL);
      fprintf(frst, "\t%.3f\t%.6f\t%.6f\n", lnL0-lnL, x[4], x[5]);
    }
    fflush(frst);
	
    printf("\nTime used: %s\n", printtime(timestr));
    if(com.model==0)  x[3] /= x[2];     /* tau1 -> xtau1, as initial for M1 */
  }

  free(com.pDclass);
  free(com.bp0124[0]);
  free(com.state); //zzz
  return(lnL0-lnL);
}

int complexroots (double Q[], double P[], double t)
{
	
  int i,j,k, ii;
  double rr[DIM], ri[DIM], vr[DIM*DIM], vi[DIM*DIM], space[DIM];
  double now;
  complex cU[DIM*DIM], cV[DIM*DIM], cRoot[DIM], cP[DIM*DIM];
  int l;
	
	
  FOR (ii, 1) {
		
    eigen (1, Q, DIM, rr, ri, vr, vi, space) ;
    FOR (i,DIM)   { cRoot[i].re=rr[i]; cRoot[i].im=ri[i]; }
    FOR (i,DIM*DIM) { cU[i].re=cV[i].re=vr[i]; cU[i].im=cV[i].im=vi[i]; }
    cmatinv (cV, DIM, DIM, space);
    cmatby (cU, cV, cP, DIM, DIM, DIM);
    cPMat (P,t, cU, cV, cRoot);	   
  }
  return (0);
}

int cPMat (double P[],double t,complex cU[],complex cV[],complex cRoot[])
{
  /*P(t) = cU * exp{cRoot*t} * cV */
  int i,j,k, status=0;
  complex cUd[DIM*DIM], cP, cY;
  double sum;
	
  FOR (i,DIM) cUd[i*DIM+0]=cU[i*DIM+0];
  for (j=1; j<DIM; j++) {
    cY.re=cRoot[j].re*t; cY.im=cRoot[j].im*t; cY=my_cexp(cY);
    for (i=0; i<DIM; i++)  cUd[i*DIM+j]=cby(cU[i*DIM+j],cY);
  }
  FOR (i,DIM)   {
    for (j=0,sum=0; j<DIM; j++) {
      for (k=0,cP=compl(0,0); k<DIM; k++) {
	cY = cby(cUd[i*DIM+k],cV[k*DIM+j]);
	cP.re+=cY.re;  cP.im+=cY.im;
      }
      P[i*DIM+j]=cP.re;
      sum+=P[i*DIM+j];
      if (P[i*DIM+j]<=0 || fabs(cP.im)>1e-4) status=-1;
    }
    if (fabs(sum-1)>1e-4) status=-1;
  }
  /*if (status)
    { printf ("\nerr cPMat.."); getchar(); matout (F0, P, DIM, DIM); }*/
	
  return (0);
}

int ReadSeqData(FILE*fout, char seqfile[], char ratefile[], int cleandata)
{
  /* Read sequences at each locus and count sites. 
     All sites with ambiguities are deleted right now.  cleandata is ignored.
  */
  FILE *fin = gfopen(seqfile,"r");
  int i,j, h, *n, nstate[8]={0}; //zzz
  double mr, mNij[5]={0};
  int lspname, k; /* 2012-4-11 by ljf add  */
  int *m;  /* 2012-4-11 by ljf add  */
  char spacename1='1', spacename2='2', spacename3='3'; //zzz The seqfile requires the specious named 1,2 and 3
	
  lspname=LSPNAME; /* 2012-4-11 by ljf add  */
	
  printf("\nReading sequence data..  %d loci\n", com.ndata);
  if((com.Nij=(int*)malloc(com.ndata*5*sizeof(int)))==NULL) error2("oom");
  memset(com.Nij, 0, com.ndata*5*sizeof(int));
	
  /* 2012-4-11 by ljf add  */
  if((com.state=(int*)malloc(com.ndata*sizeof(int)))==NULL) error2("oom");
  memset(com.state, 0, com.ndata*sizeof(int));
  /* 2012-4-11 by ljf add  */
	
  if((com.lnLmax=(double*)malloc(com.ndata*(1+com.fix_locusrate)*sizeof(double)))==NULL) 
    error2("oom lnLmax");
  if(com.fix_locusrate) com.locusrate = com.lnLmax + com.ndata;

  n = com.Nij;
  m = com.state;
  for(i=0; i<com.ndata; i++) {
    fprintf(fout, "\n\n*** Locus %d ***\n", i+1);
    ReadSeq (NULL, fin, cleandata);
    PatternWeightJC69like (fout);
		
    for(h=0; h<com.npatt; h++) {
      if(com.z[0][h]==com.z[1][h] && com.z[0][h]==com.z[2][h])
	n[0] = (int)com.fpatt[h];
      else if(com.z[0][h]==com.z[1][h] && com.z[0][h]!=com.z[2][h])
	n[1] = (int)com.fpatt[h];
      else if(com.z[0][h]!=com.z[1][h] && com.z[1][h]==com.z[2][h])
	n[2] = (int)com.fpatt[h];
      else if(com.z[0][h]==com.z[2][h] && com.z[0][h]!=com.z[1][h])
	n[3] = (int)com.fpatt[h];
      else
	n[4] = (int)com.fpatt[h];
    }
    printf("%5d: %5d%5d%5d%5d%5d%9d\n", i+1, n[0],n[1],n[2],n[3],n[4], n[0]+n[1]+n[2]+n[3]+n[4]);
		
    /* 2012-4-11 by ljf add  
       if(i==0) 
       {
       spacename1=(char*)malloc((lspname+1)*sizeof(char));
       for(k=0; k<lspname+1; k++) spacename1[k]=0;
		
       spacename3=(char*)malloc((lspname+1)*sizeof(char));
       for(k=0; k<lspname+1; k++) spacename3[k]=0;
		  
       for(k=1; k<lspname+1; k++) strncpy (spacename1,com.spname[0], k);
			
       for(k=1; k<lspname+1; k++) strncpy (spacename3,com.spname[2], k);
			  
       }
				
       if(*com.spname[0]!=*com.spname[1]&&spacename2==NULL)
       {
       spacename2=(char*)malloc((lspname+1)*sizeof(char));
       for(k=0; k<lspname+1; k++) spacename2[k]=0;
				  
       for(k=1; k<lspname+1; k++) strncpy (spacename2,com.spname[1], k);
       }
       /* 2012-4-11 by ljf add  */
		
    /* 2012-4-11 by zzz add  */
    if(*com.spname[0]==spacename1 && *com.spname[1]==spacename1 && *com.spname[2]==spacename1){
      *m=1; k=111;
    }
    else if(*com.spname[0]==spacename1 && *com.spname[1]==spacename1 && *com.spname[2]==spacename2){
      *m=2; k=112;
    }
    else if(*com.spname[0]==spacename1 && *com.spname[1]==spacename2 && *com.spname[2]==spacename1){
      *m=3; k=121;
    }
    else if(*com.spname[0]==spacename1 && *com.spname[1]==spacename2 && *com.spname[2]==spacename2){
      *m=4; k=122;
    }
    else if(*com.spname[0]==spacename2 && *com.spname[1]==spacename1 && *com.spname[2]==spacename1){
      *m=5; k=211;
    }
    else if(*com.spname[0]==spacename2 && *com.spname[1]==spacename1 && *com.spname[2]==spacename2){
      *m=6; k=212;
    }
    else if(*com.spname[0]==spacename2 && *com.spname[1]==spacename2 && *com.spname[2]==spacename1){
      *m=7; k=221;
    }
    else if(*com.spname[0]==spacename2 && *com.spname[1]==spacename2 && *com.spname[2]==spacename2){
      *m=8; k=222;	  
    }
    else error2("with sequence 3");
    nstate[*m-1]++;
    fprintf(fout, "\n/* The initial state of locus %d is %d */\n", i+1,k);
		
    /* 2012-4-11 by zzz add  */
		
    for(j=0; j<5; j++)
      mNij[j] += (double)n[j]/(com.ndata*com.ls);

    n = n + 5;
    m = m + 1;
  }
  free(com.fpatt);
	
  /* 2012-4-18 by zzz add  */
	
  fprintf(fout, "\n\n/* The counts of the initial states are*/\n");
  fprintf(fout, "111  112  121  122  211  212  221  222\n");
  for(i=0;i<8;i++)
    fprintf(fout, "%5d", nstate[i]);
	
  fprintf(fout,"\n");
  /* 2012-4-18 by zzz add  */
	
  for(i=0; i<com.ns; i++) {
    free(com.spname[i]); 
    free(com.z[i]);
  }
  fclose(fin);
	
  printf("\n\nmean Nij: %8.4f %8.4f %8.4f %8.4f %8.4f\n", mNij[0],mNij[1],mNij[2],mNij[3],mNij[4]);

  if(com.fix_locusrate) {
    if((fin = gfopen(ratefile, "r")) == NULL)
      error2("ratefile open error");
    for(i=0,mr=0; i<com.ndata; i++) {
      if(fscanf(fin, "%lf", &com.locusrate[i]) != 1) 
	error2("rate file..");
      mr = (mr*i + com.locusrate[i])/(i+1.0);
    }
    fclose(fin);
    for(i=0; i<com.ndata; i++)  com.locusrate[i] /= mr;
    printf("\nRelative rates for %d loci, scaled to have mean 1, will be used as constants\n", com.ndata);
    fprintf(fout, "\n\nRelative rates for %d loci, scaled to have mean 1, will be used as constants\n", com.ndata);
    fprintf(fout, "theta's & tau's are defined using the average rate\n");
  }
  return(0);
}

static int locus_save;
static double f0124_locus[5];

void p0124Fromb0b1 (double p[5], double b[2])
{
  /* This calculates p0,p1,p2,p3,p4 for the 5 site patterns for 3 species, 
     given branch lengths b0 and b1.
  */
  double e1, e2, e3;
	
  e1 = exp(-4./3*b[1]);
  e2 = exp(-8./3*(b[0]+b[1]));
  e3 = e1*e2;
  e1 = e1*e1;
  p[0]      = (1 + 3*e1 +  6*e2 +  6*e3)/16;
  p[1]      = (3 + 9*e1 -  6*e2 -  6*e3)/16;
  p[2]=p[3] = (3 - 3*e1 +  6*e2 -  6*e3)/16;
  p[4]      = (6 - 6*e1 - 12*e2 + 12*e3)/16;
}

double lnLb0b1 (double b[], int np)
{
  /* this calculates the lnL for the gene tree for branch lengths b0 and b1.
   */
  int *n = com.Nij + locus_save*5;
  double lnL=0, *f=f0124_locus, p[5];
	
  p0124Fromb0b1 (p, b);
  lnL = f[0]*log(p[0]);
  if(f[1]) lnL += f[1]*log(p[1]);
  if(f[2]) lnL += f[2]*log(p[2]);
  if(n[4]) lnL += f[4]*log(p[4]);
  return(-lnL);
}

int Initialize3s(double space[])
{
  int i, *n=com.Nij, n123max;
  double nt, lnL, *f=f0124_locus, b[2], bb[2][2]={{1e-6,1},{1e-6,1}}, e=1e-7;
  double din, dout;
	
  for(locus_save=0; locus_save<com.ndata; locus_save++,n+=5)  {
    for(i=0,nt=0; i<5; i++)  
      nt += n[i];
		
    /* this max may be too large and is not used. */
    for(i=0,com.lnLmax[locus_save]=0; i<5; i++) {
      if(n[i]) 
	com.lnLmax[locus_save] += n[i] * log(n[i]/(double)nt);
    }
		
    n123max = max2(n[1], n[2]);
    n123max = max2(n123max, n[3]);
    f[0] = n[0]                     / nt;
    f[1] = n123max                  / nt;
    f[2] = (n[1]+n[2]+n[3]-n123max) / nt;
    f[4] = n[4]                     / nt;
    din = f[2]+f[4]; dout=f[1]+f[4]+f[2]/2;
    b[1] = din/2;
    b[0] = (dout-din)/2;
    ming2(NULL, &lnL, lnLb0b1, NULL, b, bb, space, e, 2);
    com.lnLmax[locus_save] = -lnL*nt - 200;
  }
  return(0);
}


double lnpD_locus (int locus, int state);
int GetUVrootIM3s(double U[5*5], double V[5*5], double Root[5]);
int GetPMatIM3s(double Pt[5*5], double t, double U[5*5], double V[5*5], double Root[5]);

/* Quadrature using Gauss-Legendre rule.  The following transforms are used.
   x: (-1, 1) <== y: (a=0, b)  <==  t: (0, inf || taugap)
   y = (b+a)/2 + (b-a)/2*x; 
   t = y/(1-y)
   Three G trees are ordered G1c23, G1b, G1a
*/
int get_ij (int i, int j, int k) /*calculate the location of P_{ij}*/
{
  return ((i-1)*k+j-1);
}

double is2 (int i, double P[], int k)
{
  double g=0;
  int j;
  for (j=9;j<21;j++)
    g += P[get_ij(i,j,k)];
  return(g);
}

double is3 (int i, double P[],int k)
{
  double g=0;
  int j;
  for (j=0;j<9;j++)
    g += P[get_ij(i,j,k)];
  return(g);
}

void debug_print(double *pt, int dim){
  int i;
  for(i=0;i<dim*dim;i++){
    printf("%0.4f ", pt[i]);
    if((i + 1) % dim == 0) printf("\n");
  }
}


double lfun (double x[], int np)
{
  int *n, Gtree, locus, itau1, igrid, i,j, ixw[2], K=com.npoints, ij, s,k;
  double theta4, theta1, theta2, theta5,tau0, tau1, taugap;
  double p, q, mbeta=-1;  /* paras in beta For tau1 under M1DiscreteBeta */
  double *xI=NULL, *wI=NULL;  /* points & weights from GaussLegendre */
  /* yup[] is y upper bound for y0 & y1 */
  double yup[2]={1,1}, y[2], t[2], b[2], r[2], coeff, Li, lnL=0, z, t1, f[8]={0},g,p1,p2,p3,p4;
  double Pt0[DIM*DIM]={0}, Pt1[DIM*DIM]={0}, Ptau1[DIM*DIM];
  double M12, c1,c2, m1,m2, a;
  double Q[DIM*DIM]={0},Q1[DIM*DIM]={0};
  double theta[8]={0}; //zzz
	
  double my_coeff;
	
  xtoy(x, para, np);
  if(!LASTROUND) para[3] = para[3]*para[2];  /* tau1 */
  for(j=0; j<6; j++) para[j] /= multiplier;  /* qbeta under M1DiscreteBeta is not scaled */
  if(com.model==M2SIM3s) {
    M12 = para[6];
  }
  theta4=para[0];
  theta5=para[1];
  tau0=para[2];
  tau1=para[3];
  theta1=para[4];
  theta2=para[5];
  c1=2/theta1; 
  c2=2/theta2; 
  m1=4*M12/theta1; 
  m2=4*M12/theta2;
  if(tau0<tau1) error2("tau0<tau1");
  //printf("theta1=%f, theta2=%f, theta4= %f, theta5=%f tau0=%f, tau1=%f, M=%f\n",theta1,theta2,theta4,theta5,tau0, tau1,M12); //hhh
	
  if(com.model==M2SIM3s) {
    for (i=0;i<DIM*DIM;i++){
      if(mc1[i]) Q[i]= c1 * mc1[i];
      if(mc2[i]) Q[i]+=c2 * mc2[i];
      if(mw1[i]) Q[i]+=m1 * mw1[i];
      if(mw2[i]) Q[i]+=m2 * mw2[i];
    }
  }
		
  for(s=0;s<8;s++){
    if(s<3 || s==4){
      theta[s]=theta1;
    }
    else{ 
      theta[s]=theta2;
    }
  }
		
  GaussLegendreRule(&xI, &wI, com.npoints);
		
  for (Gtree=0; Gtree<com.nGtree;Gtree++){
    if (Gtree==G1a || Gtree==G1b || Gtree==G1c){
      yup[1]=1;
      yup[0]=tau1;
    }
    else if (Gtree==G2a || Gtree==G2b || Gtree==G2c){
      yup[1]=tau1;
      yup[0]=2*(tau0-tau1)/theta5;
    }
    else if (Gtree==G3a || Gtree==G3b || Gtree==G3c){
      yup[1]=tau1;
      yup[0]=-1; /*-1 means the up is infty*/
    }
    else if (Gtree==G4a || Gtree==G4b || Gtree==G4c){
      yup[1]=1;
      yup[0]=2*(tau0-tau1)/theta5;
    }
    else if (Gtree==G5a || Gtree==G5b || Gtree==G5c){
      yup[1]=2*(tau0-tau1)/theta5;
      yup[0]=-1;
    }
    else if (Gtree==G6a || Gtree==G6b || Gtree==G6c){
      yup[1]=-1;
      yup[0]=-1;
    }
    for(i=0;i<2;i++){
      if(yup[i]==-1) yup[i]=1; /* x->x/(1+x)*/
      else yup[i]=yup[i]/(yup[i]+1);
    }
    coeff=yup[1]*yup[0]*0.25; /*(b-a)/2*(b-a)/2 */
    for(igrid=0; igrid<K*K; igrid++) 
      {
	ixw[0] = igrid/K; ixw[1] = igrid%K;
	for(j=0; j<2; j++) {  /* t0 (y0) and t1 (y1) */
	  if(ixw[j]<K/2) { ixw[j] = K/2-1-ixw[j];  r[j]=-1; }
	  else           { ixw[j] = ixw[j]-K/2;    r[j]=1;  }
	  y[j] = yup[j]*(1 + r[j]*xI[ixw[j]])/2;
	  t[j] = y[j]/(1-y[j]);
	}
	my_coeff = coeff*wI[ixw[0]]*wI[ixw[1]]/square((1-y[0])*(1-y[1])); 
				
	if(Gtree==G1a||Gtree==G1b || Gtree==G1c){
	  if(com.model==M2SIM3s){
	    xtoy(Q,Q1,DIM*DIM);
	    complexroots (Q1, Pt1, t[0]*t[1]);//changed
						
	    if(0 && igrid==0){
	      debug_print(Pt1, DIM);
	    }
	    xtoy(Q,Q1,DIM*DIM);
	    complexroots (Q1, Pt0, t[0]-t[0]*t[1]); /*ÔÚpt.cÀï£¬Ö±½ÓÓÃ¿ÉÄÜÓÐÎÊÌâ*/ //change position
	    if(0 && igrid==0){
	      debug_print(Pt0, DIM);
	    }
	    /*??? this only need to be calculated once, now three times are calculated*/
	    if (Gtree==G1c){
	      p1=2/theta1*Pt0[get_ij(11,11,DIM)]+2/theta2*Pt0[get_ij(11,14,DIM)];
	      p2=2/theta1*Pt0[get_ij(20,11,DIM)]+2/theta2*Pt0[get_ij(20,14,DIM)];
	      p3=2/theta1*Pt0[get_ij(17,11,DIM)]+2/theta2*Pt0[get_ij(17,14,DIM)];
	      p4=2/theta1*Pt0[get_ij(14,11,DIM)]+2/theta2*Pt0[get_ij(14,14,DIM)];
	      for(s=1;s<=8;s++)
		f[s-1]=t[0]*(Pt1[get_ij(s,1,DIM)]*2/theta1* p1 + Pt1[get_ij(s,2,DIM)]*2/theta1*p2
			     +Pt1[get_ij(s,7,DIM)]*2/theta2*p3 + Pt1[get_ij(s,8,DIM)]*2/theta2*p4);
	    }
	    else if (Gtree==G1b){
	      p1 = 2/theta1*Pt0[get_ij(10,10,DIM)]+2/theta2*Pt0[get_ij(10,13,DIM)];
	      p2 = 2/theta1*Pt0[get_ij(16,10,DIM)]+2/theta2*Pt0[get_ij(16,13,DIM)];
	      p3 = 2/theta1*Pt0[get_ij(19,10,DIM)]+2/theta2*Pt0[get_ij(19,13,DIM)];
	      p4 = 2/theta1*Pt0[get_ij(13,10,DIM)]+2/theta2*Pt0[get_ij(13,13,DIM)];
	      for(s=1;s<=8;s++)
		f[s-1]=t[0]*(Pt1[get_ij(s,1,DIM)]*2/theta1*p1 + Pt1[get_ij(s,6,DIM)]*2/theta2*p2
			     +Pt1[get_ij(s,3,DIM)]*2/theta1*p3 + Pt1[get_ij(s,8,DIM)]*2/theta2*p4);
	    }
	    else if (Gtree==G1a){
	      p1=2/theta1*Pt0[get_ij(9, 9,DIM)]+2/theta2*Pt0[get_ij(9, 12,DIM)];
	      p2=2/theta1*Pt0[get_ij(15,9,DIM)]+2/theta2*Pt0[get_ij(15,12,DIM)];
	      p3=2/theta1*Pt0[get_ij(18,9,DIM)]+2/theta2*Pt0[get_ij(18,12,DIM)];
	      p4=2/theta1*Pt0[get_ij(12,9,DIM)]+2/theta2*Pt0[get_ij(12,12,DIM)];
	      for(s=1;s<=8;s++)
		f[s-1]=t[0]*(Pt1[get_ij(s,1,DIM)]*2/theta1*p1 + Pt1[get_ij(s,4,DIM)]*2/theta2*p2
			     +Pt1[get_ij(s,5,DIM)]*2/theta1*p3 + Pt1[get_ij(s,8,DIM)]*2/theta2*p4);	
	    }// if gtree==G1a
	  } /*com.model=M2)*/
					
	  else if (com.model==M0){
	    for(s=0;s<8;s++) f[s]=0;
	    f[0]= t[0]*2/theta[0] * exp(-6/theta[0] * t[0]*t[1]) * 2/theta[0] * exp(-2/theta[0] * (t[0]-t[0]*t[1]));
	    f[7]= t[0]*2/theta[7] * exp(-6/theta[7] * t[0]*t[1]) * 2/theta[7] * exp(-2/theta[7] * (t[0]-t[0]*t[1]));
	  }
	  b[1]=t[0]*t[1];	
	  b[0]=t[0]-t[0]*t[1]; /*correct?*/	
	}
	
	
	else if(Gtree==G2a || Gtree==G2b || Gtree==G2c || Gtree== G3a || Gtree ==G3b || Gtree ==G3c){
	  if(com.model==M2SIM3s)
	    {
	      xtoy(Q,Q1,DIM*DIM);
	      complexroots (Q1, Pt0, tau1-t[1]); /*ÔÚpt.cÀï£¬Ö±½ÓÓÃ¿ÉÄÜÓÐÎÊÌâ*/
	      xtoy(Q,Q1,DIM*DIM);
	      complexroots (Q1, Pt1, t[1]);
	      //printf("\n In lfun after end2\n");
	      for(s=1;s<=8;s++){ //problems about s-1
		if (Gtree==G2c || Gtree==G3c)
		  f[s-1]=(Pt1[get_ij(s,1,DIM)] * is2(11,Pt0,DIM) + Pt1[s,2,DIM]*is2(20,Pt0,DIM)) * 2/theta1
		    +      (Pt1[get_ij(s,7,DIM)] * is2(17,Pt0,DIM) + Pt1[s,8,DIM]*is2(14,Pt0,DIM)) * 2/theta2;
		else if (Gtree==G2b || Gtree==G3b)
		  f[s-1]=(Pt1[get_ij(s,1,DIM)] * is2(10,Pt0,DIM) + Pt1[s,3,DIM]*is2(19,Pt0,DIM))*2/theta1
		    +(Pt1[get_ij(s,6,DIM)] * is2(16,Pt0,DIM) + Pt1[s,8,DIM]*is2(13,Pt0,DIM))*2/theta2;
		else if (Gtree==G2a || Gtree==G3a)
		  f[s-1]=(Pt1[get_ij(s,1,DIM)] * is2( 9,Pt0,DIM) + Pt1[s,5,DIM]*is2(18,Pt0,DIM)) *2/theta1
		    +(Pt1[get_ij(s,4,DIM)] * is2(15,Pt0,DIM) + Pt1[s,8,DIM]*is2(12,Pt0,DIM)) *2/theta2;
	      } 
	    }
	  else if (com.model==M0){
	    for(s=0;s<8;s++) f[s]=0;
	    for(s=0;s<8;s++){
	      if(s==0 || s==7)
		f[s]= 2/theta[s] * exp(-6/theta[s] * t[1]) * exp(2/theta[s]*(tau1-t[1]));
	      else if ( ((s==1 || s==6) && (Gtree == G2c || Gtree==G3c)) ||
			((s==2 || s==5) && (Gtree == G2b || Gtree==G3b)) ||
			((s==3 || s==4) && (Gtree == G2a || Gtree==G3a)) )
		f[s]=2/theta[s] * exp(-6/theta[s] * t[1]);
	      //printf("inside lfun: s %d, f[s] %f\n", s, f[s]);

	    }
	  }
					
	  if(Gtree==G2a || Gtree==G2b || Gtree==G2c){
	    g=exp(-t[0]);
	    for(s=1;s<=8;s++)
	      f[s-1]*=g;
	    b[1]= t[1];
	    b[0]= theta5*t[0]/2.0+tau1-t[1];
	  }
	  else if(Gtree==G3a || Gtree==G3b || Gtree==G3c){
	    g= exp(-t[0])*exp(-2*(tau0-tau1)/theta5);
	    for(s=1;s<=8;s++)
	      f[s-1]*=g;
	    b[1]=t[1];
	    b[0]= theta4*t[0]/2.0+tau0-t[1];
	  } /* else if*/
	} // else if(Gtree==G2a || Gtree==G2b || Gtree==G2c || Gtree== G3a || Gtree ==G3b || Gtree ==G3c)
				
	else{
	  //printf("\n In lfun after end3\n");
	  if (Gtree==G4a || Gtree==G4b || Gtree==G4c){
	    g= exp(-3*t[0]*t[1]-t[0]*(1-t[1]));
	    b[1]=theta5/2*t[0]*t[1]+tau1;
	    b[0]=theta5/2*t[0]+tau1-b[1];
	  }//above checked
	  else if (Gtree==G5a || Gtree==G5b || Gtree==G5c){
	    g=exp(-2*t[1]-t[0])* exp(-2.0*(tau0-tau1)/theta5);
	    b[1]=theta5/2*t[1]+tau1;
	    b[0]=theta4/2*t[0]+tau0-b[1];
	  }
	  else if (Gtree==G6a || Gtree==G6b || Gtree==G6c){
	    g= exp(-3*t[1]-t[0])* exp(-6.0*(tau0-tau1)/theta5);
	    b[1]=theta4/2*t[1]+tau0;
	    b[0]=theta4/2*(t[0]+t[1])+tau0-b[1];
	  } /*end else if*/
	  if(com.model==M2SIM3s){
	    xtoy(Q,Q1,DIM*DIM);
	    complexroots (Q1, Ptau1, tau1);	
	    for(s=1;s<=8;s++)
	      f[s-1]= is3(s,Ptau1,DIM)*g; 
	  }
	  else if (com.model==M0){
	    for(s=0;s<8;s++){
	      if(s==0 || s==7)
		f[s] = exp(-6/theta[s]*tau1)*g;
	      else f[s] = exp(-2/theta[s]*tau1)*g;
	    }
	  }
	} /* end else */
	for(s=0;s<8;s++){
	  com.wwprior[Gtree][igrid][s] = f[s] * my_coeff;
	}
	p0124Fromb0b1(com.bp0124[Gtree]+igrid*5, b);
      } //  for(igrid)
  } //for gtree
		
  for(locus=0; locus<com.ndata; locus++) {
    Li = lnpD_locus(locus, com.state[locus]);
    lnL += Li;
    n = com.Nij + locus*5;
    //printf("%d\t%d\t%d\t%d\t%d\t%.6f\n", n[0],n[1],n[2],n[3],n[4], Li);
  }
  return(-lnL);
}




double lnpD_locus (int locus, int state)
{
  int *n = com.Nij + locus*5, K=com.npoints, n123max, igrid, k, Gtree, y, error=0,flag;
  double lnL=0, pD=0, f, p[5], p12, b[2], sump12;
	
  for(Gtree=0; Gtree<com.nGtree; Gtree++) {
    flag=1;
    if(com.model==M0){
      if(state>1 && state<8 ){
	if(Gtree==G1a || Gtree==G1b || Gtree==G1c) flag=0;
	else if((state==2 ||state==7) && (Gtree==G2a || Gtree ==G2b || Gtree==G3a ||Gtree==G3b)) flag=0;
	else if((state==3 ||state==6) && (Gtree==G2a || Gtree ==G2c || Gtree==G3a ||Gtree==G3c)) flag=0;
	else if((state==4 ||state==5) && (Gtree==G2b || Gtree ==G2c || Gtree==G3b ||Gtree==G3c)) flag=0;
      }
    }
    for(igrid=0; flag==1 && igrid<K*K; igrid++) { //??? correctly use flag? while flag=0, the related genetree does not exist
      if(com.fix_locusrate) {  /* p0124 is calculated only if locus rates are fixed. */
	b[0] = com.bp0124[Gtree][igrid*2+0] * com.locusrate[locus];
	b[1] = com.bp0124[Gtree][igrid*2+1] * com.locusrate[locus];
	p0124Fromb0b1 (p, b);
      }
      else {
	for(k=0; k<5; k++)
	  p[k] = com.bp0124[Gtree][igrid*5+k];
      }
      f = -com.lnLmax[locus];
      if(n[0]) f += n[0]*log(p[0]);
      if(n[4]) f += n[4]*log(p[4]);
      if(Gtree==G1c||Gtree==G2c||Gtree==G3c||Gtree==G4c||Gtree==G5c||Gtree==G6c){
	if(n[1]) f+= n[1]*log(p[1]);
	if(n[2]+n[3]) f+=(n[2]+n[3])*log(p[2]);
      }
      else if(Gtree==G1a||Gtree==G2a||Gtree==G3a||Gtree==G4a||Gtree==G5a||Gtree==G6a){
	if(n[2]) f+= n[2]*log(p[1]);
	if(n[3]+n[1]) f+=(n[3]+n[1])*log(p[2]);
      }
      else{
	if(n[3]) f+= n[3]*log(p[1]);
	if(n[1]+n[2]) f+=(n[1]+n[2])*log(p[2]);
      }
      f = (f<-200 ? 0 : exp(f));
      //if (Gtree == 5) {
      //printf("%-3s %-6d f %9.3e %9.3e\n", GtreeStr[Gtree], igrid+1, f, com.wwprior[Gtree][igrid][state]);
      //}
      pD += f*com.wwprior[Gtree][igrid][state];
			
    }  /* for(igrid) over the grid */
    //printf("%-3s %9.3e\n", GtreeStr[Gtree], pD);
  }     /* for(Gtree) */

  if(pD < 1e-300) {
    error = -1;
    printf("at locus %2d, pD = %.6g\n", locus+1, pD);
    lnL += -1e100 + com.lnLmax[locus];
  }
  else
    lnL += log(pD) + com.lnLmax[locus];
  if(error) {
    matout(F0, para, 1, 4+(com.model==M1DiscreteBeta)+(com.model==M2SIM3s));
    error2("floating point problem");
  }
  return(lnL);
}


int GetInitials (int np, double x[], double xb[][2])
{
  int i;
  double thetaU=499, MU=0.15;  /* MU = 0.125 should be fine */
	
  for(i=0; i<np; i++)  { xb[i][0]=0.001;  xb[i][1]=thetaU; }
  xb[3][1] = 0.999;  /* xtau */
	
  if(com.model==M0) {
    //if(x[0]<0 || x[0]>9) {
    for(i=0; i<np; i++) x[i] = 0.9+0.5*rndu(); 
    x[3] = 0.5+0.2*rndu();   /*  xtau  */
    //}
  }
  else if(com.model==M1DiscreteBeta) {
    for(i=0; i<3; i++) {
      x[i] *= 0.9+0.2*rndu();
      x[i] = max2(0.001, x[i]);
    }
    x[3] *= 0.8+0.2*rndu();  /* xtau1 */
    x[4] = 1 + 5*rndu();     /* qbeta */
    xb[4][0] = 0.1; /* q_beta */
    xb[4][1] = 499; /* q_beta */
  }
  else if(com.model==M2SIM3s) {
    for(i=0; i<3; i++) {
      x[i] *= 0.8+0.4*rndu();
      x[i] = max2(0.001, x[i]);
    }
    x[3] *= 0.8+0.2*rndu();    /* xtau1 */
    for(i=4;i<6; i++){
      x[i] = (x[0] + x[1])/2*(0.8+0.4*rndu());  /* theta1 and theta2 */
      xb[i][0] = 0.001;  xb[i][1] = thetaU;     /* theta1 and theta2 */ 
    }
    x[6] = 0.05; // + 0.1*rndu();                 /* M12 */
    xb[6][0] = 0.00001;                       /* M12 */ 
    xb[6][1] = MU;
  }
	
  printf("\nInitials & bounds\n    ");
  if(com.model==M0)     printf("theta4    theta5   tau0  xtau1  theta1   theta2 (x100)\n");
  if(com.model==M1DiscreteBeta) printf("theta0    theta1   tau0 (x100) xtau1     qbeta\n");
  if(com.model==M2SIM3s) printf("theta4    thet5   tau0      tau1   theta1   theta2 (x100)     M1&2 \n"); 
  FOR(i,np) printf(" %9.6f", x[i]); FPN(F0);
  FOR(i,np) printf(" %9.5f", xb[i][0]);  FPN(F0);
  FOR(i,np) printf(" %9.5f", xb[i][1]);  FPN(F0);
  return(0);
}
