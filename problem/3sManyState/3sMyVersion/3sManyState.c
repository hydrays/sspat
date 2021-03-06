/* 3s.c

  Ziheng Yang, 2002, 2009, 2011
  
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
int ReadSiteCounts(char *datafile);
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
	int i;
	double x[7]={1,1,1,1,1,1,1}, a,space[1000]; //zzz
	FILE *fhu; //zzz
	fhu=gfopen("lastline.txt", "a"); //zzz
	
	strcpy(com.ctlf, "3s.ctl");
	if(argc>1) strcpy(com.ctlf, argv[1]);
	starttimer();
	GetOptions (com.ctlf);
	
	if(com.seed<=0) com.seed = abs((2*(int)time(NULL)+1));
	SetSeed(com.seed, 0);
	
	fout=gfopen(com.outf, "w");
	frst=gfopen("rst", "w");
	frub=(FILE*)gfopen("rub","w");
	printf("3s (%s)\n", VerStr);
	fprintf(fout, "3s (%s)\n", VerStr);
	
	
	//zzz
	//test loop
	//int index;
	
	
	
	//ReadSiteCounts(datafile); 
#if(1)
	printf("here1\n"); //hhh
	ReadSeqData(fout, com.seqf, com.ratef, com.cleandata);
	
	/* Debugging
	printf("\n------"); 
	printf("nloci=%4d\n",com.ndata);
	for(i=0;i<com.ndata*5;i++){
	printf("%3d  ",com.Nij[i]);
	if(i%5==4) printf("\n");
	}
	*/
	
	noisy = 3;
	for(i=0; i<7; i++) x[i] = 0.5+0.5*rndu(); //zzz
	printf("here2\n"); //hhh
	a=Models0123(fout, frub, frst, x, space);
	printf("here3\n"); //hhh
	printf("2*Dln=%f\n",a);
	fprintf(fhu, "2*DlnL = %f\n", a); //
	free(com.Nij);
	free(com.state); //zzz
#else
	Simulation(fout, frub, space);
#endif
	
	fclose(frst);
	fclose(frub);
	fclose(fout);
	
	fclose(fhu); //zzz
	
	return 0;
}



int GetOptions (char *ctlf)
{
	int iopt,i, nopt=9, lline=4096;
	char line[4096],*pline, opt[32], *comment="*#", *seqerrstr="0EF";
	char *optstr[] = {"seed", "nloci", "outfile", "seqfile", "ratefile", 
		"cleandata", "npoints", "getSE", "Small_Diff"};
	double t=1;
	FILE  *fctl=gfopen (ctlf, "r");
	
	com.ncode = 4;
	com.npoints = 16;
	com.ncatBeta = 5;
	com.UseMedianBeta = 0;   /* 1 to use the median */
	com.cleandata = 1;
	com.fix_locusrate = 0;
	if (fctl) {
		if (noisy) printf ("\nReading options from %s..\n", ctlf);
		for (;;) {
			if(fgets(line, lline, fctl) == NULL) break;
			if(line[0]=='/' && line[1]=='/') 
				break;
			for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
				if (isalnum(line[i]))  { t=1; break; }
				else if (strchr(comment,line[i])) break;
				if (t==0) continue;
				sscanf (line, "%s%*s%lf", opt, &t);
				if ((pline=strstr(line, "="))==NULL)
					continue;
				
				for (iopt=0; iopt<nopt; iopt++) {
					if (strncmp(opt, optstr[iopt], 8)==0)  {
						if (noisy>=9)
							printf ("\n%3d %15s | %-20s %6.2f", iopt+1,optstr[iopt],opt,t);
						switch (iopt) {
						case ( 0): com.seed=(int)t;                    break;
						case ( 1): com.ndata=(int)t;                   
							break;
						case ( 2): sscanf(pline+1, "%s", com.outf);    break;
						case ( 3): sscanf(pline+1, "%s", com.seqf);    break;
						case ( 4): sscanf(pline+1, "%s", com.ratef);   com.fix_locusrate = 1; break;
						case ( 5):
							com.cleandata=(int)t;
							if(com.cleandata!=1) error2("use cleandata = 1");
							break;
						case ( 6): 
							sscanf(pline+1, "%d%d%d", &com.npoints, &com.ncatBeta, &com.UseMedianBeta); 
							break;
						case ( 7): com.getSE=(int)t;                   break;
						case ( 8): Small_Diff=t;                       break;
						}
						break;
					}
				}
				if (iopt==nopt)
				{ printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
		}
		fclose(fctl);
	}
	else
		if (noisy) error2("\nno ctl file..");
		
		return(0);
}


#if 0  /* this is from the old verion, now not in use. */

int ReadSiteCounts (char *datafile)
{
	FILE *fin;
	int i;
	
	fin=(FILE*)gfopen(datafile,"r");
	fscanf(fin, "%d", &com.ndata);
	if((com.Nij=(int*)malloc(com.ndata*5*sizeof(int)))==NULL) error2("oom");
	for(i=0; i<com.ndata*5; i++)
		fscanf(fin, "%d", &com.Nij[i]);  
	
	fclose(fin);
	printf("%d loci\n", com.ndata);
	return(0);
}

#endif

int ReadSeqData(FILE*fout, char seqfile[], char ratefile[], int cleandata)
{
/* Read sequences at each locus and count sites. 
All sites with ambiguities are deleted right now.  cleandata is ignored.
	*/
	FILE *fin = gfopen(seqfile,"r");
	int i,j, h, *n,nstate[8]={0}; //zzz
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
	
	for(i=0,n=com.Nij, m=com.state; i<com.ndata; i++,n+=5,m++) {  /* 2012-4-11 by ljf add m=Lij and m++  */
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
		/* printf("%5d: %5d%5d%5d%5d%5d\n", i+1, n[0],n[1],n[2],n[3],n[4]); */
		printf("%5d: %5d%5d%5d%5d%5d%9d\r", i+1, n[0],n[1],n[2],n[3],n[4], n[0]+n[1]+n[2]+n[3]+n[4]);
		
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
		
		if(noisy>=3) printf("\n");
		
		for(j=0; j<5; j++)
			mNij[j] += (double)n[j]/(com.ndata*com.ls);
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
	double my_coeff; //zzz
	double Pt[5*5], U[5*5], V[5*5], Root[5];
	double Pt0[DIM*DIM]={0}, Pt1[DIM*DIM]={1,2,3,4,5}, Ptau1[DIM*DIM];
	double M12, c1,c2, m1,m2, a;
	double Q[DIM*DIM]={0},Q1[DIM*DIM]={0};
	double theta[8]={0}; //zzz
	
	//printf("\n HDebug: Now I'm in lfun, haha.\n"); //hhh
	
	//for(k=0;Pt1[k]&&k<=DIM*DIM;k++) printf("Pt1[%d]=%f ",k,Pt1[k]);
	
	xtoy(x, para, np);
	if(!LASTROUND) para[3] = para[3]*para[2];  /* tau1 */
	for(j=0; j<6; j++) para[j] /= multiplier;  /* qbeta under M1DiscreteBeta is not scaled */
	if(com.model==M2SIM3s) {
		M12 = para[6];
	}
	theta4=para[0], theta5=para[1], tau0=para[2], tau1=para[3], theta1=para[4],theta2=para[5]; //zzz is this ok? 5 paras in m0, but read 6 paras in, use 5 of them
	c1=2/theta1; c2=2/theta2; m1=4*M12/theta1; m2=4*M12/theta2; //zzz
	if(tau0<tau1) error2("tau0<tau1");
	//printf("theta1=%f, theta2=%f, theta4= %f, theta5=%f tau0=%f, tau1=%f, M=%f\n",theta1,theta2,theta4,theta5,tau0, tau1,M12); //hhh
	
	if(com.model==M1DiscreteBeta) {
		mbeta = para[3]/para[2];
		q = para[4] = x[4];
		p = mbeta/(1-mbeta)*q;
		
		DiscreteBeta(Pt, com.tau1beta, p, q, com.ncatBeta, com.UseMedianBeta);
		for(j=0; j<com.ncatBeta; j++)
			com.tau1beta[j] *= tau0;
	}
	
	//printf("\n here is 2.5.2\n"); //hhh
	for(itau1=0; itau1<(com.model!=M1DiscreteBeta ? 1 : com.ncatBeta); itau1++) { /* for M1DiscreteBeta */
		if(com.model==M0){ //zzz
			for (i=0;i<DIM*DIM;i++){
				if(mc1[i]) Q[i]= c1 * mc1[i];
				if(mc2[i]) Q[i]+=c2 * mc2[i];
			}
		}
		//printf("com.model=%d\n",com.model);//zzz
		/*for(i=0,j=0;i<DIM*DIM;i++){
		printf("%6.2f  ",Q[i]);
		if(Q[i]) j++;
		if(i%DIM==DIM-1) printf("\n");} *///hhh
		//printf("the number of none zero=%d\n",j); //hhh
		
		if(com.model==M1DiscreteBeta){
			tau1 = para[3] = com.tau1beta[itau1];
			//taugap = 2*(tau0 - tau1)/theta1;
		}
		
		if(com.model==M2SIM3s) {
			//a=sqrt(c*c + 16*m*m);
			/*PG1a = (c + 4*m + a)*exp(-(c + 4*m - a)*tau1/2)
			- (c + 4*m - a)*exp(-(c + 4*m + a)*tau1/2);
			PG1a = 1 - PG1a/(2*a);
			
			  GetUVrootIM3s(U, V, Root);
			  GetPMatIM3s(Pt, tau1, U, V, Root);
			  PG1a = Pt[1*5+3] + Pt[1*5+4];
			PSG = (1 - PG1a)*exp(-taugap)*2/3.0; /* this seems no need for sim3smanystates*/
			
			//zzz: compute Q  
			for (i=0;i<DIM*DIM;i++){
				if(mc1[i]) Q[i]= c1 * mc1[i];
				if(mc2[i]) Q[i]+=c2 * mc2[i];
				if(mw1[i]) Q[i]+=m1 * mw1[i];
				if(mw2[i]) Q[i]+=m2 * mw2[i];
			}
			//for(i=0;i<DIM*DIM;i++){printf("4.3f",Q[i]); if(i%DIM==DIM-1) printf("\n");}
		}
		
		for(s=0;s<8;s++){
			if(s<3 || s==4){
				theta[s]=theta1;
				//printf("s=%d, theta=%f  \n",s, theta[s]); //hhh
			}
			else{ 
				theta[s]=theta2;
				//printf("s=%d, theta=%f  \n",s, theta[s]); //hhh			
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
			//printf("coeff=%f",coeff);
			//printf("\nhere is 2.5.2.1\n");
			for(igrid=0; igrid<K*K; igrid++) 
			{
				//printf("\n I'm at aaa\n"); //hhh
				ixw[0] = igrid/K; ixw[1] = igrid%K;
				for(j=0; j<2; j++) {  /* t0 (y0) and t1 (y1) */
					if(ixw[j]<K/2) { ixw[j] = K/2-1-ixw[j];  r[j]=-1; }
					else           { ixw[j] = ixw[j]-K/2;    r[j]=1;  }
					y[j] = yup[j]*(1 + r[j]*xI[ixw[j]])/2;
					t[j] = y[j]/(1-y[j]);
				}
				my_coeff = coeff*wI[ixw[0]]*wI[ixw[1]]/square((1-y[0])*(1-y[1])); 
				
				//printf("\n I'm at bbb\n"); //hhh
				
				if(Gtree==G1a||Gtree==G1b || Gtree==G1c)
				{
					//printf("\n I'm at ccc\n"); //hhh
					
					if(com.model==M2SIM3s){
						//printf("\n I'm at ccc: case 1\n");
						
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
								f[s-1]=t[0] * (Pt1[get_ij(s,1,DIM)]*2/theta1*p1 + Pt1[get_ij(s,6,DIM)]*2/theta2*p2
								+Pt1[get_ij(s,3,DIM)]*2/theta1*p3 + Pt1[get_ij(s,8,DIM)]*2/theta2*p4);
						}
						else if (Gtree==G1a){
							p1=2/theta1*Pt0[get_ij(9, 9,DIM)]+2/theta2*Pt0[get_ij(9, 12,DIM)];
							p2=2/theta1*Pt0[get_ij(15,9,DIM)]+2/theta2*Pt0[get_ij(15,12,DIM)];
							p3=2/theta1*Pt0[get_ij(18,9,DIM)]+2/theta2*Pt0[get_ij(18,12,DIM)];
							p4=2/theta1*Pt0[get_ij(12,9,DIM)]+2/theta2*Pt0[get_ij(12,12,DIM)];
							for(s=1;s<=8;s++)
								f[s-1]= t[0] * (Pt1[get_ij(s,1,DIM)]*2/theta1*p1 + Pt1[get_ij(s,4,DIM)]*2/theta2*p2
								+Pt1[get_ij(s,5,DIM)]*2/theta1*p3 + Pt1[get_ij(s,8,DIM)]*2/theta2*p4);	
							//hhh
							//if(igrid==0) {printf("Pt1=\n"); for(k=0;k<DIM*DIM; k++){printf("%f   ",Pt1[i]); if (k%DIM==DIM-1) printf("\n");}}
							//printf("p1=%f, p2=%f, p3=%f, p4=%f\n",p1,p2,p3,p4);//hhh
							//printf(" pt11=%f pt12=%f pt13=%f, pt14=%f\n",Pt1[get_ij(s,1,DIM)], Pt1[get_ij(s,4,DIM)],Pt1[get_ij(s,5,DIM)],Pt1[get_ij(s,5,DIM)]);
						}// if gtree==G1a
					} /*com.model=M2)*/
					
					else if (com.model==M0){
						//printf("\n I'm at ccc: case 2\n");
						for(s=0;s<8;s++) f[s]=0;
						//printf("\n I'm at ccc: case 2 step 2\n");
						f[0]= t[0]* 2/theta[0] * exp(-6/theta[0] * t[0] * t[1]) * 2/theta[0] * exp(-2/theta[0] * (t[0]-t[0]*t[1]));
						//printf("\n I'm at ccc: case 2 step 3\n");
						f[7]= t[0]* 2/theta[7] * exp(-6/theta[7] * t[0] * t[1]) * 2/theta[7] * exp(-2/theta[7] * (t[0]-t[0]*t[1]));
						//printf("\n I'm at ccc: case 2 step 4\n");		
					}
					
					//printf("\n I'm at ddd\n");
					
					b[1]=t[0]*t[1];	
					b[0]=t[0]-t[0]*t[1]; /*correct?*/
					//printf("test=>%d %d %f %f %f %f \n",K, igrid,yup[0],yup[1],y[0],y[1]);
					//getchar();
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
								f[s]= 2/theta[s] * exp(-6/theta[s] * t[1]) * exp(-2/theta[s]*(tau1-t[1]));
							else if ( ((s==1 || s==6) && (Gtree == G2c || Gtree==G3c)) ||
								((s==2 || s==5) && (Gtree == G2b || Gtree==G3b)) ||
								((s==3 || s==4) && (Gtree == G2a || Gtree==G3a)) )
								f[s]=2/theta[s] * exp(-6/theta[s] * t[1]);
						}
					}
					
					if(Gtree==G2a || Gtree==G2b || Gtree==G2c){
						g=exp(-t[0]);
						for(s=1;s<=8;s++)
							f[s-1]*=g;
						b[1]=t[1];
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
				//printf("In lfun after igrid=0 if 3");
				
				for(s=0;s<8;s++)
				{
					com.wwprior[Gtree][igrid][s] = f[s] * my_coeff;
					//printf("wwprior[%d][%d][%d]=%f  ",Gtree, igrid,s,com.wwprior[Gtree][igrid][s]);
				}
				//printf("\n"); //hhh
				if(!com.fix_locusrate)
					p0124Fromb0b1(com.bp0124[Gtree]+igrid*5, b);
				else {
					com.bp0124[Gtree][igrid*2+0] = b[0];
					com.bp0124[Gtree][igrid*2+1] = b[1];
				}
				// printf("%-3s %5d y%9.5f%9.5f t%9.5f%9.5f b%9.5f%9.5f\n",
				//GtreeStr[Gtree], igrid+1, y[0],y[1], t[0],t[1], b[0],b[1]); 
				//printf("In lfun after igrid loop !!!!\n"); zzz
			} //  for(igrid)
        } //for gtree

		for(i=0;i<20;i++)
		printf("  %f \n",com.wwprior[5][i][5]);
		getchar();	
		
		for(locus=0; locus<com.ndata; locus++) {
			Li = lnpD_locus(locus,com.state[locus]);
			if(com.model==M1DiscreteBeta)
				com.pDclass[itau1*com.ndata+locus] = Li;
			else
				lnL += Li;
			n = com.Nij + locus*5;
			if(debug) printf("%d\t%d\t%d\t%d\t%d\t%.6f\n", n[0],n[1],n[2],n[3],n[4], Li);
		}
	}  /* for(itau1) */


	if(com.model==M1DiscreteBeta) {
		for(locus=0,lnL=0; locus<com.ndata; locus++) {
			for(itau1=0,z=-1e300; itau1<com.ncatBeta; itau1++)
				z = max2(z, com.pDclass[itau1*com.ndata+locus]);
			for(itau1=0,Li=0; itau1<com.ncatBeta; itau1++)
				Li += 1.0/com.ncatBeta * exp(com.pDclass[itau1*com.ndata+locus] - z);
			lnL += z + log(Li);
		}
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
			printf("%-3s %-6d f %9.3e\n", GtreeStr[Gtree], igrid+1, f);
			//getchar();
			pD += f*com.wwprior[Gtree][igrid][state];
			
		}  /* for(igrid) over the grid */
		printf("gtree %d  %d  %d  pd %e\n", Gtree, flag, state, pD);
		getchar();
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
		if(x[0]<0 || x[0]>9) {
			for(i=0; i<np; i++) x[i] = 0.9+0.5*rndu(); 
			x[3] = 0.5+0.2*rndu();   /*  xtau  */
		}
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
		x[6] = 0.05 + 0.1*rndu();                 /* M12 */
		xb[6][0] = 0.00001;                       /* M12 */ 
		xb[6][1] = MU;
	}
	
	if(noisy) {
		printf("\nInitials & bounds\n    ");
		if(com.model==M0)     printf("theta4    theta5   tau0  xtau1  theta1   theta2 (x100)\n");
		if(com.model==M1DiscreteBeta) printf("theta0    theta1   tau0 (x100) xtau1     qbeta\n");
		if(com.model==M2SIM3s) printf("theta4    thet5   tau0      tau1   theta1   theta2 (x100)     M1&2 \n"); 
		FOR(i,np) printf(" %9.6f", x[i]); FPN(F0);
		FOR(i,np) printf(" %9.5f", xb[i][0]);  FPN(F0);
		FOR(i,np) printf(" %9.5f", xb[i][1]);  FPN(F0);
	}
	return(0);
}

double Models0123 (FILE *fout, FILE *frub, FILE *frst, double x[], double space[])
{
	int np, i,j, s, noisy0=noisy, K=com.npoints,nG=com.nGtree;
	char timestr[96];
	double *var, lnL, lnL0=0, e=1e-8;
	double xb[7][2];
	
	printf("\n here is 2.1\n"); //hhh
	com.pDclass = (double*)malloc((com.ncatBeta*com.ndata+com.ncatBeta)*sizeof(double));
	if(com.pDclass==NULL) error2("oom Models01231");
	com.tau1beta = com.pDclass + com.ncatBeta*com.ndata;
	s = (com.fix_locusrate ? K*K*2 : K*K*5);  /* 2 for b0 & b1; 5 for p0124 */
	com.bp0124[0] = (double*)malloc(18*s*sizeof(double));
	for(i=1;i<18;i++) com.bp0124[i] = com.bp0124[i - 1] + s;
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
	printf("\n here is 2.2\n");
	noisy = 0;
	Initialize3s(space);
	printf("\n here is 2.3\n");
	noisy = noisy0;
	
	// zzz change the com.model=0 to com.model=2 to test the model 2 only.
	for(com.model=0; com.model<3; com.model+=2) {
		if(com.model==M0){   np=6; com.nGtree=18; }
		else if(com.model==M1DiscreteBeta)   { np=5;  com.nGtree=2;  }
		else if(com.model==M2SIM3s)  { np=7;  com.nGtree=18;  }
		
		LASTROUND = 0;
		if(noisy) printf("\n\n*** Model %d (%s) ***\n", com.model, ModelStr[com.model]);
		if(fout) {
			fprintf(fout, "\n\n*** Model %d (%s) ***\n", com.model, ModelStr[com.model]);
			fprintf(frub, "\n\n*** Model %d (%s) ***\n", com.model, ModelStr[com.model]);
		}
		
		GetInitials (np, x, xb);
		x[0] = 2.0;
		x[1] =3.0;
		x[2] =4.0;
		x[3] =0.5;
		x[4] =2.0;
		x[5] =2.0;
		lnL =lfun(x,np);
		printf("%f/n",lnL);
	    getchar();
		
		printf("\n here is 2.4\n");
		
		/*
		printf("input initials? ");
		for(i=0; i<np; i++)
		scanf("%lf", &x[i]);
		*/
		var = space + np;
		
#if(0)
		printf("\nTesting lfun at fixed parameter values using ChenLi data of 53 loci\n");  
		x[0]=0.30620; x[1]=0.09868; x[2]=0.62814; x[3]=0.82706; x[4]=1.1; x[5]=2.2;
		/* x[0]=0.35895; x[1]=0.42946; x[2]=0.66026; x[3]=0.65449; x[4]=2; */
		printf("  0.30620  0.09868  0.62814  0.82706  lnL = %12.6f (K=inf)\n", -3099.411263);
		
		for(com.npoints=32; com.npoints<=64; com.npoints*=2) {
			for(x[4]=0.5; x[4]<1000; x[4] *= 2) {
				for(i=0; i<np; i++) printf(" %8.5f", x[i]);  
				lnL = lfun(x, np);
				printf("  lnL = %12.6f (%2d) %s\n", -lnL, com.npoints, printtime(timestr));
				if(com.model==M0 || com.model==M2SIM3s) break;
			}
			break;
		}
		exit(0);
#elif(0)
		printf("\nContour for the hominoid data (model %s)\n", ModelStr[com.model]);
		LASTROUND = 1;
		{
			int ii, jj, nii=22, njj=22;
			double x0[]={0.361874, 0.369049, 0.658947, 0.459227, 2.601699, 12.508765};
			double theta12Set[] = {0.5, 0.7, 0.8, 1, 1.5, 2, 2.5, 2.60, 2.7, 2.8, 
				3, 3.2, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10};
			double Mset[]       = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1, 2, 
				4, 5, 7, 8, 10, 12, 15, 16, 18, 20, 22, 25};
			
			xtoy(x0, x, np);
			for(ii=0; ii<nii; ii++) 
				x[4] = theta12Set[ii];
			for(jj=0; jj<njj; jj++) {
				x[5] = Mset[jj];
				lnL = lfun(x, np);
				printf("%.6f\t%.6f\t%.4f\n", x[4], x[5], -lnL);
			}
		}
		exit(0);
	}
#endif
	
	noisy=3;
	printf("here is 2.5\n"); //hhh
	lnL = lfun(x, np);
	printf("here is 2.6\n"); //hhh
	if(noisy) printf("\nlnL0 = %12.6f\n", -lnL);
	
	NFunCall=0;
	ming2(frub, &lnL, lfun, NULL, x, xb, space, e, np);
	
	if(com.model==M0) lnL0 = lnL;
	LASTROUND = 1;
	x[3] *= x[2];     /* xtau1 -> tau1 */
	if(noisy) {
		printf("\nlnL  = %12.6f %+12.6f\nMLEs\n    ", -lnL, lnL0-lnL);
		if(com.model==M0)      printf("theta4    theta5      tau0      tau1     theta1      theta2(x100)\n");
		if(com.model==M1DiscreteBeta)  printf("theta0    theta1      tau0      tau1 (x100) qbeta\n");
		if(com.model==M2SIM3s) printf("theta4    theta5      tau0      tau1     theta1      theta2(x100) M1&2\n"); 
		for(i=0; i<np; i++)    printf(" %9.6f", x[i]);
	}
	if(fout) {
		fprintf (fout, "\nlnL  = %12.6f %+12.6f\nMLEs\n    ", -lnL, lnL0-lnL);
		if(com.model==M0)      fprintf(fout, "theta4    theta5      tau0      tau1     theta1      theta2(x100)\n");
		if(com.model==M1DiscreteBeta)  fprintf(fout, "theta0    theta1      tau0      tau1 (x100) qbeta\n");
		if(com.model==M2SIM3s) fprintf(fout, "theta0    theta1      tau0      tau1  theta1 (x100) theta2 (x100) M1&2\n");
		for(i=0; i<np; i++)    fprintf(fout, " %9.6f", x[i]);
		fprintf(fout, "\n");
		fflush(fout);
	}
	
	if(fout && com.getSE) {
		Hessian (np, x, lnL, space, var, lfun, var+np*np);
		matinv(var, np, np, var+np*np);
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
	if(com.model==M1DiscreteBeta) {
		fprintf(fout, "\ntau1 from the discrete beta (ncat = %d), using %s\ntau1: ", com.ncatBeta, (com.UseMedianBeta?"median":"mean"));
		for(i=0; i<com.ncatBeta; i++) fprintf(fout, " %9.6f", com.tau1beta[i]);
		fprintf(fout, "\nfreq: ");
		for(i=0; i<com.ncatBeta; i++) fprintf(fout, " %9.6f", 1.0/com.ncatBeta);
		fprintf(fout, "\n");
		for(i=0; i<np; i++) fprintf(frst, " %9.6f", x[i]);  
		fprintf(frst, " %8.3f\n", lnL0-lnL);
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
}  /* for(com.model) */

free(com.pDclass);
free(com.bp0124[0]);
free(com.state); //zzz
return(lnL0-lnL);
}

int Simulation (FILE *fout, FILE *frub, double space[])
{
	char timestr[96];
	double theta0, theta1, tau0, tau1, mbeta, pbeta, qbeta, tau1beta[5];
	double t[2], b[2], p[5], y, pG0, *dlnL;
	double md12, md13, md23, d12, d13, d23, Ed12, Ed13, EtMRCA, mNij[5], x[5];
#if(0)
	double x0[5] = {0.04, 0.06, 0.06, 0.04, 2.0};
#elif(0)
	double x0[5] = {0.0035, 0.0060, 0.0066, 0.0041, 0};
#elif(1)
	double x0[5] = {0.005, 0.005, 0.006, 0.004, 0}; /* hominoid (BY08) */
#elif(0)
	double x0[5] = {0.01, 0.01, 0.02, 0.01, 0};    /* mangroves (Zhou et al. 2007) */
#endif
	int model=1, nii=1, nloci[]={1000, 100, 10}, ii, i,j, nr=200, ir, locus;
	
	/*
	printf("input model theta0 theta1 tau0 tau1 q? ");
	scanf("%d%lf%lf%lf%lf%lf", &model, &x0[0], &x0[1], &x0[2], &x0[3], &x0[4]);
	*/
	
	com.ls = 500;
	noisy = 2;
	if(com.fix_locusrate) error2("fix_locusrate in Simulation()?");
	if((dlnL=(double*)malloc(nr*sizeof(double))) == NULL) error2("oom dlnL");
	memset(dlnL, 0, nr*sizeof(double));
	for(ii=0; ii<nii; ii++) {
		com.ndata = nloci[ii];
		printf("\n\nnloci = %d  length = %d\nParameters", com.ndata, com.ls);
		matout(F0, x0, 1, 4+model);
		fprintf(fout, "\n\nnloci = %d  length = %d\nParameters", com.ndata, com.ls);
		matout(fout, x0, 1, 4+model);
		fprintf(frst, "\n\nnloci = %d  length = %d  K = %d\nParameters", com.ndata, com.ls, com.npoints);
		matout(frst, x0, 1, 4+model);
		if((com.lnLmax=(double*)malloc(com.ndata*1*sizeof(double)))==NULL) error2("oom lnLmax");
		if((com.Nij=(int*)malloc(com.ndata*5*sizeof(int)))==NULL) error2("oom");
		theta0=x0[0]; theta1=x0[1]; tau0=x0[2]; tau1=x0[3]; qbeta=x0[4];
		pG0 = 1 - exp(-2*(tau0-tau1)/theta1);
		Ed12 = tau1 + theta1/2 + (1-pG0)*(theta0-theta1)/2;
		Ed13 = tau0 + theta0/2;
		EtMRCA = tau0 + theta0/2 * (1 + (1-pG0)*1./3);
		printf("pG0 = %9.6f  E{tMRCA} = %9.6f\n", pG0, EtMRCA);
		
		md12=0; md13=0; md23=0;      
		for(i=0; i<5; i++) mNij[i]=0;
		for(ir=0; ir<nr; ir++) {
			memset(com.Nij, 0, com.ndata*5*sizeof(int));
			for(locus=0; locus<com.ndata; locus++) {
				if(model==M1DiscreteBeta) {
					mbeta = x0[3]/x0[2];  pbeta = mbeta/(1-mbeta)*qbeta;
					if(0)    /* continuous beta */ 
						tau1 = tau0*rndbeta(pbeta, qbeta);
					else {   /* discrete beta */
						DiscreteBeta(p, tau1beta, pbeta, qbeta, com.ncatBeta, 0);
						tau1 = tau0*tau1beta[(int)(com.ncatBeta*rndu())];
					}               
				}
				t[0] = rndexp(1.0);
				t[1] = rndexp(1.0);
				if(t[1]<2*(tau0-tau1)/theta1) {   /* G0 */
					b[0] = tau0 - tau1 - theta1*t[1]/2 + theta0*t[0]/2;
					b[1] = tau1 + theta1*t[1]/2;
					p0124Fromb0b1 (p, b);
				}
				else {                            /* G123 */
					t[1] = rndexp(1.0/3);
					b[0] = t[0]*theta0/2;
					b[1] = tau0 + t[1]*theta0/2;
					p0124Fromb0b1 (p, b);
					y = rndu();
					if(y<1.0/3)         { y=p[1]; p[1]=p[2]; p[2]=y; }
					else if (y<2.0/3)   { y=p[1]; p[1]=p[3]; p[3]=y; }
				}
				
				for(j=0; j<4; j++)  p[j+1] += p[j];
				for(i=0; i<com.ls; i++) {
					for (j=0,y=rndu(); j<5-1; j++) 
						if (y<p[j]) break;
						com.Nij[locus*5+j] ++;
				}
				
				d12 = (com.Nij[locus*5+2]+com.Nij[locus*5+3]+com.Nij[locus*5+4])/(double)com.ls;
				d13 = (com.Nij[locus*5+1]+com.Nij[locus*5+2]+com.Nij[locus*5+4])/(double)com.ls;
				d23 = (com.Nij[locus*5+1]+com.Nij[locus*5+3]+com.Nij[locus*5+4])/(double)com.ls;
				d12 = -3/4.*log(1 - 4./3*d12);
				d13 = -3/4.*log(1 - 4./3*d13);
				d23 = -3/4.*log(1 - 4./3*d23);
				md12 += d12/(2.0*nr*com.ndata);
				md13 += d13/(2.0*nr*com.ndata);
				md23 += d23/(2.0*nr*com.ndata);
				for(i=0; i<5; i++)
					mNij[i] += (double)com.Nij[locus*5+i]/(com.ndata*com.ls);
			}
			
			printf("\nReplicate %3d\n", ir+1);
			fprintf(fout, "\nReplicate %3d\n", ir+1);
			fprintf(frub, "\nReplicate %3d\n", ir+1);
			
			printf("\nmean Nij: %8.4f %8.4f %8.4f %8.4f %8.4f\n", mNij[0],mNij[1],mNij[2],mNij[3],mNij[4]);
			for(i=0; i<4; i++)
				x[i] = x0[i]*multiplier*(0.8+0.4*rndu());
			x[4] = x0[4]*(0.8+0.4*rndu());
			x[3] = x0[3]/x0[2];
			dlnL[ir] = Models0123(fout, frub, frst, x, space);
			
			printf("%3d/%3d %9.3f %ss\n", ir+1, nr, dlnL[ir], printtime(timestr));
			for(i=0; i<5; i++) mNij[i] = 0;
		}
		
		printf("\nd12 %9.5f = %9.5f d13 d23: %9.5f = %9.5f = %9.5f\n", Ed12, md12, Ed13, md13, md23);
		
		fprintf(fout, "\n\nList of DlnL\n");
		for(i=0; i<nr; i++) fprintf(fout, "%9.5f\n", dlnL[i]);
		free(com.Nij);  free(com.lnLmax);
	}
	printf("\nTime used: %s\n", printtime(timestr));
	return 0;
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
		cY.re=cRoot[j].re*t; cY.im=cRoot[j].im*t; cY=cexp(cY);
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
