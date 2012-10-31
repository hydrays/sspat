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
#define NS            3
#define NGENE         1  /* required by ReadSeq, but not really in use */
#define LSPNAME       50
#define NCODE         4

extern double Small_Diff;
extern int noisy, NFunCall;

struct CommonInfo {
   char *z[3], *spname[3], outf[128], seqf[128], ratef[128], ctlf[128], fix_locusrate;
   int model, ncode, cleandata, seed, npoints, ncatBeta, UseMedianBeta, getSE;
   int ndata, ngene, seqtype, ns, ls, posG[1+1], lgene[1], *pose, npatt, readpattern;
   int *Nij, nGtree;
   double *fpatt, kappa, alpha, rho, rgene[1], pi[4], piG[1][4];
   double *lnLmax, *locusrate;
   double *pDclass, *tau1beta, *bp0124[3], *wwprior[3];
}  com ;

int LASTROUND, multiplier=100;
double para[8];  /* theta0 theta1 tau0 tau1 qbeta */
int debug = 0;

enum {M0, M1DiscreteBeta, M2SIM3s} MODELS;
enum {G1a, G1b, G1c, G2a, G2b, G2c, G3a, G3b, G3c,  G4a, G4b, G4c, G5a, G5b, G5c,  G6a, G6b, G6c} GTREES;   
char *ModelStr[3] = {"M0", "DiscreteBeta", "SIM3s"};
char *GtreeStr[18] = {"G1a", "G1b", "G1c", "G2a", "G2b", "G2c","G3a", "G3b", "G3c","G4a", "G4b", "G4c","G5a", "G5b", "G5c","G6a", "G6b", "G6c"};

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
   char VerStr[32] = "Version 2.1, July 2012";
   int i;
   double x[8]={1,1,1,1,1}, space[1000];

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
   /*
   ReadSiteCounts(datafile); 
   */
#if(1)
   ReadSeqData(fout, com.seqf, com.ratef, com.cleandata);
   noisy = 3;
   for(i=0; i<4; i++) x[i] = 0.5+0.5*rndu();
   Models0123(fout, frub, frst, x, space);
   free(com.Nij);
#else
   Simulation(fout, frub, space);
#endif

   fclose(frst);
   fclose(frub);
   fclose(fout);
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
   int i,j, h, *n;
   double mr, mNij[5]={0};
   char *names[2]={"ABC", "123"}, *pz[3], ch, newOrder[3];

   printf("\nReading sequence data..  %d loci\n", com.ndata);
   if((com.Nij=(int*)malloc(com.ndata*5*sizeof(int)))==NULL) error2("oom");
   memset(com.Nij, 0, com.ndata*5*sizeof(int));

   if((com.lnLmax=(double*)malloc(com.ndata*(1+com.fix_locusrate)*sizeof(double)))==NULL) 
      error2("oom lnLmax");
   if(com.fix_locusrate) com.locusrate = com.lnLmax + com.ndata;

   for(i=0,n=com.Nij; i<com.ndata; i++,n+=5) {
      fprintf(fout, "\n\n*** Locus %d ***\n", i+1);
      ReadSeq (NULL, fin, cleandata);
      PatternWeightJC69like (fout);

      /* process sequence names */
      for(j=0; j<3; j++) newOrder[j] = j;
      for(j=0; j<3; j++) pz[j] = com.z[j];
      for(j=0; j<3; j++) {
         ch = toupper(com.spname[j][0]);
         newOrder[j] = strchr(names[0], ch) - names[0];
         if(newOrder[j]<0 || newOrder[j]>2)
            newOrder[j] = strchr(names[1], ch) - names[1];
         if(newOrder[j]<0 || newOrder[j]>2)
            newOrder[j] = -9;
      }
      if(newOrder[0]!=newOrder[1] && newOrder[0]!=newOrder[2] && newOrder[1]!=newOrder[2]) {
         for(j=0; j<3; j++)
            com.z[newOrder[j]] = pz[j];
      }

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
      if(noisy>=3) printf("\n");

      for(j=0; j<5; j++)
         mNij[j] += (double)n[j]/(com.ndata*com.ls);
   }
   free(com.fpatt);
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
      com.lnLmax[locus_save] = -lnL*nt - 300;   /* log(10^300) = 690.8 */
   }
   return(0);
}


double lnpD_locus (int locus);
int GetUVrootIM3s(double U[5*5], double V[5*5], double Root[5]);
int GetPMatIM3s(double Pt[5*5], double t, double U[5*5], double V[5*5], double Root[5]);

/* Quadrature using Gauss-Legendre rule.  The following transforms are used.
   x: (-1, 1) <== y: (a=0, b)  <==  t: (0, inf || taugap)
   y = (b+a)/2 + (b-a)/2*x; 
   t = y/(1-y)
   Three G trees are ordered G1c23, G1b, G1a
*/
double lfun (double x[], int np)
{
   int *n, Gtree, locus, itau1, igrid, j, ixw[2], K=com.npoints;
   double theta0, theta1, tau0, tau1, theta12, taugap;
   double p, q, mbeta=-1;  /* paras in beta for tau1 under M1DiscreteBeta */
   double *xI=NULL, *wI=NULL;  /* points & weights from GaussLegendre */
   /* yup[] is y upper bound for y0 & y1 */
   double yup[2]={1,1}, y[2], t[2], b[2], s[2], coeff, Li, lnL=0, z, t1;
   double Pt[5*5], U[5*5], V[5*5], Root[5], PG1a=0;
   double M12, c, m, a, PSG;

   xtoy(x, para, np);
   if(!LASTROUND) para[3] = para[3]*para[2];  /* tau1 */
   for(j=0; j<4; j++) para[j] /= multiplier;  /* qbeta under M1DiscreteBeta is not scaled */
   if(com.model==M2SIM3s) {
      theta12 = para[4] /= multiplier; /* theta1&2 */
      M12 = para[5];
   }
   theta0=para[0], theta1=para[1], tau0=para[2], tau1=para[3]; 
   if(tau0<tau1) error2("tau0<tau1");

   if(com.model==M1DiscreteBeta) {
      mbeta = para[3]/para[2];
      q = para[4] = x[4];
      p = mbeta/(1-mbeta)*q;

      DiscreteBeta(Pt, com.tau1beta, p, q, com.ncatBeta, com.UseMedianBeta);
      for(j=0; j<com.ncatBeta; j++)
         com.tau1beta[j] *= tau0;
   }

   for(itau1=0; itau1<(com.model!=M1DiscreteBeta ? 1 : com.ncatBeta); itau1++) { /* for M1DiscreteBeta */
      if(com.model==M1DiscreteBeta)
         tau1 = para[3] = com.tau1beta[itau1];
      taugap = 2*(tau0 - tau1)/theta1;
      if(com.model==M2SIM3s) {
         c=2/theta12, m=4*M12/theta12, a=sqrt(c*c + 16*m*m);
         PG1a = (c + 4*m + a)*exp(-(c + 4*m - a)*tau1/2)
              - (c + 4*m - a)*exp(-(c + 4*m + a)*tau1/2);
         PG1a = 1 - PG1a/(2*a);

         GetUVrootIM3s(U, V, Root);
         GetPMatIM3s(Pt, tau1, U, V, Root);
         PG1a = Pt[1*5+3] + Pt[1*5+4];
         PSG = (1 - PG1a)*exp(-taugap)*2/3.0;
         /* the code above is for error checking, as PG1a is calculated using 2 methods to 
            confirm that they are identical.  This computation is not expensive
         */
      }

      GaussLegendreRule(&xI, &wI, com.npoints);
      for(Gtree=0; Gtree<com.nGtree; Gtree++) {
         /* coeff includes (b-a)/2 * (b-a)/2 for y0 & y1. */
         if     (Gtree==Gk)   coeff = (1-PG1a)*exp(-taugap)*0.5*0.5;
         else if(Gtree==G1b)  coeff = (1-PG1a)*taugap/(taugap+1)*0.5*0.5; 
         else if(Gtree==G1a)  coeff = (2*tau1/theta12)/(2*tau1/theta12+1)*0.5*0.5;
         if     (Gtree==Gk)   yup[1] = 1;                                    /* for y1 in G1c23 */
         else if(Gtree==G1b)  yup[1] = taugap/(taugap+1);                    /* for y1 in G1b */
         else if(Gtree==G1a)  yup[1] = (2*tau1/theta12)/(2*tau1/theta12+1);  /* for y1 in G1a */
         for(igrid=0; igrid<K*K; igrid++) {
            ixw[0] = igrid/K; ixw[1] = igrid%K;
            for(j=0; j<2; j++) {  /* t0 (y0) and t1 (y1) */
               if(ixw[j]<K/2) { ixw[j] = K/2-1-ixw[j];  s[j]=-1; }
               else           { ixw[j] = ixw[j]-K/2;    s[j]=1;  }
               y[j] = yup[j]*(1 + s[j]*xI[ixw[j]])/2;
               t[j] = y[j]/(1-y[j]);
            }
            com.wwprior[Gtree][igrid] = wI[ixw[0]]*wI[ixw[1]]/square((1-y[0])*(1-y[1]))*coeff;

            /* set up blengths[] & p0124 */
            if(Gtree == Gk) {            /* G1c23 */
               com.wwprior[Gtree][igrid] *= exp(-3*t[1]-t[0]);
               b[0] = theta0*t[0]/2;
               b[1] = tau0 + theta0*t[1]/2;
            }
            else if(Gtree == G1b) {      /* G1b */
               com.wwprior[Gtree][igrid] *= exp(-t[1]-t[0]);
               b[0] = tau0 - tau1 - theta1*t[1]/2 + theta0*t[0]/2;
               b[1] = tau1 + theta1*t[1]/2;
            }
            else if(Gtree == G1a) {       /* G1a */
               t1 = theta12*t[1]/2;
               z = 2*exp(-(c+4*m+a)*t1/2) *(exp(a*t1) - 1) * m / a;  /* 2*P_{123,113}(t1) */

               /*
               GetPMatIM3s(Pt, theta12*t[1]/2, U, V, Root);
               if(fabs(Pt[1*5+0] - Pt[1*5+2]) > 1e-9)
                  error2("the two terms should be the same?");
               if(fabs(Pt[1*5+0] + Pt[1*5+2] - z) > 1e-9) {
                  matout(F0, x, np, np);
                  error2("the two ways should be the same?");
               }
               */
               com.wwprior[Gtree][igrid] *= z*exp(-t[0]);
               b[0] = tau0 - theta12*t[1]/2 + theta0*t[0]/2;
               b[1] = theta12*t[1]/2;
            }

            if(!com.fix_locusrate)
               p0124Fromb0b1(com.bp0124[Gtree]+igrid*5, b);
            else {
               com.bp0124[Gtree][igrid*2+0] = b[0];
               com.bp0124[Gtree][igrid*2+1] = b[1];
            }

            /* printf("%-3s %5d y%9.5f%9.5f t%9.5f%9.5f b%9.5f%9.5f\n",
               GtreeStr[Gtree], igrid+1, y[0],y[1], t[0],t[1], b[0],b[1]); */

         }  /*  for(igrid) */
      }

      for(locus=0; locus<com.ndata; locus++) {
         Li = lnpD_locus(locus);
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

double lnpD_locus (int locus)
{
   int *n = com.Nij + locus*5, K=com.npoints, n123max, igrid, k, Gtree, y, error=0;
   double lnL=0, lmax=com.lnLmax[locus], pD=0, f, p[5], p12, b[2], sump12;
   /* lmax can be dynamically adjusted. */

   for(Gtree=0; Gtree<com.nGtree; Gtree++) {
      for(igrid=0; igrid<K*K; igrid++) {
         if(com.fix_locusrate) {  /* p0124 is calculated only if locus rates are fixed. */
            b[0] = com.bp0124[Gtree][igrid*2+0] * com.locusrate[locus];
            b[1] = com.bp0124[Gtree][igrid*2+1] * com.locusrate[locus];
            p0124Fromb0b1 (p, b);
         }
         else 
            for(k=0; k<5; k++)
               p[k] = com.bp0124[Gtree][igrid*5+k];

         f = -lmax;
         if(n[0]) f += n[0]*log(p[0]);
         if(n[4]) f += n[4]*log(p[4]);
         if(Gtree == Gk) {  /* G1c23 */
            if((k = n[1]+n[2]+n[3]) )
               f += k * log(p[2]);
            n123max = n[1];
            if(n123max < n[2]) n123max = n[2];
            if(n123max < n[3]) n123max = n[3];
            p12 = p[1]/p[2];
            f += n123max*log(p12);
            sump12  = ((y=n[1]-n123max) ? pow(p12,y) : 1);
            sump12 += ((y=n[2]-n123max) ? pow(p12,y) : 1);
            sump12 += ((y=n[3]-n123max) ? pow(p12,y) : 1);
            f = (f < -500 ? 0 : exp(f)*sump12);
         }
         else {      /* G1a & G1b */
            if(n[1]) f += n[1]*log(p[1]);
            if(n[2]+n[3]) f += (n[2]+n[3])*log(p[2]);
            f = (f < -500 ? 0 : exp(f));
         }
         if(debug==1) printf("%-3s %-6d f %9.3e\n", GtreeStr[Gtree], igrid+1, f);
         pD += f*com.wwprior[Gtree][igrid];
      }  /* for(igrid) over the grid */
   }     /* for(Gtree) */
   if(pD < 1e-300) {
      error = -1;
      printf("\nat locus %2d, pD = %.6g\n", locus+1, pD);
      lnL += -1e100 + lmax;
   }
   else
      lnL += log(pD) + lmax;
   if(error) {
       matout(F0, para, 1, 4+(com.model==M1DiscreteBeta)+2*(com.model==M2SIM3s));
       printf("locus %d", locus+1);
       error2("floating point problem");
   }
   return(lnL);
}

int GetUVrootIM3s(double U[5*5], double V[5*5], double Root[5])
{
/* This generates U, V, Root for Q for epoch E1.
*/
   double theta12=para[4], M12=para[5], c=2/theta12, m=4*M12/theta12;
   double a=sqrt(c*c + 16*m*m);
   char *states[] = {"11", "12", "22", "1", "2"};

   U[0*5+0] = 1;  U[0*5+1] = -1;  U[0*5+2] = -1;  U[0*5+3] = 1;  U[0*5+4] = 1; 
   U[1*5+0] = 1;  U[1*5+1] =  0;  U[1*5+2] =  0;  U[1*5+3] = (c-a)/(4*m);  U[1*5+4] = (c+a)/(4*m); 
   U[2*5+0] = 1;  U[2*5+1] =  1;  U[2*5+2] =  1;  U[2*5+3] = 1;  U[2*5+4] = 1; 
   U[3*5+0] = 1;  U[3*5+1] =  0;  U[3*5+2] = -1;  U[3*5+3] = 0;  U[3*5+4] = 0; 
   U[4*5+0] = 1;  U[4*5+1] =  0;  U[4*5+2] =  1;  U[4*5+3] = 0;  U[4*5+4] = 0; 

   V[0*5+0] = 0;          V[0*5+1] =  0;      V[0*5+2] = 0;          V[0*5+3] = 0.5;               V[0*5+4] =  0.5;
   V[1*5+0] = -0.5;       V[1*5+1] =  0;      V[1*5+2] = 0.5;        V[1*5+3] = 0.5;               V[1*5+4] = -0.5;
   V[2*5+0] = 0;          V[2*5+1] =  0;      V[2*5+2] = 0;          V[2*5+3] = -0.5;              V[2*5+4] =  0.5; 
   V[3*5+0] = (1+c/a)/4;  V[3*5+1] = -2*m/a;  V[3*5+2] = (1+c/a)/4;  V[3*5+3] = -(c+a-4*m)/(4*a);  V[3*5+4] = -(c+a-4*m)/(4*a); 
   V[4*5+0] = (1-c/a)/4;  V[4*5+1] =  2*m/a;  V[4*5+2] = (1-c/a)/4;  V[4*5+3] =  (c-a-4*m)/(4*a);  V[4*5+4] =  (c-a-4*m)/(4*a); 
 
   Root[0] = 0;
   Root[1] = -(2*m+c);
   Root[2] = -2*m;
   Root[3] = -(4*m + c + a)/2;
   Root[4] = -(4*m + c - a)/2;
   return(0);

#if(0)
   zero(Q, s*s);
   Q[0*s+1] = 2*m1;     Q[0*s+3] = c1;
   Q[1*s+0] = m2;       Q[1*s+2] = m1;
   Q[2*s+1] = 2*m2;     Q[2*s+4] = c2;
   Q[3*s+4] = m1;
   Q[4*s+3] = m2;
   for(i=0; i<s; i++)
      Q[i*s+i] = -sum(Q+i*s, s);
   matout(F0, Q, s, s);

   if(debug==9) {
      printf("\nQ & P(%.5f)\n", t);
      matout2(F0, Q, s, s, 10, 6);
   }
   matexp(Q, t, s, (t>0.1 ? 31 : 7), space);
   matout(F0, Q, s, s);
#endif

   return(0);
}

int GetPMatIM3s(double Pt[5*5], double t, double U[5*5], double V[5*5], double Root[5])
{
/* P(t) = U * exp(Root*t) U.
*/
   int i, j, k, n=5;
   double expt, uexpt, *pP;

   for (k=0,zero(Pt,n*n); k<n; k++)
      for (i=0,pP=Pt,expt=exp(t*Root[k]); i<n; i++)
         for (j=0,uexpt=U[i*n+k]*expt; j<n; j++)
            *pP++ += uexpt*V[k*n+j];

   return(0);
}


#if(0)
int GenerateQIM3s(double Q[])
{
/* This is not used right now since we only need the Q for a 5-state chain.  
   This generates Q1 and Q2 for epochs E1 and E2.
*/
   double tau0=0.05, tau1=0.025, q1=0.01, q2=0.02, q3=0.04, q4=0.04, q5=0.05, M1=0.01, M2=0.05;
   double m1=M1, m2=M2, c1=2/q1, c2=2/q2, c3=2/q3, c5=2/q5;
   int s1=17, s2=9, i, j;
   char *statesE1[] = {"111", "112", "113", "122", "123", "133", "222", "223", 
      "233", "333", "11", "12", "13", "22", "23", "33", "1", "2", "3"};
   char *statesE2[] = {"333", "335", "355", "555", "33", "35", "55", "3", "5"};
   
   zero(Q1, s1);
   zero(Q2, s2);

   Q1[ 0*s1+ 1] = 3*m1;     Q1[ 0*s1+10] = 3*c1;
   Q1[ 1*s1+ 0] = m2;       Q1[ 1*s1+ 3] = 2*m1;    Q1[ 1*s1+11] = c1;
   Q1[ 2*s1+ 4] = 2*m1;     Q1[ 2*s1+12] = c1;
   Q1[ 3*s1+ 1] = 2*m1;     Q1[ 3*s1+ 6] = m1;      Q1[ 3*s1+11] = c2;
   Q1[ 4*s1+ 2] = m2;       Q1[ 4*s1+ 7] = m1;
   Q1[ 5*s1+ 8] = m1;       Q1[ 5*s1+12] = c3;
   Q1[ 6*s1+ 3] = 3*m2;     Q1[ 6*s1+13] = 3*c2;
   Q1[ 7*s1+ 4] = 2*m2;     Q1[ 7*s1+14] = c2;
   Q1[ 8*s1+ 5] = 2*m2;     Q1[ 8*s1+14] = c3;
   Q1[ 9*s1+15] = 3*c3;
   Q1[10*s1+11] = 2*m1;     Q1[10*s1+16] = c1;
   Q1[11*s1+10] = m2;       Q1[11*s1+13] = m1;
   Q1[12*s1+14] = m1;
   Q1[13*s1+11] = 2*m2;     Q1[13*s1+17] = c2;
   Q1[14*s1+12] = m2;
   Q1[15*s1+18] = c3;
   Q1[16*s1+17] = m1;
   Q1[17*s1+16] = m2;

   for(i=0; i<s1; i++)
      Q1[i*s1+i] = -sum(Q1+i*s1, s1);
   matout2(F0, Q1, s1, s1, 8, 3);

   Q2[0*s2+4] = 3*c3;
   Q2[1*s2+5] = c3;
   Q2[2*s2+5] = c5;
   Q2[3*s2+6] = 3*c5;
   Q2[4*s2+7] = c3;
   Q2[6*s2+8] = c5;

   for(i=0; i<s2; i++)
      Q2[i*s2+i] = -sum(Q2+i*s2, s2);
   matout2(F0, Q2, s2, s2, 8, 3);

   return(0);
}

#endif


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


int GetInitials (int np, double x[], double xb[][2])
{
   int i;
   double thetaU=499, MU=0.15;  /* MU = 0.125 should be fine */

   for(i=0; i<np; i++)  { xb[i][0]=0.001;  xb[i][1]=thetaU; }
   xb[3][1] = 0.999;  /* xtau */

   if(com.model==M0) {
      if(x[0]<0 || x[0]>9) {
         for(i=0; i<3; i++) x[i] = 0.8 + 0.2*rndu();
         x[3] = 0.5 + 0.2*rndu();   /*  xtau  */
      }
   }
   else if(com.model==M1DiscreteBeta) {
      for(i=0; i<3; i++) {
         x[i] *= 0.95 + 0.1*rndu();
         x[i] = max2(0.001, x[i]);
      }
      x[3] *= 0.95 + 0.1*rndu();  /* xtau1 */
      x[4] = 1 + 5*rndu();        /* qbeta */
      xb[4][0] = 0.1;             /* q_beta */
      xb[4][1] = 499;             /* q_beta */
   }
   else if(com.model==M2SIM3s) {
      for(i=0; i<3; i++) {
         x[i] *= 0.95 + 0.1*rndu();
         x[i] = max2(0.001, x[i]);
      }
      x[3] *= 0.95 + 0.1*rndu();                  /* xtau1 */
      x[4] = (x[0] + x[1])/2*(0.9 + 0.2*rndu());  /* theta12 */
      xb[4][0] = 0.001;  xb[4][1] = thetaU;       /* theta12 */ 
      x[5] = 0.005 + 0.1*rndu();                  /* M12 */
      xb[5][0] = 0.0001;                          /* M12 */
      xb[5][1] = MU;
   }

   if(noisy) {
      printf("\nInitials & bounds\n    ");
      if(com.model==M0)     printf("theta0    theta1   tau0 (x100) xtau1\n");
      if(com.model==M1DiscreteBeta) printf("theta0    theta1   tau0 (x100) xtau1     qbeta\n");
      if(com.model==M2SIM3s) printf("theta0    theta1   tau0      tau1  theta1&2      M1&2 (x100)\n"); 
      FOR(i,np) printf(" %9.6f", x[i]); FPN(F0);
      FOR(i,np) printf(" %9.5f", xb[i][0]);  FPN(F0);
      FOR(i,np) printf(" %9.5f", xb[i][1]);  FPN(F0);
   }
   return(0);
}

double Models0123 (FILE *fout, FILE *frub, FILE *frst, double x[], double space[])
{
   int np, i,j, s, noisy0=noisy, K=com.npoints;
   char timestr[96];
   double *var, lnL, lnL0=0, e=1e-8;
   double xb[6][2];

   com.pDclass = (double*)malloc((com.ncatBeta*com.ndata+com.ncatBeta)*sizeof(double));
   if(com.pDclass==NULL) error2("oom Models01231");
   com.tau1beta = com.pDclass + com.ncatBeta*com.ndata;
   s = (com.fix_locusrate ? K*K*2 : K*K*5);  /* 2 for b0 & b1; 5 for p0124 */
   com.bp0124[0] = (double*)malloc(3*s*sizeof(double));
   com.bp0124[1] = com.bp0124[0] + s;
   com.bp0124[2] = com.bp0124[1] + s;
   com.wwprior[0] = (double*)malloc(3*K*K*sizeof(double));
   com.wwprior[1] = com.wwprior[0] + K*K;
   com.wwprior[2] = com.wwprior[1] + K*K;
   if(com.bp0124[0]==NULL || com.wwprior[0]==NULL) 
      error2("oom Models01232");

   noisy = 0;
   Initialize3s(space);
   noisy = noisy0;

   for(com.model=0; com.model<3; com.model++) {
      if(com.model==M0)                    { np=4;  com.nGtree=2;  }
      else if(com.model==M1DiscreteBeta)   { np=5;  com.nGtree=2;  }
      else if(com.model==M2SIM3s)          { np=6;  com.nGtree=3;  }

      LASTROUND = 0;
      if(noisy) printf("\n\n*** Model %d (%s) ***\n", com.model, ModelStr[com.model]);
      if(fout) {
         fprintf(fout, "\n\n*** Model %d (%s) ***\n", com.model, ModelStr[com.model]);
         fprintf(frub, "\n\n*** Model %d (%s) ***\n", com.model, ModelStr[com.model]);
      }

      GetInitials (np, x, xb);

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
         for(ii=0; ii<nii; ii++) {
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
      lnL = lfun(x, np);
      if(noisy) printf("\nlnL0 = %12.6f\n", -lnL);

      NFunCall=0;
      ming2(frub, &lnL, lfun, NULL, x, xb, space, e, np);
      
      if(com.model==M0) lnL0 = lnL;
      LASTROUND = 1;
      x[3] *= x[2];     /* xtau1 -> tau1 */
      if(noisy) {
         printf("\nlnL  = %12.6f %+12.6f\nMLEs\n    ", -lnL, lnL0-lnL);
         if(com.model==M0)      printf("theta0    theta1      tau0      tau1 (x100)\n");
         if(com.model==M1DiscreteBeta)  printf("theta0    theta1      tau0      tau1 (x100) qbeta\n");
         if(com.model==M2SIM3s) printf("theta0    theta1      tau0      tau1  theta1&2 (x100) M1&2\n"); 
         for(i=0; i<np; i++)    printf(" %9.6f", x[i]);
      }
      if(fout) {
         fprintf (fout, "\nlnL  = %12.6f %+12.6f\nMLEs\n    ", -lnL, lnL0-lnL);
         if(com.model==M0)      fprintf(fout, "theta0    theta1      tau0      tau1 (x100)\n");
         if(com.model==M1DiscreteBeta)  fprintf(fout, "theta0    theta1      tau0      tau1 (x100) qbeta\n");
         if(com.model==M2SIM3s) fprintf(fout, "theta0    theta1      tau0      tau1  theta1&2 (x100) M1&2\n");
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
         fprintf(frst, " %9.6f\n", lnL0-lnL);
      }

      if(com.model==M2SIM3s) {
         for(i=0; i<np; i++) fprintf(frst, "\t%.6f", x[i]);  
         printf("\nThe other local peak is at\n");
         fprintf(fout, "\nThe other local peak is at\n");
         x[4] = x[4]/(8*x[5]);  /*  theta/(8M)  */
         x[5] = 1/(64*x[5]);    /*  1/(64M)  */
         lnL = lfun(x, np);
         for(i=0; i<np; i++) printf(" %9.6f", x[i]);  
         printf(" %12.6f\n", -lnL);
         for(i=0; i<np; i++) fprintf(fout, " %9.6f", x[i]);  
         fprintf(fout, " %12.6f\n", -lnL);
         fprintf(frst, "\t%.6f\t%.6f\t%.6f\n", lnL0-lnL, x[4], x[5]);
      }
      fflush(frst);

      printf("\nTime used: %s\n", printtime(timestr));
      if(com.model==0)  x[3] /= x[2];     /* tau1 -> xtau1, as initial for M1 */
   }  /* for(com.model) */

   free(com.pDclass);
   free(com.bp0124[0]);
   free(com.wwprior[0]);
   return(lnL0-lnL);
}
