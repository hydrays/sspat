
double lfun (double x[], int np)
{
	
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
	
  //printf("\n here is 2.5.2\n"); //hhh
  //zzz: compute Q  
  for (i=0;i<DIM*DIM;i++){
    if(mc1[i]) Q[i]= c1 * mc1[i];
    if(mc2[i]) Q[i]+=c2 * mc2[i];
    if(mw1[i]) Q[i]+=m1 * mw1[i];
    if(mw2[i]) Q[i]+=m2 * mw2[i];
  }
  //for(i=0;i<DIM*DIM;i++){printf("4.3f",Q[i]); if(i%DIM==DIM-1) printf("\n");}
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
	com.wwprior[Gtree][igrid][0] = coeff*wI[ixw[0]]*wI[ixw[1]]/square((1-y[0])*(1-y[1])); 
			
	//printf("\n I'm at bbb\n"); //hhh
			
	if(Gtree==G1a||Gtree==G1b || Gtree==G1c){
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
		f[s-1]=Pt1[get_ij(s,1,DIM)]*2/theta1* p1 + Pt1[get_ij(s,2,DIM)]*2/theta1*p2
		  +Pt1[get_ij(s,7,DIM)]*2/theta2*p3 + Pt1[get_ij(s,8,DIM)]*2/theta2*p4;
	    }
	    else if (Gtree==G1b){
	      p1 = 2/theta1*Pt0[get_ij(10,10,DIM)]+2/theta2*Pt0[get_ij(10,13,DIM)];
	      p2 = 2/theta1*Pt0[get_ij(16,10,DIM)]+2/theta2*Pt0[get_ij(16,13,DIM)];
	      p3 = 2/theta1*Pt0[get_ij(19,10,DIM)]+2/theta2*Pt0[get_ij(19,13,DIM)];
	      p4 = 2/theta1*Pt0[get_ij(13,10,DIM)]+2/theta2*Pt0[get_ij(13,13,DIM)];
	      for(s=1;s<=8;s++)
		f[s-1]=Pt1[get_ij(s,1,DIM)]*2/theta1*p1 + Pt1[get_ij(s,6,DIM)]*2/theta2*p2
		  +Pt1[get_ij(s,3,DIM)]*2/theta1*p3 + Pt1[get_ij(s,8,DIM)]*2/theta2*p4;
	    }
	    else if (Gtree==G1a){
	      p1=2/theta1*Pt0[get_ij(9, 9,DIM)]+2/theta2*Pt0[get_ij(9, 12,DIM)];
	      p2=2/theta1*Pt0[get_ij(15,9,DIM)]+2/theta2*Pt0[get_ij(15,12,DIM)];
	      p3=2/theta1*Pt0[get_ij(18,9,DIM)]+2/theta2*Pt0[get_ij(18,12,DIM)];
	      p4=2/theta1*Pt0[get_ij(12,9,DIM)]+2/theta2*Pt0[get_ij(12,12,DIM)];
	      for(s=1;s<=8;s++)
		f[s-1]=Pt1[get_ij(s,1,DIM)]*2/theta1*p1 + Pt1[get_ij(s,4,DIM)]*2/theta2*p2
		  +Pt1[get_ij(s,5,DIM)]*2/theta1*p3 + Pt1[get_ij(s,8,DIM)]*2/theta2*p4;
	    }// if gtree==G1a
	  } /*com.model=M0)*/
				
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
	  if(Gtree==G2a || Gtree==G2b || Gtree==G2c){
	    g=exp(-t[0]);
	    for(s=1;s<=8;s++)
	      f[s-1]*=g;
	    b[1]=t[1];
	    b[0]= 2.0*t[0]/theta5+tau1-t[1];
	  }
	  else if(Gtree==G3a || Gtree==G3b || Gtree==G3c){
	    g= exp(-t[0])*exp(-2*(tau0-tau1)/theta5);
	    for(s=1;s<=8;s++)
	      f[s-1]*=g;
	    b[1]=t[1];
	    b[0]= 2.0*t[0]/theta4+tau0-t[1];
	  }
	}

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
	}
    
	for(s=0;s<8;s++)
	  {
	    com.wwprior[Gtree][igrid][s] = f[s] * com.wwprior[Gtree][igrid][0];
	    //printf("wwprior[%d][%d][%d]=%f  ",Gtree, igrid,s,com.wwprior[Gtree][igrid][s]);
	  }
	//printf("\n"); //hhh
	if(!com.fix_locusrate)
	  p0124Fromb0b1(com.bp0124[Gtree]+igrid*5, b);
	else {
	  com.bp0124[Gtree][igrid*2+0] = b[0];
	  com.bp0124[Gtree][igrid*2+1] = b[1];
	}
	/* printf("%-3s %5d y%9.5f%9.5f t%9.5f%9.5f b%9.5f%9.5f\n",
	   GtreeStr[Gtree], igrid+1, y[0],y[1], t[0],t[1], b[0],b[1]); */
	//printf("In lfun after igrid loop !!!!\n"); zzz
      }
	
    for(locus=0; locus<com.ndata; locus++) {
      Li = lnpD_locus(locus,com.state[locus]);
      if(com.model==M1DiscreteBeta)
	com.pDclass[itau1*com.ndata+locus] = Li;
      else
	lnL += Li;
      n = com.Nij + locus*5;
      if(debug) printf("%d\t%d\t%d\t%d\t%d\t%.6f\n", n[0],n[1],n[2],n[3],n[4], Li);
    }
	
    return(-lnL);
  }
