
double lfun (double x[], int np)
{
  GaussLegendreRule(&xI, &wI, com.npoints);
  for (Gtree=0; Gtree<com.nGtree;Gtree++){
    for(igrid=0; igrid<K*K; igrid++){
      if(Gtree==G1a||Gtree==G1b || Gtree==G1c){
	if(com.model==M2SIM3s){
	  xtoy(Q,Q1,DIM*DIM);
	  complexroots (Q1, Pt1, t[0]*t[1]);
	  if(0 && igrid==0){
	    debug_print(Pt1, DIM);
	  }
	  xtoy(Q,Q1,DIM*DIM);
	  complexroots (Q1, Pt0, t[0]-t[0]*t[1]);
	  if(0 && igrid==0){
	    debug_print(Pt0, DIM);
	  }
	  if (Gtree==G1c){
	  }
	  else if (Gtree==G1b){
	  }
	  else if (Gtree==G1a){
	  }
	}
	b[1]=t[0]*t[1];	
	b[0]=t[0]-t[0]*t[1];
      }
      else if(Gtree==G2a || Gtree==G2b || Gtree==G2c || Gtree== G3a || Gtree ==G3b || Gtree ==G3c){
	if(com.model==M2SIM3s){
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
