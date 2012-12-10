double lfun (double x[], int np)
{
  GaussLegendreRule(&xI, &wI, com.npoints);
  for (Gtree=0; Gtree<com.nGtree;Gtree++){
    for(igrid=0; igrid<K*K; igrid++){

      if(Gtree==G1a||Gtree==G1b || Gtree==G1c){
	if(com.model==M2SIM3s){
	  if(0 && igrid==0){	  }
	  if(0 && igrid==0){	  }
	  if (Gtree==G1c){	  }
	  else if (Gtree==G1b){	  }
	  else if (Gtree==G1a){	  }
	}
      }
      else if(Gtree==G2a || Gtree==G2b || Gtree==G2c || Gtree== G3a || Gtree ==G3b || Gtree ==G3c){
	if(com.model==M2SIM3s){
	  for(s=1;s<=8;s++){	  }
	}
	if(Gtree==G2a || Gtree==G2b || Gtree==G2c){	}
	else if(Gtree==G3a || Gtree==G3b || Gtree==G3c){	}
      }
      else{
	if (Gtree==G4a || Gtree==G4b || Gtree==G4c){      }
	else if (Gtree==G5a || Gtree==G5b || Gtree==G5c){	}
	else if (Gtree==G6a || Gtree==G6b || Gtree==G6c){	}
	if(com.model==M2SIM3s){	}
      }

      for(s=0;s<8;s++){      }

      if(!com.fix_locusrate)
	p0124Fromb0b1(com.bp0124[Gtree]+igrid*5, b);
      else{      }
    }
  }
  for(locus=0; locus<com.ndata; locus++) {    }
  return(-lnL);
}
