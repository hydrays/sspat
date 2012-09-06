'twotypes' <- function(T, r1, r2, p, w){
	   ngood <- 1
	   nbad <- 0
	   for (i in seq(T)){
	       ngood <- ngood*2
	       r <- runif(ngood)
	       ngood_new <- ngood
	       for (j in seq(ngood)){
	       	   if(ngood>0){
			if ( r[j] < p ){ 
		           ngood_new <- ngood_new -1
		           nbad <- nbad + 1
			}
		   }		      
	       }
	       ngood <- ngood_new
	       cat(i, ngood, nbad, '\n')
	   }
	   return(ngood+nbad)
}     
	       