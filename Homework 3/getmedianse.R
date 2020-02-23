
### Below is an R function for computing a standard error estimate of the median survival time (using right-censored data)
### based on bootstrap resampling.
###
### Input:   'survobj' = survival object
###          'bootruns' = number of replications performed in bootstrap scheme
###
### Output:  numerical element 'se' = estimate of standard error;


library(survival)

getmedianse = function(survobj,bootruns=2000){
	
	n = length(survobj[,1]);
	
	computemedian = function(dat,indices){
	
		fit = survfit(dat[indices,]~1)
		med = as.numeric(summary(fit)$table[7]);
		obs = 1;
		if(is.na(med)){
			med = max(fit$time);
			obs = 0;
		}	
		return(c(med,obs))
	
	}
	
	medvals = NULL; obsvals = NULL
	for(i in 1:bootruns){
	
		indices = sample(1:n,replace=TRUE);
		c = computemedian(survobj,indices);
		medvals[i] = c[1]; obsvals[i] = c[2];
		
	}
	
	medsurv = survfit(Surv(medvals,obsvals)~1);
	jumplocs = medsurv$time; jumpsizes = c(1,medsurv$surv[-length(medsurv$surv)])-medsurv$surv;
	stderr = sqrt(sum(jumplocs^2*jumpsizes)-(sum(jumplocs*jumpsizes))^2);
	
	return(stderr)
	
}
