esemblePrediction <-
function (listpred, y=NULL, method = "average", prop=NULL){
	
mat <- sapply(listpred, cbind)
	
if(method == "average"){
		avgpred <- apply(mat, 1, mean)
}

if(method == "proportions"){
		if(length(prop) ==0){
		# find best prop from the data 
		if(length(y)==0){
			stop("need actual y values to determine proportions, if you do not have use method = average")
		}
		prop.s <- seq(0,10)
		n.mod <- ncol(mat)
		n.proprow <- 1000000
        pgridmat <- matrix(sample(prop.s, n.mod*n.proprow, replace = TRUE), ncol=n.mod)
		maxp <- apply(pgridmat,1, sum)
		prop.grid <- round (pgridmat/maxp,1)
		prop.grid[,n.mod] <- 1- apply(prop.grid[,-n.mod], 1, sum)
		RMSE.esemble.p <- rep (NA, n.proprow)
		for (i in 1:n.proprow){
			propi <- prop.grid[i,]
			avgpredi <- mat %*% propi 
			yhat = as.vector (avgpredi)
			y <- as.vector (y)
			residual_i <- y -  yhat 
			RMSE.esemble.p[i] = sqrt( 1/length(y) *(sum((residual_i)^2)))
		}
		min.rmse.row <- which.min(RMSE.esemble.p)
		prop <- prop.grid[min.rmse.row, ]
		}
		avgpred <- mat %*% prop 
}
	
if(length(y) != 0){	
	yhat <- as.vector(avgpred)
	residual.all <- mat-yhat 
	
	RMSE.all <- rep(NA, ncol(residual.all))
	for (j in 1:ncol(residual.all)){
	RMSE.all[j] = sqrt( 1/length(y) *(sum((residual.all[,j])^2)))
}

yhat = avgpred 
residual <- y -  yhat 
RMSE.esemble = sqrt( 1/length(y) *(sum((residual)^2)))
return(list(RMSE.models = RMSE.all, RMSE.esemble=RMSE.esemble, prop=prop))
} else{
return(avgpred)	
}
}
