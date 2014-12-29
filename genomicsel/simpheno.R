#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),200,1000)
for (i in 1:200) {
  M[i,] <- ifelse(runif(1000)<0.5,-1,1)
}

simpheno <- function(geno, grandmean = 0, markeffdist = rnorm(n=ncol(geno), 0,0.05), residualdist = rnorm(n=nrow(geno), 0,1), seed=123, plot=TRUE ){
#random phenotypes
set.seed(seed)
u <- markeffdist
g <- as.vector(crossprod(t(geno),u))
y <- grandmean + g + residualdist
genovar <- var(g)
print(paste("genotypic variance = ", round(genovar,2), sep=""))
print(paste("error variance = ", round(var(residualdist),2), sep=""))
print(paste("heritability = ", round(var(g)/ var(y),2)))
if(plot){
  par(mfrow=c(2,2))
  hist(u, main = "histogram of marker effects", col="pink")
  hist(g, main = "histogram of genotypic effects", col="red")
  hist(residualdist, main="histogram of residual", col="blue")
  hist(y, main= "histogram of resulting phenotype", col="green")
}
return(y)
}
simpheno(M, 50, markeffdist= rnorm(n=ncol(M), 0,0.05))