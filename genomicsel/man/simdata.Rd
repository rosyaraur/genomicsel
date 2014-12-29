\name{simdata}
\alias{simdata}
\docType{data}
\title{
Simulated dataset 
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("simdata")}
\format{

The dataset consists of the following components: 

List of 5

 $ Ymat: num [1:200, 1] y variable with 200 individual values 
 
  ..- attr(*, "dimnames")=List of 2
  
  .. ..$ : chr [1:200] individual ids
  
  .. ..$ : chr "y1" y value 
  
 $ Mmat: num [1:200, 1:1000] marker matrix of 200 individuals and 1000 markers 
 
  ..- attr(*, "dimnames")=List of 2
  
  .. ..$ : chr [1:200] individuals "id1" "id2" "id3" "id4" ...
  
  .. ..$ : chr [1:1000] markers "mrk1" "mrk2" "mrk3" "mrk4" ...
  
 $ Xmat: num [1:200, 1:2] covariate matrix of two covariates and 200 individuals 
 
  ..- attr(*, "dimnames")=List of 2
  
  .. ..$ : chr [1:200] indiviudal ids "id1" "id2" "id3" "id4" ...
  
  .. ..$ : chr [1:2]  two covariate variables "x1" "x2"
  
 $ Kmat: num [1:200, 1:200] K matrix  200 individuals x 200 individuals ...
 
  ..- attr(*, "dimnames")=List of 2
  
  .. ..$ : chr [1:200] individuals "id1" "id2" "id3" "id4" ...
  
  .. ..$ : chr [1:200] individuals "id1" "id2" "id3" "id4" ...
  
 $ map :'data.frame':	1000 obs. of  3 variables: map with 3 coumns with marker namode1l <- fit.model (data_tr, method = "lasso", module = "elesticnet",  argument_list= NULL)
plot(mode1l)
trainp <- doTrain(data_tr, method = "lasso", module = "enet",  validation = "cv", k = 3, times =1, seed=NULL, argument_list= NULL)
plot(trainp)me, chromsome, position in genome 
 
  ..$ mark   : Factor w/ 1000 levels "mrk1","mrk10",..:  marker names 
  
  ..$ chr    : int [1:1000] chromsomes...
  
  ..$ postion: num [1:1000] position in genome in cM ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
rm(list=ls())
#source("C:/Users/umesh.rosyara/SkyDrive/Documents/genomicsel/readallfiles.r")
ls()
data(simdata)
str(simdata)

# rows and columns of evary dataset must be named.
# phenotypic data 
rownames(simdata$Ymat)
colnames(simdata$Ymat)

# genotypic data 
rownames(simdata$Mmat)
colnames(simdata$Mmat)[1:10]

# covariates 
rownames(simdata$Xmat)
colnames(simdata$Xmat)

# Kinship / relatedness matrix 
rownames(simdata$Kmat)
colnames(simdata$Kmat)

# creat GPO object 
out2 <- creatGPO (Ymat=simdata$Ymat, Mmat=simdata$Mmat, Xmat=NULL,Kmat=simdata$Kmat, 
map = simdata$map,  minMAF=0.0, minMISS = 0.1, thinINT = NULL, pRandSample = NULL, 
seed=NULL, expected.geno = -1:1)

is(out2)

# split data into training and test data 
lst <- doDataPartition(out2, times = 1, partitions=2, p = c(0.8), seed = c(45))
data_tr <- sepPartition (lst, rep=1, part=1)
test_data <- sepPartition (lst, rep=1, part=2)

# fitting model using rrBLUP module and rrBLUP method 
out <- fit.model(gpobj=data_tr, module = "rrBLUP", method = "rrBLUP", argument_list=NULL)
predts <- predict(out, newM=test_data$Mmat, newX=NULL)
cor(predts,test_data$Ymat)
plot(predts,test_data$Ymat)

trainRR <- doTrain(data_tr, module = "rrBLUP", method = "rrBLUP",
validation = "cv", k = 10, times =1, seed=NULL, argument_list= NULL)
plot(trainRR)


# random forest 
out <- fit.model(gpobj=data_tr, module = "randomForest", method = "randomForest", argument_list=NULL)
pdts <- predict(out, test_data$Mmat)
cor(pdts,test_data$Ymat)
plot(pdts,test_data$Ymat)

# perform training 
outtrain <- doTrain(data_tr, method = "randomForest", module = "randomForest",  validation = "cv", k = 10, times =1, seed=NULL, argument_list= NULL)
plot(outtrain)

# fit BGLR model 
mod_brr <- fit.model (data_tr, module = "BGLR", method = "BRR", argument_list=NULL)
mod_brr$r.squared 

mod1_baya <- fit.model (data_tr, module = "BGLR", method = "BayesA", argument_list=NULL)
mod1_baya$r.squared

mod2_bayc <- fit.model (data_tr, module = "BGLR", method = "BayesC", argument_list=NULL)
mod2_bayc$r.squared 
pdts <- predict (mod2_bayc, test_data$Mmat)
pdts

cor(pdts,test_data$Ymat)
plot(pdts,test_data$Ymat)

results_bayc <- doTrain(data_tr, method = "BayesC", module = "BGLR",  validation = "cv", k = 3, times =1, seed=NULL, argument_list= NULL)
plot(results_bayc)


# fit elestric net 
#pcr method 
mod_pcr <- fit.model (data_tr, method = "pcr", module = "pls",  argument_list= NULL)
plot(mod_pcr)
str(mod_pcr)

pc_model <- doTrain (data_tr, method = "pcr", module = "pls", validation = "cv", k = 10, times =1, seed=NULL, argument_list= NULL)
#Error in predict.mvr(fitobj, ncomp = ncomp, newdata = xmat) : 
#  'newdata' does not have the correct number of columns  
plot(pc_model)

tm_pls <- doTrain (data_tr, method = "pls", module = "pls", validation = "cv", k = 3, times =1, seed=NULL, argument_list= NULL)
plot(tm_pls)

tm_cppls <- doTrain (data_tr, method = "cppls", module = "pls", validation = "cv", k = 3, times =1, seed=NULL, argument_list= NULL)
plot(tm_cppls)

# Error in solve.default(crossprod(P[, 1:a, drop = FALSE], W[, 1:a, drop = FALSE])) : system is computationally singular: reciprocal condition number = 6.17308e-17

# method Foba 
mod_foba <- fit.model (data_tr, method = "foba", module = "foba",  argument_list= NULL)
tm_foba <- doTrain (data_tr, method = "foba", module = "foba", validation = "cv", k = 4, times =1, seed=NULL, argument_list= NULL)
plot(tm_foba)


# svm in kernlab package 
  
mod_svm <- fit.model (data_tr, method = "svm", module ="kernlab" ,  argument_list= NULL)
mod_svm

out_svm <- doTrain (data_tr, method = "svm", module ="kernlab", validation = "cv", k = 3, times =1, seed=NULL, argument_list= NULL)
plot(out_svm)




}
\keyword{datasets}
