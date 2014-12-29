\name{doPreProcess}
\alias{doPreProcess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preprocessing of the dataset 
}
\description{
The function can be useful for transforming and imputing marker data, transforming y variables and covariates. 
}
\usage{
doPreProcess(gpobj, impute.M = "meanImpute", impute.X = FALSE, transform.X = c("scale", "center", "pca"), transform.M = FALSE, pca.threshold = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpobj}{
data set gpo object 
}
  \item{impute.M}{
Marker imputing method - available methods: "knnImpute", "rfImpute", "medianImpute", "meanImpute","modeImpute"
}
  \item{impute.X}{
logical if we need impute X ????  
}
  \item{transform.X}{
Transform covariate - available methods: "scale", "center", "pca"
}
  \item{transform.M}{
logical, Transforming marker data using PCA
}
  \item{pca.threshold}{
PCA threshold is transform.M is TRUE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Umesh R Rosyara 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(wheat_tri)
#source("C:/Users/umesh.rosyara/SkyDrive/Documents/genomicsel/readallfiles.r")

pheno <- wheat_tri$pheno
y.name <- wheat_tri$pheno[,1]
pheno <- as.matrix(pheno[,c(-1,-2, -4)], ncol=1)
rownames(pheno) <- y.name
colnames(pheno) <- "pht"

geno <- wheat_tri$geno[,-2860]
rownames(geno) <- wheat_tri$geno[,1]
geno <- geno[,-1]


wheatTT <- creatGPO (Ymat=pheno, Mmat=geno, Xmat=NULL,Kmat=NULL, map = NULL,  minMAF=0, minMISS = 0, thinINT = NULL, pRandSample = NULL, seed=NULL, 
expected.geno = -1:1)


wheatTT1 <- doPreProcess (wheatTT, impute.M = "meanImpute", transform.Y=NULL,  transform.X = NULL)

lst <- doDataPartition(wheatTT1, times = 1, partitions=2, p = c(0.8), seed = c(45))

data_tr <- sepPartition (lst, rep=1, part=1)
test_data <- sepPartition (lst, rep=1, part=2)

outTT <- fit.model(gpobj=data_tr, module = "rrBLUP", method = "rrBLUP", argument_list=NULL)
predtsTT <- predict(outTT, newM=test_data$Mmat, newX=NULL)
cor(predtsTT,test_data$Ymat)
plot(predtsTT,test_data$Ymat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
