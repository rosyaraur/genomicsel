\name{creatGPO}
\alias{creatGPO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creating dataset of class gpo 
}
\description{
As first step the function creates a dataset combining all of the input data required or available to perform further computations. 
}
\usage{
creatGPO(Ymat, Mmat, Xmat = NULL, Kmat = NULL, map = NULL, minMAF = 0, maxMISS = 0, thinINT = NULL, pRandSample = NULL, seed = NULL, expected.geno = -4:4, expected.homo = c("AA", "BB", "CC", "GG", "TT"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ymat}{
Y dataset matrix (should be matrix with named rows (indivudals / subjects / genotypes / accessions) and columns (phenotypic data / trait data name) - required 
}
  \item{Mmat}{
Marker dataset matrix (should be matrix with named rows (indivudals / subjects / genotypes / accessions) and columns (marker data)
The marker dataset can be both numeric or character. If character dataset is provided the function will convert it to numeric value - required 
}
  \item{Xmat}{
X dataset matrix - additional covariates (should be matrix with named rows (indivudals / subjects / genotypes / accessions) and columns (phenotypic data / trait data name) - optional argument 
}
  \item{Kmat}{
Relatedness / kinship matrix, should be named, the name of individuals in Xmat, Mmat and Ymat should match if available (rownames(Ymat)==rownames(Xmat)== rownames(Mmat) ) - optional argument 
}
  \item{map}{
The map should have at least three columns in the following order - Marker name, chromosome (integer), position (integer).
The column names of marker matrix should match with marker names in the map. Optional argument, only required in some cases for example if thinning of markers need to be done by distance. 
}
  \item{minMAF}{
Minimum Minor Allele Frequency (MAF) - marker filetering crieteria
}
  \item{maxMISS}{
Maximum missining value allowed - marker filtering criteria 
}
  \item{thinINT}{
Marker thinning interval - default is NULL means no thining of marker will be done and all markers will be used (that meets the other criteria as minMAF, maxMISS) 
}
  \item{pRandSample}{
proportional of marker to be Randomly selected - default is NULL means no selection 
}
  \item{seed}{
set seed 
}
  \item{expected.geno}{
Expected genotype (integer) - if the provided matrix contains numerically coded for example AA=-1, AB=0, BB=1 in diploids or AAAA=-2, AAAB=-1, AABB=0, ABBB=1, BBBB=2. 
}
  \item{expected.homo}{
The expected homozygote genotypes - if the recoding of marker matrix need to be done. 
}
}
\details{
Naming of row and column of matrix is very import to work. 
}
\value{
The function creat a data object of class "gpo" which is required by most of other functions. The object consists of 
Ymat, Mmat, Xmat, Kmat, map and other arguments supplied. 
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
rm(list=ls())
#source("C:/Users/umesh.rosyara/SkyDrive/Documents/genomicsel/readallfiles.r")
#source("/Users/umesh.rosyara/Dropbox/genomicsel/readallfiles1.r")
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
simdataGPO <- creatGPO (Ymat=simdata$Ymat, Mmat=simdata$Mmat, Xmat=NULL,Kmat=simdata$Kmat, 
map = simdata$map,  minMAF=0.0, maxMISS = 0, thinINT = NULL, pRandSample = NULL, 
seed=NULL, expected.geno = -1:1)
is(simdataGPO)

}

