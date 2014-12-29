creatGPO <-
function (Ymat, Mmat, Xmat=NULL,Kmat=NULL, map = NULL,  minMAF=0, 
maxMISS = 0, thinINT = NULL, pRandSample = NULL, seed=NULL, impute.M = NULL, 
expected.geno = -4:4, expected.homo= c("AA", "BB", "CC", "GG",  "TT")){

if(length(Ymat)==0 | length(Mmat) ==0){
    stop("Please provide valid Ymat and Mmat")
	}
if(length(rownames(Ymat)) ==0  | length(rownames(Mmat))==0 | length(rownames(Mmat)) == 0){
     stop("Ymat, Mmat - row and columns must be named")
	 }
if(length(Xmat) != 0){
    if(length(rownames(Xmat)) ==0 | length(colnames(Xmat))==0){
	stop("please provide rownames and column names to Xmat")
	}
	}
if(length(Kmat) != 0){
    if(length(rownames(Kmat)) ==0 | length(colnames(Kmat)) ==0){
	stop("please provide rownames and column names to Kmat")
	}
	}	
#Mmat should have markers in rows
# ids in the column
trait.n <- colnames(Ymat)
# input marker data check 
all.genotype <- unique(c(Mmat))
M <- t(Mmat) 
if(is.numeric(all.genotype)){
num.genotypes = c(expected.geno, NA)
if(!all(all.genotype %in% num.genotypes)){
print(paste("Only allowed genotype characters (unless imputed with mean) :"))
print(expected.geno)
}
} else{
tmp1 <- matrix(NA,nrow(M),ncol(M))
for (i in 1:ncol (M)) {
all.gn <- unique (M[,i]) 
homo.geno <- expected.homo 
homo.ref <- homo.geno[which( homo.geno %in% all.gn)]
ref.al1 = as.character (substr(homo.ref, 1,1)[1])
count.al <- function(ref.al, s) {
    s[s=="NA"] <- NA
	s[s=="-"] <- NA
	s[s=="--"] <- NA
	s[s=="."] <- NA
	si <- s
	si[!is.na(si)] <- 1
    s2 <- gsub(ref.al,"",s, fixed = TRUE)
	s3 = nchar(s) - nchar(s2)
    return (s3 * as.numeric(si))
  }
x <- M[,i]
tmp1[,i] <- count.al(ref.al = ref.al1, s=x)
}
M <- tmp1
}

rownames(M) <- colnames(Mmat)
colnames(M) <- rownames(Mmat)

Mmat <- M

# remove missing value from Ymat 
Ymat <- na.omit(Ymat)

if(length(Xmat) != 0){
Xmat <- na.omit(Xmat)
}
print("missing values from Ymat or Xmat are removed")

trait.n <- colnames(Ymat)
y.id.name <- rownames (Ymat)
m.id.name <- colnames (Mmat)

k.id.name <- rownames(Kmat)
y.m <- y.id.name[y.id.name%in%m.id.name]

# ids both in M and Y 
Ymat <- matrix(Ymat[y.m,],nrow=length(y.m))
rownames(Ymat)<- y.m
colnames(Ymat) <- trait.n 
Mmat <- t(Mmat)
Mmat <- Mmat[y.m,]

if(length(Xmat)!=0){
y.id.name <- rownames (Ymat)
x.id.name <- rownames (Xmat)
x.m <- y.id.name[y.id.name%in%x.id.name]
# ids both in M and Y and X
Ymat <- matrix(Ymat[x.m,], nrow=length(x.m))
rownames(Ymat) <- x.m
colnames(Ymat) <- trait.n 
Mmat <- Mmat[x.m,]
Xmat <- Xmat[x.m,]
}

if(length(Kmat)!=0){
y.id.name <- rownames (Ymat)
k.id.name <- rownames (Kmat)
# ids both in M, Y, X and K 
k.m <- y.id.name[y.id.name%in%k.id.name]
Ymat <- matrix (Ymat[k.m,], nrow=length(k.m))
rownames(Ymat) <- k.m
colnames(Ymat) <- trait.n 
Mmat <- Mmat[k.m,]
Xmat <- Xmat[k.m,]
Kmat <- Kmat[k.m,k.m]
}
if(maxMISS > 0){
# filting markers with high missing value 
miss.v <- apply(Mmat, 2, function(x) sum(is.na(x))) / nrow(Mmat) 
miss.col <- c(which(miss.v > maxMISS))

if(length(miss.col) != 0){
Mmat <- Mmat[,-miss.col]
}
}
if(minMAF > 0){
# filtering markers with low minor allele frequency 
maf.fun <- function(x){
freq <- mean(x + 1, na.rm = TRUE)/2
 MAF <- min(freq, 1 - freq)
 return(MAF)
 }
lmaf.v <- apply(Mmat, 2, maf.fun)
lmaf.col <- c(which(lmaf.v < minMAF))

if(length(lmaf.col) != 0){
Mmat <- Mmat[,-lmaf.col]
}
}
 
# thinning 
if(length(map) !=0){
if(length(thinINT)!=0){
if(nrow(map) < 3){
   stop("Physical / linkage map is required and should have three columns: marker, chr, position")
   }
# sort the map 
map <- map[order(map[, 2], map[, 3]), ]

# check and tallary map with the M matrix 
m.map.name <- map[,1]
m.mat.name <- colnames (Mmat)   
mp.m <- m.map.name[m.map.name%in%m.mat.name]
Mmat <- Mmat[,mp.m]
map <- map [map[,1]== mp.m,] 
nchr <- length (unique (map[,2]))
# thinning the map
for ( i in 1:nchr){
mapi <- map [map[,2]==i,]
mapii <- mapi
slI <- rep (NA, nrow(mapi))
slI[1] <- mapi[1,3]
d1 <- mapi[1,3]
wl = nrow(mapi)

for (j in 1:wl){
d2 <- d1 + thinINT
mapip <- mapi[mapi[,3] >= d2,]
slI[j+1] <-  mapip[1,3]
mapi <- mapip 
d1 <- mapi[1,3]
}
rna.s <- slI[!is.na(slI)]
chrP <- rep(i, length(rna.s))
if(i ==1){
disn <- rna.s
chrn <- chrP 
} else{
disn <- c(rna.s, disn)
chrn <- c(chrn,chrP)
}
}
n.map <- cbind(chrn, disn)
n.map.i <- paste(n.map[,1],n.map[,2], sep=".")
map$map.i <- paste(map[,2],map[,3], sep=".")
select.mark <- as.character (map[map[,4]%in%n.map.i,][,1])
Mmat <- Mmat[,select.mark]
}
}

if(length(pRandSample)!=0){
nmark <- ncol(Mmat)
n.mark.sample <- round (nmark*pRandSample,0)
if(length(seed) !=0){
set.seed(seed)
}
id.samp <- sample(1:nmark, n.mark.sample)
Mmat <- Mmat[, id.samp]
}
if( length(impute.M) !=0){ 
impute  <- c( "knnImpute", "rfImpute", "medianImpute", "meanImpute", "modeImpute")

if (any(!(impute.M %in% impute))) 
        stop(paste("'impute method' should be one of:", paste(impute, collapse = ", ")))

if(impute.M== "meanImpute"){
Mmat <- apply(Mmat, 2, function(x){ y = mean(x, na.rm = TRUE); x[which(is.na(x))] <- y; x })
print("Note: missing values within markers are imputed with each marker genotype means")
}		

if(impute.M== "medianImpute"){
Mmat <- apply(Mmat, 2, function(x){ y = median(x, na.rm = FALSE)(x, na.rm = TRUE); x[which(is.na(x))] <- y; x })
print("Note: missing values within markers are imputed with each marker genotype median")
}
			
if(impute.M== "modeImpute"){
mode <- function(x){
temp <- table(as.vector(x))
mod <- names(temp)[temp == max(temp)]
return(mod)
}
Mmat <- apply(Mmat, 2, function(x){ y = mode(x)(x); x[which(is.na(x))] <- y; x })
print("Note: missing values within markers are imputed with each marker genotype mode")
}
}
colnames(Ymat) <- trait.n 
outobj <- list (Ymat=Ymat, Mmat=Mmat, Xmat=Xmat,Kmat=Kmat, map=map,  
minMAF=minMAF, maxMISS=maxMISS, thinINT=thinINT, pRandSample=pRandSample, seed=seed)

class(outobj) <- c("gpo", class(outobj))
return(outobj)
}
