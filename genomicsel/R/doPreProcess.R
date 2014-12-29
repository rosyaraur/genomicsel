doPreProcess <-
function (gpobj, impute.M = "meanImpute", transform.Y = c("boxcox", "scale", "center"),
 transform.X = c("scale", "center", "pca"),  pca.threshold = 0.8, map = NULL,  minMAF=0, 
maxMISS = 0, thinINT = NULL, pRandSample = NULL, seed=NULL){

# checking if the input object is gpo object 
if(is(gpobj) != "gpo") {
                    stop ("The input dataset must be gpo object, use function creatGPO to create one ")
} 
Ymat <- gpobj$Ymat
Mmat <- gpobj$Mmat
Xmat <- gpobj$Xmat

# check if names consistant 
if(length(transform.Y) != 0){
transformY  <- c("center", "scale",  "pca","boxcox")
if (any(!(transformY %in% transform.Y ))) 
        stop(paste("'transform method' should be one of:", paste(tranform.Y, collapse = ", ")))

if (any(transform.Y %in% c("scale", "center"))){
Ymat <- scale(Ymat, scale = TRUE)
}
if (any(transform.Y %in% c("center"))){
Ymat <- scale(Ymat, scale = FALSE)
}		
if(transform.Y == "boxcox"){
boxcoxf <- function(y){
ymin <- min(y)
if(ymin < 0){
y = (y - ymin + 1)
}
out <- boxcox(y ~ 1, lambda = seq(-5, 5, length = 100),plotit = FALSE )
df1 <- data.frame (out$x, out$y)
lambda1 <- df1[which.max(df1[,2]), 1]
y1 <- y^lambda1
return(y1)
}
Ymat <- apply(Ymat, 2, boxcoxf)
}
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

# pre-processing 
if(length(impute.M) != 0){
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
gpobj$Mmat <- Mmat
gpobj$impute.M <- impute.M
}

if(length(transform.X) !=0){
transform  <- c("center", "scale",  "pca")
if (any(!(transform %in% transform.X ))) 
        stop(paste("'transform method' should be one of:", paste(tranform.X, collapse = ", ")))
		
if (any(transform.X %in% c("scale", "center"))){
Xmat <- scale(Xmat, scale = TRUE)
}

if (any(transform.X %in% c("center"))){
Xmat <- scale(Xmat, scale = FALSE)
}

if (any(transform.X %in% c("scale", "center", "pca"))){
prout <- prcomp(Xmat, scale = TRUE, center = TRUE)
psum <- summary(prout)
var.explain <- cumsum (as.numeric (psum$importance[2,]))
pcth <- which(var.explain > pca.threshold)
p.scores <- prout$x
Xmat <- p.scores[,1:pcth]
}
gpobj$Xmat <- Xmat
gpobj$transform.X <- transform.X
}
gpobj$preprocessed <- TRUE
return(gpobj)
}
