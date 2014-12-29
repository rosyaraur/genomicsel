doDataPartition <-
function (y, times = 1, partitions=3, p = c(0.6,0.3), seed = NULL){
   if(any(p > 1) | any(p <0) ){
      stop("any value in p should be between 0 and 1")
	  }
   if(sum(p) > 1){
      stop("the sum of the p should be less than equal to 1")
	  }
	if(length(p) + 1 != partitions){
      stop( "number of p value provided must be equal to partitions-1")
}	  
 pI <- 1-sum(p)
 if(pI ==0){
    stop(" the last undefined partition has p=0")
	}
	
  if(any(class(y) %in% "gpo")){
	n = nrow(y$Ymat)	
	   }
   if(is.vector(y)){
		n = length(y)
	}
   if(is.matrix(y)){
		n = nrow(y)
	}
n1 = round(n*p,0)
npL <- n-sum(n1)
npL1 <- c(n1, npL)

 
     if(length(seed) !=0){
	 set.seed(seed)
	 }
	 all.list <- list()
	 
	 for ( j in 1:times){
	 if(length(seed) !=0){
	 	if(length(seed) < times){
	 		stop("Number of seeds should be equal to number of times")
	  }
	  	 set.seed(seed[j])
	  }
	  partions.id <- list()
	  partdata <- list() 
	  for ( i in 1:length(npL1)){
	  if(i ==1){
	  spool <- 1:n
	  smp <- sample(spool, size = npL1[i], replace = FALSE, prob = NULL)
	  partions.id[[i]] <- smp
	  } else {
	  spool <- setdiff(spool, smp)
	  smp <- sample(spool, size = npL1[i], replace = FALSE, prob = NULL)
	  partions.id[[i]] <- smp
	 }
	 if(any(class(y) %in% "gpo")){
	 tS <- partions.id[[i]]
	 Ymat1 <- matrix(y$Ymat[tS,], ncol=1)
	 rownames(Ymat1) <- unlist(rownames(y$Ymat)[tS])
	 Mmat1 <- y$Mmat
	 Mmat1 <- Mmat1[tS,]
	 if(length(y$Xmat) != 0){
	 Xmat1 <- y$Xmat[tS,]
	 rownames(Xmat1) <- rownames(y$Xmat)[tS]
	 	 } else{
	 Xmat1 <- y$Xmat
	 }
	  if(length(y$Kmat) != 0){
	 Kmat1 <- y$Kmat[tS,tS]
	 nms <- rownames(y$Kmat)[tS]
	 rownames(Kmat1) <- nms
	 colnames(Kmat1) <- nms
	 } else {
	 Kmat1 <- y$Kmat
	 }
	 map1 <- y$map
	 minMAF <- y$minMAF
	 minMISS <- y$minMISS
	 thinINT <- y$thinINT
	 pRandSample <- y$pRandSample
	 seed <- y$seed
	 outobj <- list (Ymat=Ymat1, Mmat=Mmat1, Xmat=Xmat1,Kmat=Kmat1, map=map1,  minMAF=minMAF, minMISS=minMISS, thinINT=thinINT, pRandSample=pRandSample, seed=seed)
     class(outobj) <- c("gpo", class(outobj))
	 partdata[[i]] <- outobj
	 }
	 }
	 if(any(class(y) %in% "gpo")){
	 all.list[[j]] <- partdata
	 } else{
	 all.list[[j]] <- partions.id
    }
}	
	class(all.list) <- c("partobj", class(all.list))
	return(all.list)
}
