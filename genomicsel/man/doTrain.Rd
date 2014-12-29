\name{doTrain}
\alias{doTrain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doTrain(gpobj, module = "rrBLUP", method = "rrBLUP", validation = "cv", k = 10, times = 1, seed = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpobj}{
%%     ~~Describe \code{gpobj} here~~
}
  \item{module}{
%%     ~~Describe \code{module} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{validation}{
%%     ~~Describe \code{validation} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (gpobj, module = "rrBLUP", method = "rrBLUP", validation = "cv", 
    k = 10, times = 1, seed = NULL, ...) 
{
    if (is(gpobj) != "gpo") {
        stop("The input dataset must be gpo object, use function creatGPO to create one ")
    }
    RMSE_in <- rep(NA, k)
    RMSE_out <- rep(NA, k)
    Rsqr_in <- rep(NA, k)
    Rsqr_out <- rep(NA, k)
    train.fold <- doFoldCreate(gpobj, k, seed = seed)
    for (i in 1:k) {
        test_temp <- train.fold[[i]]
        allv <- 1:nrow(gpobj$Ymat)
        train_temp <- setdiff(allv, test_temp)
        yt <- matrix(gpobj$Ymat[train_temp, ], ncol = 1)
        rownames(yt) <- rownames(gpobj$Ymat)[train_temp]
        colnames(yt) <- colnames(gpobj$Ymat)
        Mt = gpobj$Mmat[train_temp, ]
        if (length(gpobj$Xmat) != 0) {
            Xt = gpobj$Xmat[train_temp, ]
        }
        else {
            Xt = NULL
        }
        if (length(gpobj$kmat) != 0) {
            Kt = gpobj$kmat[train_temp, train_temp]
        }
        else {
            Kt = NULL
        }
        tmp_gpo <- creatGPO(Ymat = yt, Mmat = Mt, Xmat = Xt, 
            Kmat = Kt)
        yts <- matrix(gpobj$Ymat[test_temp, ], ncol = 1)
        rownames(yts) <- rownames(gpobj$Ymat)[test_temp]
        colnames(yts) <- colnames(gpobj$Ymat)
        Mts = gpobj$Mmat[test_temp, ]
        if (length(gpobj$Xmat) != 0) {
            Xts = gpobj$Xmat[test_temp, ]
        }
        else {
            Xts <- NULL
        }
        if (length(gpobj$kmat) != 0) {
            Kts = gpobj$kmat[test_temp, test_temp]
        }
        else {
            Kts <- NULL
        }
        test_gpo <- creatGPO(Ymat = yts, Mmat = Mts, Xmat = Xts, 
            Kmat = Kts)
        ans <- fit.model(gpobj = tmp_gpo, module = module, method = method, 
            ...)
        if (c("ksvm") \%in\% class(ans)) {
            RMSE_in[i] <- ans@param$RMSE
            Rsqr_in[i] <- ans@param$r.squared
        }
        else {
            Rsqr_in[i] <- ans$r.squared
            RMSE_in[i] <- ans$RMSE
        }
        ans.out <- rmse.out(ans, newgobj = test_gpo, gpobj = NULL, 
            ...)
        if (c("ksvm") \%in\% class(ans.out)) {
            RMSE_out[i] <- ans.out@param$RMSE
            Rsqr_out[i] <- ans.out@param$r.squared
        }
        else {
            Rsqr_out[i] <- ans.out$r.squared
            RMSE_out[i] <- ans.out$RMSE
        }
    }
    outdf <- data.frame(k = 1:k, RMSE_in = RMSE_in, Rsqr_in = Rsqr_in, 
        RMSE_out = RMSE_out, Rsqr_out = Rsqr_out)
    class(outdf) <- c("doTrain")
    return(outdf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
