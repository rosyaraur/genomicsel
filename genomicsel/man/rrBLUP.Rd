\name{rrBLUP}
\alias{rrBLUP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rrBLUP(gpobj, method = "REML", bounds = c(1e-09, 1e+09))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gpobj}{
%%     ~~Describe \code{gpobj} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{bounds}{
%%     ~~Describe \code{bounds} here~~
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
function (gpobj, method = "REML", bounds = c(1e-09, 1e+09)) 
{
    require(rrBLUP)
    yt <- gpobj$Ymat
    Mt = gpobj$Mmat
    Xt = gpobj$Xmat
    Kt = gpobj$kmat
    ymean <- mean(as.vector(yt), na.rm = TRUE)
    outobj <- mixed.solve(y = yt, Z = Mt, K = Kt, X = Xt, bounds = bounds, 
        method = method, SE = FALSE, return.Hinv = FALSE)
    mixed.u <- outobj$u
    mixed.beta <- outobj$beta
    outobj$ymean <- ymean
    if (length(Xt) != 0) {
        yhat <- ymean + Xt \%*\% mixed.beta + Mt \%*\% mixed.u
    }
    else {
        yhat <- ymean + Mt \%*\% mixed.u
    }
    outobj$yhat <- yhat
    residual <- yt - yhat
    outobj$residual <- residual
    RMSE = sqrt(1/length(yt) * (sum((residual)^2)))
    outobj$RMSE <- RMSE
    outobj$r.squared <- summary(lm(yt ~ yhat))$r.squared
    class(outobj) <- c("rrBLUP", class(outobj))
    return(outobj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
