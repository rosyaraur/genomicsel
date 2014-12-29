\name{sepPartition}
\alias{sepPartition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
internal function - Seperate partition 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sepPartition(partobj, rep = 1, part = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{partobj}{
%%     ~~Describe \code{partobj} here~~
}
  \item{rep}{
%%     ~~Describe \code{rep} here~~
}
  \item{part}{
%%     ~~Describe \code{part} here~~
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
function (partobj, rep = 1, part = 1) 
{
    if (c("partobj") \%in\% class(partobj)) {
        data <- partobj[[rep]][[part]]
        return(data)
    }
    else {
        stop("The input dataset must be partobj, use function  doDataPartition function \n\t\t\t\t\tto build partobj")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
