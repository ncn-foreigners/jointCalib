\name{calib_quantiles}
\alias{calib_quantiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calib_quantiles(X, X_q, d, N, totals, totals_q, backend = c("sampling", "laeken"), method = c("raking", "linear", "logit"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{X_q}{
%%     ~~Describe \code{X_q} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{totals}{
%%     ~~Describe \code{totals} here~~
}
  \item{totals_q}{
%%     ~~Describe \code{totals_q} here~~
}
  \item{backend}{
%%     ~~Describe \code{backend} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (X, X_q, d, N, totals, totals_q, backend = c("sampling", 
    "laeken"), method = c("raking", "linear", "logit"), ...) 
{
    totals_q_vec <- unlist(totals_q)
    quantiles <- names(totals_q_vec)
    if (all(grepl("\%", quantiles))) {
        quantiles <- as.numeric(gsub("\%", "", quantiles))/100
    }
    T_mat <- c(N, quantiles, totals)
    A <- calib_quantiles_create_matrix(X_q, N, totals_q)
    X <- cbind(1, A, X)
    if (backend == "sampling") {
        w_res <- sampling::calib(Xs = X, d = d, total = T_mat, 
            method = method, ...)
    }
    if (backend == "laeken") {
        w_res <- laeken::calibWeights(X = X, d = d, totals = T_mat, 
            method = method, ...)
    }
    w <- w_res * d
    return(list(w = w, Xs = X, totals = T_mat, diff = colSums(X * 
        w) - T_mat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
