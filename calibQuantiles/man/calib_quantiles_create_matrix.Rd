\name{calib_quantiles_create_matrix}
\alias{calib_quantiles_create_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calib_quantiles_create_matrix(X_q, N, totals_q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X_q}{
%%     ~~Describe \code{X_q} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
  \item{totals_q}{
%%     ~~Describe \code{totals_q} here~~
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
function (X_q, N, totals_q) 
{
    stopifnot(`dim of X_q is not equal to length of totals_q` = NCOL(X_q) == 
        length(totals_q))
    A <- list()
    for (k in 1:NCOL(X_q)) {
        totals_q_k <- totals_q[[k]]
        A_q <- matrix(0, nrow = NROW(X_q), ncol = NROW(totals_q[[k]]))
        x_sorted <- sort(X_q[, k])
        for (i in 1:NROW(totals_q_k)) {
            poz <- which(x_sorted <= totals_q_k[i])
            n_poz <- NROW(poz)
            L <- x_sorted[n_poz]
            U <- x_sorted[n_poz + 1]
            B <- (totals_q_k[i] - L)/(U - L)
            A_q[X_q[, k] < L, i] <- 1/N
            A_q[X_q[, k] == U, i] <- B/N
        }
        A[[k]] <- A_q
    }
    A <- Reduce("cbind", A)
    return(A)
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
