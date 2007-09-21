%%  WARNING: This file was automatically generated from the associated 
%%  wav_xform.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_xform.mid file in the project MID directory. Once the 
%%  wav_xform.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavBestBasis function

\name{wavBestBasis}
\alias{wavBestBasis}
\title{DWPT Best basis selection}
\concept{transforms, waveletbest basis}
\usage{wavBestBasis(costs)}
\description{The discrete wavelet packet transform (DWPT) contains a multitude of
disjoint dyadic decompositions representing an ensemble of different
bases. Best basis selection is an attempt to isloate one such basis
in an \emph{optimal} way.}

\arguments{
\item{costs}{a numeric vector containing the costs for each
crystal in a DWPT in
\eqn{{\cal{C}}(W_{0,0}), {\cal{C}}(W_{1,0}), {\cal{C}}(W_{1,1}), {\cal{C}}(W_{2,0}), \ldots, {\cal{C}}(W_{J,2^J-1})}{C(W(0,0)), C(W(1,0)), C(W(1,1)), C(W(2,0)), ..., C(W(J,2^J-1))}
order where \eqn{{\cal{C}}(\cdot)}{C()} is the additive cost functional and
\eqn{W_{j,n}}{W(j,n)} is the DWPT crystal at level \eqn{j} and oscillation (local node)
index \eqn{n} for \eqn{j=1,\ldots,J}{j=1,...,J}.}

}
\references{
Ronald R. Coifman and Mladen Victor Wickerhauser,
``Entropy-Based Algorithms for Best Basis Selection",
\emph{IEEE Transactions on Information Theory},
\bold{38}(2), pp. 713--718, 1992.

D. B. Percival and  A. T. Walden, \emph{Wavelet Methods for
Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavDWPT}}.}

\examples{
## calculate the DWPT of the difference of the 
## atomic clock series 
W <- wavDWPT(diff(atomclock), n.level=6)

## define an entropy cost functional 
"entropy" <- function(x){
    iz <- which(x==0)
    z <- -x^2 * log(x^2)
    if (length(iz))
       z[iz] <- 0
    sum(z)
}

## create the cost vector 
C <- unlist(lapply(W$data, entropy))

## calculate the bets basis 
z <- wavBestBasis(C)

## print the crystals of the best basis 
paste("W(", z$level, ",", z$osc, ")", sep="")
}
\keyword{univar}

