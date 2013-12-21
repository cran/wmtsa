%%  WARNING: This file was automatically generated from the associated 
%%  wav_var.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_var.mid file in the project MID directory. Once the 
%%  wav_var.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavVarConfidence function

\name{wavVarConfidence}
\alias{wavVarConfidence}
\title{Wavelet variance confidence intervals}
\concept{wavelet}
\usage{wavVarConfidence(wvar, edof, probability=0.95)}
\description{This function calculates wavelet variance confidence intervals
for the unbiased and block averaged discrete wavelet variance estimates.
Given \eqn{\hat\nu_X^2(\tau_j)}{var{tau(j)}} are the
time independent unbiased wavelet variance estimates at scales
\eqn{\tau_j\equiv 2^{j-1}}{tau(j)} where
\eqn{j} are the
decomposition levels, the approximate
\eqn{100(1-2p)}\% confidence interval is given by
\deqn{ \biggl[{ n\hat\nu_X^2(\tau_j) \over Q_n(1-p) } , { n\hat\nu_X^2(\tau_j) \over Q_n(p) } \biggr] }{[n * var{tau(j)} / Q(1-p), n * var{tau(j)} / Q(p)]}
where \eqn{Q_n(p)}{Q(p)} is the
\eqn{p \times 100}{p x 100} percentage
point for a chi-squared distribution with n degrees of freedom  distribution.}

\arguments{
\item{wvar}{a vector containing the block-averaged unbiased wavelet variance estimates.}

\item{edof}{a vector containing the equivalent degrees of freedom estimates. See \code{\link{wavEDOF}} for details.}
\item{probability}{the probability desired for the confidence
intervals. Supported probabilities are 0.005, .025, .05, .95, .975, and .995. Default: \code{0.95}.}
}

\value{
a list of the low and high confidence interval limits for
levels \eqn{1,\ldots, J}{1,..., J}.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavVar}}, \code{\link{wavEDOF}}.}

\examples{
## first calculate the EDOF for the ocean series 
edof <- wavEDOF(ocean)

## calculate the 95% confidence intervals for EDOF 
## mode 1 
wavVarConfidence(edof$variance.unbiased, edof$EDOF1)
}
\keyword{distribution}

