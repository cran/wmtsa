%%  WARNING: This file was automatically generated from the associated 
%%  wav_fdp.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_fdp.mid file in the project MID directory. Once the 
%%  wav_fdp.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavFDPTime function

\name{wavFDPTime}
\alias{wavFDPTime}
\title{Instantaneous estimation of fractionally differenced model parameters}
\concept{waveletfractionally differenced (FD) process}
\usage{wavFDPTime(x, wavelet="s8", levels=NULL,
    biased=FALSE, estimator="mle",
    dof.order=0, delta.range=c(-10.0,10.0),
    position=list(from=1,by=1,units=character()), units=character(),
    title.data=character(), documentation=character(), keep.series=FALSE)}
\description{The MODWT is used to calculate instantaneous estimates of the FD parameter,
the variance of the FD parameter and the innovations variance.
The user can select between maximum likelihood and least squares
estimators. Localized estimates may also be formed by using multiple
chi-squared degrees of freedom in estimating the FD model parameters.}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series.}
\item{biased}{a logical flag used to choose between denoting biased or unbiased
estimates. Biased estimates are those which
use all available levels in calculating
the FD model parameters. Unbiased estimates
are calculated with only those wavelet
coefficients not subject to circular filter
operations, i.e. only the interior wavelet
coefficients are used in calculating unbiased
estimates. Default: \code{TRUE}.}

\item{delta.range}{a two-element vector containing the search range for the FD parameter.
Typically, the range \eqn{[-10,10]} is suitable for all physical systems. Default: \code{c(-10, 10)}.}

\item{documentation}{a character string used to describe the input
\code{x}. Default: \code{character()}.}

\item{dof.order}{the degree of freedom (DOF) order. The number of
chi-square DOFs used in estimating the FD
parameters is equal to \eqn{2 \times \mbox{\tt dof.order} + 1}{2 * dof.order + 1}
where necessarily \eqn{\mbox{\tt dof.order} > 0}{dof.order > 0}. As the order
increases, the estimates will become smoother but less localized in time. Default: \code{0}.}

\item{estimator}{a character string denoting the estimation method.
Use "lse" for least squares estimates and "mle" for maximum likelihood estimates.
Default: \code{"lse"}.}

\item{keep.series}{a logical value. If \code{TRUE}, the original series
is preserved in the output object. Default: \code{FALSE}.}

\item{levels}{a vector containing the decomposition levels. The levels may be given
in any order but must be positive. Default: \code{1:J} where \eqn{J} is the maximum wavelet decomposition
level at which there exists at least one interior wavelet coefficient.}

\item{position}{a \code{list} containing the arguments
\code{from, by} and \code{to} which describe the position(s) of the input
\code{x}. All position arguments need not be specified as missing members
will be filled in by their default values. Default: \code{list(from=1, by=1, units=character())}.}

\item{title.data}{a character string representing the name of the input
\code{x}. Default: \code{character()}.}

\item{units}{a string denoting the units of the time series. Default: \code{character()} (no units).}

\item{wavelet}{a character string denoting the filter type.
See \code{wavDaubechies} for details.
Default: \code{"s8"}.}
}

\value{
an object of class  \code{wavFDP}.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000, 340--92.

W. Constantine, D. B. Percival and  P. G. Reinhall,
\emph{Inertial Range Determination for Aerothermal Turbulence Using
Fractionally Differenced Processes and Wavelets},
Physical Review E, 2001,  64(036301), 12 pages.
}

\seealso{
\code{\link{wavFDP}}, \code{\link{wavFDPBlock}}, \code{\link{wavFDPBand}}, \code{\link{wavFDPSDF}}.}

\examples{
## perform a unbiased instantaneous LSE of FD 
## parameters for an FD(0.45, 1) realization 
## over levels 1 through 6 using Daubechies 
## least asymmetric 8-tap filters. Use a zeroth 
## order DOF (equivalent to 1 chi-square DOF) 
z <- wavFDPTime(fdp045, levels=1:6, wavelet="s8", est="lse", biased=FALSE)

## display the results 
print(z)

## plot the results 
plot(z)

## plot the results with the confidence intervals 
## centered about the mean (known) value of the 
## the FD parameter 
plot(z, mean.delta=0.45)
}
\keyword{models}

