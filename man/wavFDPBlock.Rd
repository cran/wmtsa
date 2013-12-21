%%  WARNING: This file was automatically generated from the associated 
%%  wav_fdp.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_fdp.mid file in the project MID directory. Once the 
%%  wav_fdp.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavFDPBlock function

\name{wavFDPBlock}
\alias{wavFDPBlock}
\title{Block-dependent estimation of fractionally differenced (FD) model parameters}
\concept{waveletfractionally differenced (FD) process}
\usage{wavFDPBlock(x, wavelet="s8", levels=NULL, sdf=NULL,
    boundary=NULL, edof.mode=1,
    estimator="wlse", delta.range=c(-10.0,10.0),
    position=list(from=1,by=1,units=character()), units=character(),
    title.data=character(), documentation=character(), keep.series=FALSE)}
\description{A discrete wavelet transform of the input series is used to calculate
block-dependent estimates of the
FD parameter, the variance of the FD parameter and the
innovations variance. Both a
maximum likelihood estimation (MLE) and weighted least
squares estimation (WLSE) scheme are supported.
If an MLE scheme is chosen, then the DWT
is used for its ability to de-correlate long-memory
processes. If a WLSE scheme is chosen, then the MODWT
is used for its known statistical wavelet variance properties.}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series.}
\item{boundary}{a character string representing the different methods by
which boundary wavelet coefficients and scaling coefficients
are handled in calculating the FD model parameters.
The options for this argument are dependent upon the
\code{estimator} argument.

For the \bold{MLE} case, the \code{boundary} options are:

\describe{
\item{\code{"stationary"}}{Under a stationary FD process model, boundary wavelet
and scaling coefficients are used in estimating
the FD model parameters.}

\item{\code{"nonstationary"}}{A stationary-nonstationary
FD model assumes that the governing process may fall into the
nonstationary regime and, accordingly, the boundary wavelet
coefficients and scaling coefficients are excluded in
estimating the FD model parameters.}}


For the \bold{WLSE} case, the \code{boundary} options are:

\describe{
\item{\code{"biased"}}{Boundary wavelet coefficients are included in the estimate.}

\item{\code{"unbiased"}}{Boundary wavelet coefficients are excluded in the estimate.}}


The scaling coefficients
are (always) excluded in weighted least squares estimates
of FD model parameters. Default: \code{"unbiased"}.}

\item{delta.range}{a two-element vector containing the search range for the FD parameter. Typically,
the range \eqn{[-10,10]} is suitable for all physical systems. Default: \code{c(-10 10)}.}

\item{documentation}{a character string used to describe the input
\code{data}. Default: \code{character()}.}

\item{edof.mode}{the mode by which the equivalent degrees of
freedom are calculated. This argument is
limited to 1,2, or 3 and is used only for the WLSE scheme.
See \code{wavEDOF} for details.
Default: \code{1}.}

\item{estimator}{a character string denoting the estimation method. Use \code{"wlse"} for
a weighted least squares estimate and \code{"mle"} for a maximum likelihood estimate.
Default: \code{"wlse"}.}

\item{keep.series}{a logical value. If \code{TRUE}, the original series
is preserved in the output object. Default: \code{FALSE}.}

\item{levels}{a vector containing the decomposition levels. The levels may be given
in any order but must be positive. Default: \code{1:J} where \eqn{J} is the maximum wavelet decomposition
level at which there exists at least one interior wavelet coefficient.}

\item{position}{a \code{list} containing the arguments
\code{from, by} and \code{to} which describe the position(s) of the input
\code{data}. All position arguments need not be specified as missing members
will be filled in by their default values. Default: \code{list(from=1, by=1, units=character())}.}

\item{sdf}{a vector containing a discretized approximation
of the process spectral density function (SDF). The
coefficients of this argument should correspond
exactly with the normalized Fourier frequencies
\eqn{f=[0, 1/P , 2/P, 3/P, \ldots, (M-1)/P]}{f=0, 1/P , 2/P, 3/P, ..., (M-1)/P}, where
\eqn{P=2*(M-1)} and
\eqn{M} is the number of points in the SDF
vector. For example, if the sdf vector contains five
elements, the corresponding frequencies will be
\eqn{f=[0$, $1/8$, $1/4$, $3/8$, $1/2]}{f=[0, 1/8, 1/4, 3/8, 1/2]}.
This argument is used only for the WLSE scheme when calculating EDOF mode
2 estimates. Default: \code{NULL} (EDOF mode 2 not used).}

\item{title.data}{a character string representing the name of the input
\code{data}. Default: \code{character()}.}

\item{units}{a string denoting the units of the time series. Default: \code{character()} (no units).}

\item{wavelet}{a character string denoting the filter type. See \code{wavDaubechies} for details.
Default: \code{"s8"}.}
}

\value{
an object of class \code{wavFDP}.
}
\details{
%
When \code{estimator="mle"} and
\code{boundary="stationary"},
the \code{levels} vector is forced to take on
values \eqn{{1,2,\ldots,J}}{[1,2,...,J]}
where \eqn{J} is the maximum number of levels in a full DWT.
This is done because (in this case) the scaling coefficient and all wavelet coefficients
are used to form the FD model parameter estimates.

In using the WLSE scheme it is recommended that only the unbiased
estimator be used since the confidence intervals for the
biased estimator have not been sufficiently studied.
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
\code{\link{wavEDOF}}, \code{\link{wavFDP}}, \code{\link{wavFDPTime}}, \code{\link{wavFDPBand}}, \code{\link{wavFDPSDF}}.}

\examples{
## perform a block-averaged MLE of FD parameters 
## for an FD(0.45, 1) realization over levels 1 
## through 6 using a stationary-nonstationary 
## FD model and Daubechies least asymmetric 
## 8-tap filters 
wavFDPBlock(fdp045, levels=1:6, wavelet="s8", est="mle", boundary="nonstationary")
}
\keyword{models}

