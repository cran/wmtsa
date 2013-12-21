%%  WARNING: This file was automatically generated from the associated 
%%  wav_fdp.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_fdp.mid file in the project MID directory. Once the 
%%  wav_fdp.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavFDPBand function

\name{wavFDPBand}
\alias{wavFDPBand}
\title{Mid-octave spectral density function (SDF) estimation}
\concept{waveletfractionally differenced (FD) process}
\usage{wavFDPBand(delta=1/4, method="bandpass", scaling=TRUE,
    levels=1:5, n.sample=n.sample <- 2^(max(levels)+1))}
\description{The wavelet and scaling filters used for wavelet decompositions
are nominally associated with approximate bandpass filters.
Specifically, at decomposition level \eqn{j},
the wavelet transform
coefficients correspond approximately to the normalized
frequency range of \eqn{[ 1/2^{j+1}, 1 /2^j ]}{[ 1/2^(j+1), 1 /2^j ]}. The square of the
wavelet coefficients are used to form the so-called wavelet
variance (or wavelet spectrum) which is seen as a regularization
of the SDF. Under an assumed FD process, this
function estimates the mid-octave SDF values. The estimates are
calculated assuming that the wavelet transform filters form
perfect (rectangular) passbands. Decomposition levels
\eqn{1} and \eqn{2}
are calculated using a second order Taylor series expansion
about the mid-octave frequencies while, for levels greater
than \eqn{2}, a small angle approximation
(\eqn{sin(\pi f) \approx \pi f}{sin(pi * f) ~ pi * f})
is used to develop a closed form solution which is a function of
FD model parameters as well as the mid-octave frequencies.}

\arguments{

\item{delta}{the fractional difference parameter. If the
scaling band estimates are desired (prompted by
setting \eqn{\mbox{\tt n.sample} > 0}{n.sample > 0}),
then \code{delta} must be less
than 0.5 since the formulae for calculating the
scaling band estimates implicitly assume
stationarity. Default: \code{0.4}.}

\item{levels}{a vector containing the decomposition levels.
If \eqn{\mbox{\tt n.sample} \le 0}{n.sample <= 0}, then
the levels may be given in any order and levels may be skipped.
If, however, \code{n.sample} > \eqn{0},
then \code{levels} must contain the values
\eqn{1, 2, 3,\ldots, J}{1, 2, 3, ..., J}
where \eqn{J} is the maximum wavelet
transform decomposition level. Default: \code{1:5}.}

\item{method}{a character string denoting the method to be used for estimating the average
spectral density values at the center frequency (on a log scale)
of each DWT octave. The choices are

\describe{
\item{\code{"integration"}}{Numerical integration of the standard FDP spectral density function.}

\item{\code{"bandpass"}}{A small angle approximation to the standard FDP spectral density functions for
decomposition levels \eqn{j \ge 3}{j >= 3} in combination with a Taylor series approximation for levels
\eqn{j=1,2}{j=1,2}.}}

Default: \code{"bandpass"}.}

\item{n.sample}{the number of samples in the time series.
Although no time series is actually passed to
the \code{wavFDPBand} function,
the \code{n.sample} argument is used
in estimating the mid-octave SDF value over the band of frequencies
which are nominally associated with the scaling
filter in a wavelet transform.
If \eqn{\mbox{\tt n.sample} > 0}{n.sample > 0},
this function will append the estimate of the
average SDF value over the scaling band to the
wavelet octave estimates.
If \eqn{\mbox{\tt n.sample} \le 0}{n.sample <= 0},
only the wavelet octave estimates are returned. Default: \code{1024}.}

\item{scaling}{a logical flag.  If
\code{TRUE}, the mid-octave value of the FDP SDF
octave corresponding to the scaling coefficients is also returned.
Default: \code{TRUE}.}
}

\value{
a vector containing the mid-octave SDF estimates for an FD process.
}
\details{
%
Estimates are made for the scaling filter band based upon
an implicit assumption that the FD process is stationary
(\eqn{\delta < 1/2}{delta < 0.5}).
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000, 343--54.
}

\seealso{
\code{\link{wavFDPBlock}}, \code{\link{wavFDPTime}}, \code{\link{wavVar}}, \code{\link{wavFDPSDF}}.}

\examples{
## calculate the mid-octave SDF values for an FD 
## process over various wavelet bands 
wavFDPBand(levels=c(1, 3, 5:7), delta=0.45, scaling=FALSE)
}
\keyword{models}

