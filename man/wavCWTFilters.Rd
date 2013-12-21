%%  WARNING: This file was automatically generated from the associated 
%%  wav_filt.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_filt.mid file in the project MID directory. Once the 
%%  wav_filt.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavCWTFilters function

\name{wavCWTFilters}
\alias{wavCWTFilters}
\title{Frequency response of continuous wavelet transform filters}
\concept{filters, waveletwaveletsingularity detection}
\usage{wavCWTFilters(wavelet="Gaussian2", frequency=seq(0, 2 * pi, length=1000),
    shift=3, variance=1, times=NULL)}
\description{Returns the frequency response of a continuous
wavelet filter. The choices for filters are limited to
Haar, Gaussian, and Morlet families.}

\arguments{

\item{frequency}{a numeric vector denoting the frequencies (in rad/sec) over which the
frequency response function for the specified wavelet should be evaluated.
Default: \code{seq(0, 2 * pi, length=1000)}}

\item{shift}{the frequency shift \eqn{\omega_0}{w0} of the Morlet wavelet. Default: \code{3}.}

\item{times}{a numeric vector of values corresponding to times at which the
specified filter should be evaluated. If not \code{NULL}, the impulse response of
the specified filter is returned, otherwise the freequency response is returned.
Default: \code{NULL}.}

\item{variance}{the variance of a Gaussian PDF. Used only for the (derivatives of) Guassian filters.
Default: \code{1}.}

\item{wavelet}{a character string denoting the wavelet filter.
Choices are "haar", "gaussian1", "gaussian2", and "morlet",
where \code{gaussian1} and \code{gaussian2} represent the first and second derivatives of a
Gaussian PDF. Default: \code{"gaussian2"}.}
}

\value{
the frequency response corresponding to the input frequencies.
}
\references{
D. B. Percival and  A. T. Walden, \emph{Wavelet Methods for
Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavCWT}}.}

\examples{
## set the plot grid 
old.plt <- splitplot(2,2,1)

## create a frequency vector 
w <- seq(-10,10,length=1000)

## calculate the frequency responses for various 
## mother wavelets and plot the results 

wavelets <- c("haar", "gaussian1", "gaussian2", "morlet")

for (i in seq(wavelets)){
    if (i > 1) splitplot(2,2,i)
    filt <- wavelets[i]
    plot(w, abs(wavCWTFilters(wavelet=filt, frequency=w)),
       ylab="|Psi(w)|", xlab="frequency", type="l")
title(filt)
}

par(old.plt)
}
\keyword{univar}

