%%  WARNING: This file was automatically generated from the associated 
%%  wav_xform.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_xform.mid file in the project MID directory. Once the 
%%  wav_xform.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavCWT, as.matrix.wavCWT, plot.wavCWT, print.wavCWT functions

\name{wavCWT}
\alias{wavCWT}
\alias{as.matrix.wavCWT}
\alias{plot.wavCWT}
\alias{print.wavCWT}
\title{Continuous wavelet transform}
\concept{transforms, waveletclass constructorwaveletsingularity detection}
\usage{wavCWT(x, scale.range=deltat(x) * c(1, length(x)), n.scale=100,
    wavelet="gaussian2", shift=5, variance=1)}
\description{The continuous wavelet transform (CWT) is a highly redundant
transformation of a real-valued or complex-valued function \eqn{f(x)}, mapping it
from the time domain to the so-called time-scale domain. Loosely,
speaking the CWT coefficients are proportional to the variability of a
function at a given time and scale.

The CWT is defined by a
complex correlation of a scaled and time-shifted mother wavelet with a
function \eqn{f(x)}.
Let \eqn{\psi(x)}{psi(x)} be a real- or complex-valued function representing a
\emph{mother} wavelet, i.e. a function which meets the standard
mathematical criteria for a wavelet
and one that can be used to generate
all other wavelets within the same family.
Let \eqn{\psi^*(\cdot)}{Conj{psi(w)}} be the complex conjugate of
\eqn{\psi(\cdot)}{psi(w)}. The CWT of
\eqn{f(x)} is defined as

\deqn{
W_f(a,b) \equiv \frac{1}{\sqrt{a}} \int_{-\infty}^\infty f(x) \psi^* \Bigl(\frac{x-b}{a}
\Bigr) \; dx,
}{W(a,b)=1/sqrt(a) Int{f(x) psi((x-b)/a)} dx}
for \eqn{(a,b) \in {\mathcal{R}}}{(a,b) in R} and \eqn{a > 0}{a > 0},
where \eqn{a} is the scale of the wavelet and \eqn{b} is the shift of the
wavelet in time. It can be shown that the above complex correlation
maintains a duality with the Fourier transform defined by the relation

\deqn{ W_f(a,b) \equiv \frac{1}{\sqrt{a}} \int_{-\infty}^\infty f(x) \psi^* \Bigl(\frac{x-b}{a}
\Bigr) \; dx \longleftrightarrow \sqrt{a} \, F(\omega) \,\Psi^*(a\omega) }{sqrt{a} * F(w) * Conj{Psi(a * w)}}

where \eqn{F(\cdot)}{F()} is the Fourier transform of \eqn{f(x)} and \eqn{\omega}{w} is the frequency in radians.
This function calculates the CWT in the Fourier domain followed by an inverse Fourier transform.}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series. The time series
may be of class class \code{rts}, \code{ts}, \code{cts}, or \code{signalSeries}, or be a numeric
vector.}
\item{n.scale}{the number of scales to evaluate over the \code{scale.range}. Default: \code{100}.}

\item{scale.range}{a two-element vector containing the range of scales over which to evaluate the CWT.
The smallest specified scale must be greater than or
equal to the sampling.interval of the time series.
Default: \code{deltat(x) * c(1, length(x))}.}

\item{shift}{numeric value representing the frequency shift to
use for the Morlet wavelet filter. Default: \code{5}.}

\item{variance}{if the wavelet filter
is of type \code{"gaussian1"} or
\code{"gaussian2"} then this parameter
represents the variance of the Gaussian PDF used to
scale the corresponding filters. Default: \code{1}.}

\item{wavelet}{a character string denoting the wavelet filter to use in
calculating the CWT. Choices are "haar", "gaussian1", "gaussian2", and "morlet",
where gaussian1 and gaussian2 represent the first and second derivatives of a
Gaussian PDF. Default: \code{"gaussian2"}.}
}

\value{
an object of class \code{wavCWT}.
}

\section{S3 METHODS}{
\describe{

\item{as.matrix}{returns the CWT coefficients as a complex matrix with rows and
columns representing times and scales, respectively.}

\item{plot}{plots the CWT. By default, the modulus of the CWT coefficients are plotted in
the time-scale plane.  The plot method also supports the following optional arguments:

\describe{
\item{phase}{Logical flag. If \code{TRUE}, the phase of the CWT is plotted. Default: \code{FALSE}.}

\item{xlab}{A character string used for the label on the x-axis. Default: \code{"time"}.}

\item{ylab}{A character string used for the label on the y-axis.
Default: \code{"log2(units)"} where \code{units} are the units of the time
series if available. If \code{units} are not available, \code{"log2(scale)"} is used.}

\item{power.stretch}{A numeric value used to scale the magnitude of the
CWT coefficients for display purposes only. If \code{power.stretch=0}, then image is
transformed with \code{log(abs(x)+1)}. Otherwise, the image is
transformed with \code{(abs(x))}\eqn{\mbox{\textasciicircum}}{^}\code{power.stretch}. Default: \code{0.5}.}

\item{type}{A character string denoting the type of plot to produce. Choices are \code{"image"} and
\code{"persp"} which plot the CWT as an image or as a meshed perspective plot, respectively. The perspective
plot resamples the data to contain a maximum of \code{grid.size} rows and columns.}

\item{zoom}{A four-element vector containing the ranges to zoom into the data in
c(xmin, xmax, ymin, ymax) format. Default: \code{NULL} (no zoom).}

\item{grid.size}{An integer specifying the maximum number of lines to use for each
dimension in creating meshed perspective plots. Default: \code{100}.}

\item{new}{A logical flag. If \code{TRUE}, a new plot is forced with the \code{frame()} command.}}}

\item{print}{prints a qualitiative summary of the CWT and its components.}

}
}

\references{
D. B. Percival and  A. T. Walden, \emph{Wavelet Methods for
Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavCWTFilters}}.}

\examples{
## calculate the CWT of the sunspots series using 
## a Mexican hat wavelet (gaussian2) 
sunspots.cwt <- wavCWT(sunspots)

## print the result 
print(sunspots.cwt)

## plot an image of the modulus of the CWT and the 
## time series 
plot(sunspots.cwt, series=TRUE)

## plot a coarse-scale wire-frame perspective of 
## the CWT 
plot(sunspots.cwt, type="persp")
}
\keyword{univar}

