%%  WARNING: This file was automatically generated from the associated 
%%  wav_filt.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_filt.mid file in the project MID directory. Once the 
%%  wav_filt.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavDaubechies, plot.wavDaubechies, print.wavDaubechies functions

\name{wavDaubechies}
\alias{wavDaubechies}
\alias{plot.wavDaubechies}
\alias{print.wavDaubechies}
\title{Daubechies wavelet and scaling filters}
\concept{class constructorfilterswavelet}
\usage{wavDaubechies(wavelet="s8", normalized=TRUE)}
\description{Ingrid Daubechies, a noted pioneer in wavelet theory, has
established a number of wavelet filter types, each with different
mathematical properties.  This function calculates the wavelet
and scaling coefficients for a given filter type.
The wavelet coefficients,
\eqn{h_k}{h(k)} for
\eqn{k=0,\ldots,L-1}{k=0,...,L-1} where
\eqn{L} is the filter length, are
related to the scaling coefficients through the quadrature mirror
filter (QMF) relation
\deqn{h_k=(-1)^{k-L} g_{L-1-k}. }{h(k)=(-1)^(k-L) g(L-1-k)}}

\arguments{

\item{normalized}{a logical value. If \code{TRUE},
the filters are normalized by dividing each filter
coefficient by the \eqn{\sqrt{2}}{sqrt(2)}
(useful for maximum overlap wavelet transforms).
If \code{FALSE}, no
normalization is used. Default: \code{TRUE}.}

\item{wavelet}{a character string denoting the filter type. Supported types include:
\describe{
\item{EXTREMAL PHASE (daublet):}{"haar", "d2", "d4", "d6", "d8", "d10", "d12", "d14", "d16", "d18", "d20"}

\item{LEAST ASYMMETRIC (symmlet):}{"s2","s4", "s6", "s8", "s10", "s12", "s14", "s16", "s18", "s20"}

\item{BEST LOCALIZED:}{"l2","l4", "l6", "l14", "l18", "l20"}

\item{COIFLET:}{"c6", "c12", "c18", "c24", "c30"}}


Default: \code{"s8"}.}
}

\value{
an object of class \code{wavDaubechies}.
}

\section{S3 METHODS}{
\describe{

\item{plot}{plot Daubechies filters.

Usage: plot(x, type="time")

\describe{
\item{x}{A \code{wavDaubechies} object.}

\item{type}{A character string denoting the type of plot to produce.
Choices are \code{"time"}, \code{"gain"}, and \code{"phase"} for an
impulse response, squared gain, and phase plot, respectively.
Default: \code{"time"}.}}}

\item{print}{print Daubechies filters.

Usage: print(x, verbose=TRUE)

\describe{
\item{x}{A \code{wavDaubechies} object.}

\item{verbose}{A logical value. If \code{TRUE}, the filter coefficients
are also printed. Default: \code{TRUE}.}}}

}
}

\details{
%
Only relevant for Daubechies filter types. Inconsistent ordering of the coefficients
in Daubechies' book was recognized and corrected by Percival (see references).
The "correct" order is given here.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.

I. Daubechies,
\emph{Orthonormal Bases of Compactly Supported Wavelets},
Communications on Pure and, Applied Mathematics, 41, 909--96.
}

\seealso{
\code{\link{wavGain}}, \code{\link{wavDWT}}, \code{\link{wavMODWT}}, \code{\link{wavMODWPT}}.}

\examples{
## obtain Daubechies least asymmetric 8-tap filter 
## set 
filters <- wavDaubechies("s8", normalize=TRUE)

## plot the impulse responses 
plot(filters, type="time")

## plot the gain function 
plot(filters, type="gain")
}
\keyword{univar}

