%%  WARNING: This file was automatically generated from the associated 
%%  wav_coef.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_coef.mid file in the project MID directory. Once the 
%%  wav_coef.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavBoundary, [.wavBoundary, plot.wavBoundary, print.wavBoundary, print.summary.wavBoundary, summary.wavBoundary functions

\name{wavBoundary}
\alias{wavBoundary}
\alias{[.wavBoundary}
\alias{plot.wavBoundary}
\alias{print.wavBoundary}
\alias{print.summary.wavBoundary}
\alias{summary.wavBoundary}
\title{Wavelet transform boundary coefficent identification}
\concept{wavelet}
\usage{wavBoundary(x)}
\description{A wavelet transform boundary coefficient is one subject to circular filter
operations (or other boundary treatments). Conversely, the interior transform
coefficients are those that are not affected by the imposed boundary treatment.
The \code{wavBoundary}
function separates the boundary coefficients from the interior
wavelet transform coefficients.}

\arguments{
\item{x}{a DWT or MODWT transform object with class \code{wavTransform}.}

}

\value{
an object of class \code{wavBoundary}.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.

I. Daubechies,
\emph{Orthonormal Bases of Compactly Supported Wavelets},
Communications on Pure and, Applied Mathematics, 41, 909--96.
}

\seealso{
\code{\link{wavIndex}}, \code{\link{wavDWT}}, \code{\link{wavMODWT}}, \code{\link{wavShift}}.}

\examples{
## calculate the MODWT of the sunspots series 
W <- wavMODWT(sunspots)

## identify the boundary coefficients 
z <- wavBoundary(W)

## plot the results 
plot(wavShift(z))

## obtain a summary 
summary(z)
}
\keyword{univar}

