%%  WARNING: This file was automatically generated from the associated 
%%  wav_gnrc.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_gnrc.mid file in the project MID directory. Once the 
%%  wav_gnrc.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the reconstruct function

\name{reconstruct}
\alias{reconstruct}
\title{Reconstruction (inverse transform) of various wavelet transforms}
\concept{wavelet}
\usage{reconstruct(x, ...)}
\description{Inverts a discrete wavelet transform, mapping the data back into
the time domain.}

\arguments{
\item{x}{an object of class \code{wavTransform} or \code{wavMRD}.}
\item{...}{optional arguments passed directly to the reconstruction function.}
}

\value{
a numeric vector containing the result.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.

I. Daubechies,
\emph{Orthonormal Bases of Compactly Supported Wavelets},
Communications on Pure and, Applied Mathematics, 41, 909--96.
}

\seealso{
\code{\link{wavDWT}}, \code{\link{wavMODWT}}, \code{\link{wavMRD}}.}

\examples{
## create a random series, calculate a forward 
## transform, then invert 
x       <- rnorm(1024)
x.dwt   <- reconstruct(wavDWT(x))
x.modwt <- reconstruct(wavMODWT(x))
all(c(vecnorm(x.dwt-x), vecnorm(x.modwt-x)) < .Machine$single.eps)
}
\keyword{methods}

