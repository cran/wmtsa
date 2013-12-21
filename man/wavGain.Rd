%%  WARNING: This file was automatically generated from the associated 
%%  wav_filt.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_filt.mid file in the project MID directory. Once the 
%%  wav_filt.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavGain, plot.wavGain, print.wavGain functions

\name{wavGain}
\alias{wavGain}
\alias{plot.wavGain}
\alias{print.wavGain}
\title{The gain functions for Daubechies wavelet and scaling filters}
\concept{waveletfilters}
\usage{wavGain(wavelet="s8", n.levels=5, n.fft=1024, normalize=TRUE)}
\description{Given \eqn{\{g\}}{g} and \eqn{\{h\}}{h}
are the impulse responses for the scaling
and wavelet filters, respectively, and
\eqn{G_1(f)}{G(1,1f)} and
\eqn{H_1(f)}{H(1,f)} are their
corresponding gain functions, then the gain functions for decomposition
level \eqn{j > 1}{j > 1} are calculated using the
recursive algorithm:

\deqn{ H_j(f)=H_1(2^{j-1}f) G_{j-1}(f),\] \[ G_j(f)=G_1(2^{j-1}f) G_{j-1}(f). }{H(j,f)=H(1, 2^(j-1)*f) * G(j-1, f), G(j,f)=G(1, 2^(j-1)*f) * G(j-1,f).}}

\arguments{

\item{n.fft}{the number of Fourier coefficients to use in approximating the
gain functions. Default: \code{1024}.}

\item{n.levels}{the number of decomposition levels. Default: \code{5}.}

\item{normalize}{a boolean value. If \code{TRUE},
the filters are normalized by dividing each filter
coefficient by the \eqn{\sqrt{2}}{sqrt(2)}
(used for maximal overlap wavelet transforms).
If \code{FALSE}, no normalization is used.
Default: \code{TRUE}.}

\item{wavelet}{a character string denoting the filter type.
See \code{wavDaubechies} for details. Default: \code{"s8"}.}
}

\value{
an object of class \code{wavGain}.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.

I. Daubechies,
\emph{Orthonormal Bases of Compactly Supported Wavelets},
Communications on Pure and, Applied Mathematics, 41, 909--96.
}

\seealso{
\code{\link{wavDaubechies}}.}

\examples{
## approximate the gain functions for the 
## normalized Daubechies least asymmetric 
## 20-tap filters for levels 1,...,5 using a 
## 1024 Fourier frequencies 
result <- wavGain(wavelet="s20", n.levels=5,
    norm=TRUE)

## plot the results 
plot(result)
}
\keyword{univar}

