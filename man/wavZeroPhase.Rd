%%  WARNING: This file was automatically generated from the associated 
%%  wav_filt.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_filt.mid file in the project MID directory. Once the 
%%  wav_filt.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavZeroPhase function

\name{wavZeroPhase}
\alias{wavZeroPhase}
\title{Zero phase shift factors for Daubechies symmlet and coiflet filters}
\concept{waveletzero phase shifting}
\usage{wavZeroPhase(wavelet="s8", levels=1:3)}
\description{Daubechies coiflet and symmlet filters are approximate linear phase
filters. Consequently, the wavelet and scaling coefficients of the DWT
and MODWT can be circularly shifted
for approximate zero phase alignment with the original time series. This
function calculates the circular shift factors needed to bring the wavelet
and scaling coefficients to approximate zero phase.}

\arguments{

\item{levels}{an integer vector containing the decomposition levels. Default: \code{1:3}.}

\item{wavelet}{a character string denoting the filter type. See \code{wavDaubechies} for details. Default: \code{"s8"}.}
}

\value{
a list containing the shifts for each crystal of a DWTor MODWT
for the specified decomposition \code{levels}. A negative shift factor implies
an advance (circular shift to the left) of the wavelet transform crystals.
}
\details{
%
Only relevant for DWT or MODWT definitions
as given in the above reference and is valid only for Daubechies
symmlet and coiflet filters.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavDaubechies}}, \code{\link{wavDWT}}, \code{\link{wavMODWT}}, \code{\link{wavShift}}.}

\examples{
## calculate the zero phase shift factors for 
## Daubechies coiflet 12-tap filters for levels 
## 2 and 4. 
wavZeroPhase(wavelet="c12", levels=c(2,4))
}
\keyword{univar}

