%%  WARNING: This file was automatically generated from the associated 
%%  wav_coef.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_coef.mid file in the project MID directory. Once the 
%%  wav_coef.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavIndex function

\name{wavIndex}
\alias{wavIndex}
\title{Boundary and interior wavelet coefficient identification}
\concept{wavelet}
\usage{wavIndex(x)}
\description{The boundary wavelet and scaling coefficients are those subject to circular
filtering operations. This function returns the range of indices which
span the interior (or nonboundary) wavelet and scaling coefficients. If approximate
zero phase filters are used in the wavelet trsnsform input then the shift
factors needed to bring the coefficients to (approximate) zero phase are also returned.}

\arguments{
\item{x}{an object of class \code{wavTransform} or \code{wavBoundary}. The transform type must be a DWT or MODWT.}

}

\value{
a list the indices locating the interior and boundary
coefficients as well as the the zero phase shift factors need for each level of the transform.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavDWT}}, \code{\link{wavMODWT}}.}

\examples{
## calculate the coefficient indices for a MODWT 
## of a simple time series 
wavIndex(wavMODWT(1:8, wavelet="s8"))
}
\keyword{univar}

