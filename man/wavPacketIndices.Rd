%%  WARNING: This file was automatically generated from the associated 
%%  wav_xform.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_xform.mid file in the project MID directory. Once the 
%%  wav_xform.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavPacketIndices function

\name{wavPacketIndices}
\alias{wavPacketIndices}
\title{Wavelet packet node indices}
\concept{wavelet}
\usage{wavPacketIndices(x, check.basis=TRUE)}
\description{Converts flattened wavelet packet node indices to corresponding \emph{level} and \emph{oscillation}
indices.}

\arguments{
\item{x}{a vector of flattened wavelet packet node indices.}
\item{check.basis}{a logical value. If \code{TRUE}, the set of
specified indices is checked to ensure that the corresponding
wavelet packet nodes form a legitimate basis by ensuring that
(i)	the union of all frequency ranges corresponding to the packet crystals span
the normalized frequencies [0,1/2] and (ii) the normalized frequency ranges
for all nodes do not overlap. Default: \code{TRUE}.}
}

\value{
a list of \code{flat}, \code{level}, and \code{osc} vectors containing the flattened, decomposition level,
and oscillation indices, respectively, corresponding to the input.
}
\details{
%
In general, wavelet packet crystals can be arranged in the so-called \emph{natural order}
ala \eqn{W_{0,0} , W_{1,0}, W_{1,1}, W_{2,0}, W_{2,1}, W_{2,2}, W_{2,3}, \ldots , W_{J,0}, \ldots, W_{J, N_J}}{W(0,0) , W(1,0), W(1,1), W(2,0), W(2,1), W(2,2), W(2,3), ... , W(J,0), ..., W(J, NJ)}
where \eqn{J} is the number of decomposition levels and \eqn{N_J=2^J-1}{NJ}.
By definition, \eqn{W_{0,0}}{W(0,0)} is the original time series.
A given crystal is identified in the \eqn{W_{j,n}}{W(j,n)} form above or by a flattened index.
For example, the DWPT crystal in level 2 at oscillation index 1 is identified either by {j,n}={2,1} or
by its flattened index 4 (with zero based indexing, 4 represents the fifth crystal of the wavelet packet
transform in natural order). This function converts such flattened wavelet packet indices to
the \eqn{W_{j,n}}{W(j,n)} form.
}
\seealso{
\code{\link{wavDWPT}}.}

\examples{
## specify a basis formed by the flattened indices 
## corresponding to the wavelet packet nodes 
## W(2,0:1) and W(3,4:7), but submit them in 
## arbitrary order 
wavPacketIndices(c(14,11,12,13,4,3), check.basis=TRUE)
}
\keyword{univar}

