%%  WARNING: This file was automatically generated from the associated 
%%  wav_xform.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_xform.mid file in the project MID directory. Once the 
%%  wav_xform.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavPacketBasis function

\name{wavPacketBasis}
\alias{wavPacketBasis}
\title{Extract wavelet packet basis from a DWPT}
\concept{waveletwavelet packet basistransforms, discrete wavelet packet transform}
\usage{wavPacketBasis(x, indices=0)}
\description{Returns the DWPT crystals (in a list) corresponding to the
basis specified by the indices vector. The indices
are mapped as follows:

\describe{
\item{0}{original series}

\item{1:2}{\eqn{\{W_{1,0}, W_{1,1}\}}{W(1,0), W(1,1)}, i.e., all level 1 crystals}

\item{3:6}{\eqn{\{W_{2,0},\ldots, W_{2,3}\}}{W(2,0),...,W(2,3)}, i.e., all level 2 crystals}}


and so on. If the indices do not form a basis, an error is issued.}

\arguments{
\item{x}{an object of class \code{wavTransform} associated with the output
of the \code{wavDWPT} function.}
\item{indices}{an integer vector. Each integer denotes a particular crystal of the DWPT
to extract. The set of crystals shoudl form a basis, i.e., the collective frequency ranges
associated with the set of crystals should span normalized frequencies [0, 1/2]. The indices for each
DWPT level and the corresponding (ideal) normalized frequency ranges are listed in the table below:

\describe{
\item{0}{Frequency range: [0, 1/2], associated with crystal \eqn{W_{0,0}}{W(0,0)} (the original series).}

\item{1,2}{Frequency range: [0,1/4],[1/4, 1/2], associated with crystals \eqn{W_{1,0}$, $W_{1,1}}{W(1,0), W(1,1)}, respectively.}

\item{3,4,5,6}{Frequency range: [0,1/8],[1/8, 1/4],[1/4,3/8],[3/8, 1/2], associated with crystals
\eqn{W_{2,0}$,$W_{2,1}$,$W_{2,2}$,$W_{2,3}}{W(2,0),W(2,1),W(2,2),W(2,3)}, respectively.}}


and so forth.}
}
\seealso{
\code{\link{wavDWPT}}, \code{\link{wavBestBasis}}.}

\examples{
## calculate a 3-level DWPT of the sunspots series 
W <- wavDWPT(sunspots, n.level=3)

## extract the level 1 basis 
W12 <- wavPacketBasis(W, 1:2)

## obtain the names of the crystals that were 
## extracted: "w1.0" "w1.1" 
names(W12$data)

## extract basis corresponding to crystal set: 
## "w2.0" "w2.1" "w1.1". This set comprises a 
## split-level basis 
Wsplit <- wavPacketBasis(W, c(3,4,2))
names(Wsplit$data)
}
\keyword{univar}

