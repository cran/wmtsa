%%  WARNING: This file was automatically generated from the associated 
%%  wav_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_util.mid file in the project MID directory. Once the 
%%  wav_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavMaxLevel function

\name{wavMaxLevel}
\alias{wavMaxLevel}
\title{Maximum decomposition level}
\concept{waveletutilities}
\usage{wavMaxLevel(n.taps=8, n.sample=1024, xform="modwt")}
\description{Interior wavelet coefficients are those not subject to circular filter
operations. This function calculates the maximum level for a wavelet
transform for which there exists at least one interior wavelet coefficient.}

\arguments{

\item{n.sample}{the number of points in the original time series. Default: \code{1024}.}

\item{n.taps}{the length of the wavelet filter. Default: \code{8}.}

\item{xform}{a character string denoting the transform type. Supported types
are \code{"dwt"}, \code{"dwpt"}, \code{"modwt"}, and \code{"modwpt"}. Default: \code{"modwt"}.}
}

\value{
an integer denoting the maximum decomposition level which
contains more than one interior wavelet coefficient.
}
\seealso{
\code{\link{wavDWT}}, \code{\link{wavMODWT}}, \code{\link{wavDWPT}}.}

\examples{
wavMaxLevel(n.taps=8, n.sample=1024, xform="modwt")
}
\keyword{univar}

