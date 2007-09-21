%%  WARNING: This file was automatically generated from the associated 
%%  wav_gnrc.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_gnrc.mid file in the project MID directory. Once the 
%%  wav_gnrc.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavStackPlot function

\name{wavStackPlot}
\alias{wavStackPlot}
\title{Generic function for generating stacked plots}
\concept{wavelet}
\usage{wavStackPlot(x, ...)}
\description{Stackplots stack each plot on top of one another without having to
manually divide the ploting region via \code{par("mfrow")} or other.}

\arguments{
\item{x}{any object. Missing values ( NAs) are allowed.}
\item{...}{optional arguments to be passed directly to the
inherited function without alteration and with the original names preserved.}
}
\note{An extended data analysis plot is shown.}

\seealso{
\code{\link{wavMRD}}, \code{\link{wavTransform}}, \code{\link{eda.plot}}, \code{\link{crystal.names}}.}

\examples{
methods(wavStackPlot)
}
\keyword{methods}

