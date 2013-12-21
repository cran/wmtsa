%%  WARNING: This file was automatically generated from the associated 
%%  wav_gnrc.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_gnrc.mid file in the project MID directory. Once the 
%%  wav_gnrc.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the eda.plot function

\name{eda.plot}
\alias{eda.plot}
\title{Generic function for generating extended data analysis plots}
\concept{wavelet}
\usage{eda.plot(x, ...)}
\description{Data analysis plots are used to visually summarize the salient features
of the output and typically involve a combination of plots
in a single plot frame.}

\arguments{
\item{x}{any object. Missing values ( NAs) are allowed.}
\item{...}{optional arguments to be passed directly to the
inherited function without alteration and with the original names preserved.}
}
\note{An extended data analysis plot is shown.}

\seealso{
\code{\link{wavMRD}}, \code{\link{wavTransform}}, \code{\link{wavStackPlot}}, \code{\link{crystal.names}}.}

\examples{
methods(eda.plot)
}
\keyword{methods}

