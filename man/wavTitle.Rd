%%  WARNING: This file was automatically generated from the associated 
%%  wav_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_util.mid file in the project MID directory. Once the 
%%  wav_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavTitle function

\name{wavTitle}
\alias{wavTitle}
\title{Extract the name of the data used to generate objects of various wavelet classes}
\concept{wavelet}
\usage{wavTitle(x, default="x")}
\description{Wavelet functions store the original name of the data used to create the output in various
locations within the output object. This function provides a means by which the user can
directly access data name.}

\arguments{
\item{x}{an object of class \code{mra}, \code{mrd}, \code{decompose},
\code{wavTransform} or \code{signalSeries}.}
\item{default}{a default chacter string to use if no valid time
series name is found. Default: \code{"x"}.}
}

\value{
a character string vector containing the result.
}
\seealso{
\code{\link{wavDWT}}, \code{\link{wavMODWT}}.}

\examples{
wavTitle(wavDWT(1:8))
}
\keyword{attribute}

