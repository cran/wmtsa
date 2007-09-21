%%  WARNING: This file was automatically generated from the associated 
%%  wav_data.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_data.mid file in the project MID directory. Once the 
%%  wav_data.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the atomclock function

\name{atomclock}
\alias{atomclock}
\title{Cesium Beam Atomic Clock Data}
\concept{measured time series}
\description{This series represents the difference in time between a cesium
beam atomic clock and an official time scale known as UTC(USNO)
maintained by the US Naval Observatory, Washington DC. The UTC
portion of the USNO series refers to coordinate universal time
which is used as an international time standard. Negative values
in the resulting (difference) series represents a lag in time
relative to the UTC(USNO) standard. The differences in time were
recorded at the same time for consecutive days in the 1970s
resulting in a sampling interval of 1 day. The amplitude units
are expressed in microseconds.}

\seealso{
\code{\link{nile}}, \code{\link{ocean}}, \code{\link{ecg}}, \code{\link{fdp045}}.}

\examples{
plot(atomclock)
}
\keyword{ts}

