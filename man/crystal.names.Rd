%%  WARNING: This file was automatically generated from the associated 
%%  wav_gnrc.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_gnrc.mid file in the project MID directory. Once the 
%%  wav_gnrc.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the crystal.names function

\name{crystal.names}
\alias{crystal.names}
\title{Generic function for obtaining crystal names of wavelet transform objects}
\concept{wavelet}
\usage{crystal.names(x, ...)}
\description{Returns the names assigned to crystals (vectors of wavelet transform coefficients).}

\arguments{
\item{x}{any object. Missing values (NAs) are allowed.}
\item{...}{optional arguments to be passed directly to the
inherited function without alteration and with the original names preserved.}
}

\value{
a vector of character strings.
}
\seealso{
\code{\link{wavMRD}}.}

\examples{
methods(crystal.names)
}
\keyword{methods}

