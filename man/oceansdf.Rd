%%  WARNING: This file was automatically generated from the associated 
%%  wav_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_util.mid file in the project MID directory. Once the 
%%  wav_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the oceansdf function

\name{oceansdf}
\alias{oceansdf}
\title{SDF for ocean series}
\concept{oceanutilities}
\usage{oceansdf(f)}
\description{Parametric spectral density function for \code{ocean} series.}

\arguments{
\item{f}{normalized frequencies over which the SDF function is evaluated.}

}

\value{
a numeric vector containing the values of ocean SDF function evaluated at
the specified normalized frequencies.
}
\seealso{
\code{\link{ocean}}.}

\examples{
oceansdf(c(0.25,0.375))
}
\keyword{datagen}

