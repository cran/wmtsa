%%  WARNING: This file was automatically generated from the associated 
%%  wav_sig.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_sig.mid file in the project MID directory. Once the 
%%  wav_sig.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the make.signal function

\name{make.signal}
\alias{make.signal}
\title{Test signal generation}
\concept{series generation}
\usage{make.signal(name, n=1024, snr=Inf)}
\description{Generates various test signals for wavelet transforms.}

\arguments{
\item{name}{a character string denoting the type of test signal
to create. Supported values are:
\code{"dirac"}, \code{"kronecker"}, \code{"heavisine"}, \code{"bumps"}, \code{"blocks"},
\code{"doppler"}, \code{"ramp"}, \code{"cusp"}, \code{"crease"}, \code{"sing"}, \code{"hisine"},
\code{"losine"}, \code{"linchirp"}, \code{"twochirp"}, \code{"quadchirp"},
\code{"mishmash1"}, \code{"mishmash2"}, \code{"mishmash3"}, \code{"levelshift"},
\code{"jumpsine"}, \code{"gauss"}, \code{"patches"},
\code{"linear"}, \code{"quadratic"}, and \code{"cubic"}.}
\item{n}{an integer specifying the length out the output series.
Default: \code{1024}.}

\item{snr}{a numeric value representing the approximate
signal to noise ratio of the output. Default: \code{Inf} (no noise).}
}

\value{
a vector of numeric values containing the resulting test series.
}
\seealso{
\code{\link{wavDWT}}, \code{\link{wavMODWT}}.}

\examples{
nms <- c("blocks", "linchirp", "mishmash1", "bumps")
z <- lapply(nms, make.signal)
stackPlot(x=seq(1024),y=z, ylab=nms)
}
\keyword{ts}

