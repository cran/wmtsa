%%  WARNING: This file was automatically generated from the associated 
%%  wav_fdp.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_fdp.mid file in the project MID directory. Once the 
%%  wav_fdp.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavFDPSDF function

\name{wavFDPSDF}
\alias{wavFDPSDF}
\title{Spectral density function for a fractionally differenced process}
\concept{fractionally differenced (FD) processspectral density function generation}
\usage{wavFDPSDF(f, delta=0.45, variance=1, response=NULL)}
\description{Returns the spectral density function (SDF) for a
fractionally differenced (FD) process. Given a unit sampling rate, the
SDF for an FD proces is
\deqn{ {\sigma_\varepsilon^2 \over |2 \sin(\pi f)|^{2 \delta}}, }{variance / abs(2 * sin(pi*f))^(2 * delta),}
where \eqn{\sigma_\varepsilon^2}{variance} is the
innovations variance, \eqn{\delta}{delta} is
the FD parameter, and \eqn{f}
is the normalized frequency for \eqn{|f| < 1/2}{|f| < 1/2}.}

\arguments{
\item{f}{a numeric value representing normalized frequency
where the sampling interval is unity.}
\item{delta}{the FD parameter. Default: \code{0.45}.}

\item{response}{a \code{list} containing the objects
\code{frequency} and
\code{sqrgain} which represent, respectively, a numeric
normalized frequency vector corresponding to a wavelet squared gain
response at a particular wavelet decomposition level. This argument
typically will not be set by the user. Rather, it is used internally
by FD process maximum likelihood estimators. Default: \code{NULL}.}

\item{variance}{the FD innovations variance. Default: \code{1}.}
}

\value{
the SDF values corresponding to the FD model parameters.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000, 340--92.
}

\seealso{
\code{\link{wavFDPBand}}, \code{\link{wavFDPBlock}}, \code{\link{wavFDPTime}}.}

\examples{
## create a normalized frequency vector 
f <- seq(from=1e-2, to=1/2, length=100)

## calculate the FDP SDF for delta=0.45 and unit 
## innovations variance 
S <- wavFDPSDF(f, delta=0.45, variance=1)

## plot the results 
plot(f, S,log="xy", xlab="Frequency", ylab="SDF of FDP(0.45, 1)")
}
\keyword{models}

