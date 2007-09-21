%%  WARNING: This file was automatically generated from the associated 
%%  wav_boot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_boot.mid file in the project MID directory. Once the 
%%  wav_boot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavBootstrap function

\name{wavBootstrap}
\alias{wavBootstrap}
\title{Adaptive wavelet-based bootstrapping}
\concept{waveletsimulationsurrogate data generationbootstrapping}
\usage{wavBootstrap(x, white.indices=wavDWPTWhitest(x),
    n.realization=1, wavelet="s8", n.level=NULL)}
\description{Given a set of indices which represent the whitest transform available
in a DWPT, this function randomizes the coefficients in each of the crystals
comprising the transform (via random selection with replacement) followed
by an inverse transform. The z is a bootstrapped version of the original
time series.}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series or an
object of class \code{wavTransform} as output by the \code{wavDWPT} function.}
\item{n.level}{the number of decomposition levels.  This argument is used only if
\code{x} is a time series. Default: \code{floor(logb(length(x), base=2)) - 2}.}

\item{n.realization}{the number of realizations to generate. Default: \code{1}.}

\item{wavelet}{a character string denoting the filter type.
See \code{wavDaubechies} for details. This argument is used only if
\code{x} is a time series. Default: \code{"s8"}.}

\item{white.indices}{a \code{list} containing the \code{level} and \code{osc} vectors denoting
the level and oscillation index, respectively, of the whitest transform. Default: \code{wavDWPTWhitest(x)}.}
}

\value{
a list of numeric vectors containing the bootstrapped series. If \code{n.realization=1},
the the output is a numeric vector (not packed into a \code{list}).
}
\references{
D. B. Percival, S. Sardy and A. C. Davison,
\emph{Wavestrapping Time Series: Adaptive Wavelet-Based Bootstrapping},
in W. J. Fitzgerald, R. L. Smith, A. T. Walden and P. C. Young (Eds.),
\emph{Nonlinear and Nonstationary Signal Processing},
Cambridge, England: Cambridge University Press, 2001.
}

\seealso{
\code{\link{wavDWPT}}, \code{\link{wavDWPTWhitest}}.}

\examples{
## wavestrap the sunspots series 
x <- as.numeric(sunspots)
z <- wavBootstrap(x, n.realization=1)

stackPlot(x=seq(along=sunspots),
y=data.frame(x, z, abs(z)),
ylab=list(text=c("sunspots","wavestrap","|wavestrap|")))

title("Wavelet-based bootstrapping of sunspots series", cex=0.7)
}
\keyword{datagen}

