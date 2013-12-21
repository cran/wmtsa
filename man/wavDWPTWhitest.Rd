%%  WARNING: This file was automatically generated from the associated 
%%  wav_boot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_boot.mid file in the project MID directory. Once the 
%%  wav_boot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavDWPTWhitest function

\name{wavDWPTWhitest}
\alias{wavDWPTWhitest}
\title{Seeks the whitest transform of a discrete wavelet packet transform (DWPT)}
\concept{waveletsimulationsurrogate data generationwavelet packetsbootstrapping}
\usage{wavDWPTWhitest(x, significance=0.05, test="port2", wavelet="s8", n.level=NULL)}
\description{This function seeks the whitest orthonormal transform of a DWPT.
The goal is to segment the normalized frequency interval
\eqn{f \in [0, 1/2]}{[0, 1/2]} into subintervals such that,
within each subinterval, the variability of the (corresponding)
spectral density function (SDF) is minimized, i.e., each segment of the SDF is
as flat as possible. Given an \eqn{N}-point uniformly sampled time series \bold{X},
and denoting \eqn{\mathbf{W}_{j,n}}{W(j,n)} as the
DWPT crystal at level \eqn{j} and (sequency ordered) oscillation index \eqn{n},
this optimization is achieved as follows:

\describe{
\item{1}{Perform a level \eqn{J - 2} partial DWPT of \bold{X} where
\eqn{J=\lfloor \log_2 (\mathbf{X}) \rfloor}{J=floor(log2(N))}.
By definition, \eqn{\mathbf{W}_{0,0} \equiv \mathbf{X}}{W(0,0)=X}. Begin
step \code{2} with \eqn{j=n=0}.}

\item{2}{Perform a white noise test on the current (parent) crystal: \eqn{\mathbf{W}_{j,n}}{W(j,n)}.
If it passes (or the current crystal is in the last decomposition level) retain the crystal.
Otherwise, discard the current parent crystal
and perform the white noise test on its children: \eqn{\mathbf{W}_{j+1,2n}}{W(j+1,2n)} and
\eqn{\mathbf{W}_{j+1,2n+1}}{W(j+1,2n+1)}.}

\item{3}{Repeat step \code{2} as many times as necessary until a suitable transform is found.}}}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series or an
object of class \code{wavTransform}.}
\item{n.level}{the number of decomposition levels.  This argument is used only if
\code{x} is a time series. Default: \code{floor(logb(length(x), base=2)) - 2}.}

\item{significance}{a numeric value on the interval (0,1)
which qualitatively signifies the fraction of times that the
white noise hypothesis is incorrectly rejected.
The \code{significance} is used to calculate comparative chi-square distribution
\eqn{p} x \eqn{100} percentage points where \eqn{p=1 - \mbox{significance}}{p=1 - significance}
(the chi-square degrees of freedom are estimated automatically within the specified white noise test).
Default: \code{0.05}.}

\item{test}{a character string denoting the white noise test to use.
Options are \code{"port1"}, \code{"port2"}, \code{"port3"} and \code{"cumper"} respresenting
the Portmanteau I, II, III and cumulative periodogram tests, respectively. See the
reference(s) for more details. Default: \code{"port2"}.}

\item{wavelet}{a character string denoting the filter type.
See \code{wavDaubechies} for details. This argument is used only if
\code{x} is a time series. Default: \code{"s8"}.}
}

\value{
a \code{list} containing the \code{level} and \code{osc} vectors denoting
the level and oscillation index, respectively, of the whitest transform.
}
\references{
D. B. Percival, S. Sardy and A. C. Davison,
\emph{Wavestrapping Time Series: Adaptive Wavelet-Based Bootstrapping},
in W. J. Fitzgerald, R. L. Smith, A. T. Walden and P. C. Young (Eds.),
\emph{Nonlinear and Nonstationary Signal Processing},
Cambridge, England: Cambridge University Press, 2001.
}

\seealso{
\code{\link{wavDWPT}}, \code{\link{wavBootstrap}}.}

\examples{
## calculate the DWPT of the sunspots series 
W <- wavDWPT(as.numeric(sunspots), wavelet="s8", n.levels=9)

## find the whitest transform based on the 
## Portmanteau I white noise test 
z <- wavDWPTWhitest(W, test="port1")

print(z)
}
\keyword{datagen}

