%%  WARNING: This file was automatically generated from the associated 
%%  wav_xform.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_xform.mid file in the project MID directory. Once the 
%%  wav_xform.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavDWT function

\name{wavDWT}
\alias{wavDWT}
\title{The discrete wavelet transform (DWT)}
\concept{wavelettransforms, discrete wavelet}
\usage{wavDWT(x, n.levels=ilogb(length(x), base=2),
    wavelet="s8", position=list(from=1,by=1,units=character()), units=character(),
    title.data=character(), documentation=character(), keep.series=FALSE)}
\description{The discrete wavelet transform using convolution style filtering
and periodic extension.

Let \eqn{j, t} be the decomposition level,
and time index, respectively, and
\eqn{s_{0,t}=X_{t=0}^{N-1}}{s(0,t)=X(t) for t=0,...,N-1} where
\eqn{X_t}{X(t)} is a real-valued uniformly-sampled time series. The
\eqn{j^{th}}{jth} level DWT wavelet
coefficients (\eqn{d_{j,t}}{d(j,t)})
and scaling coefficients (\eqn{s_{j,t}}{s(j,t)})
are defined as \eqn{d_{j,t} \equiv \sum_{l=0}^{L-1} h_l s_{j-1,2t+1-l \bmod N_{j-1}},
\quad t=0,\ldots, N_j -1}{d(j,t)=sum(h(l) s(j-1, t - 2t+1-l) mod N(j-1))} and
\eqn{s_{j,t} \equiv \sum_{l=0}^{L-1} g_l s_{j-1,2t+1-l \bmod N_{j-1}},
\quad t=0,\ldots, N_j -1}{s(j,t)=sum(g(l) s(j-1, t - 2t+1-l mod N(j-1)))}
for \eqn{j=1,\ldots,J}{j=1,...,J} where \eqn{\{ h_l \}}{h(l)} and \eqn{\{ g_l \}}{g(l)} are the \eqn{j^{th}}{jth} level wavelet and scaling filter, respectively, and
\eqn{N_j \equiv N / 2^j}{Nj=2^(j-1)}.  The DWT is a collection of all wavelet coefficients and the
scaling coefficients at the last level:
\eqn{\mathbf{d_1,d_2},\ldots,\mathbf{d_J,s_J}}{d(1),d(2),...,d(J),s(J)} where
\eqn{\mathbf{d_j}}{d(j)} and
\eqn{\mathbf{s_j}}{s(j)} denote a collection of wavelet
and scaling coefficients, respectively, at level \eqn{j}.}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series.}
\item{documentation}{a character string used to describe the input
\code{data}. Default: \code{character()}.}

\item{keep.series}{a logical value. If \code{TRUE}, the original series
is preserved in the output object. Default: \code{FALSE}.}

\item{n.levels}{the number of decomposition levels.
Default: \code{as.integer(floor(logb(length(x),base=2)))}.}

\item{position}{a \code{list} containing the arguments
\code{from, by} and \code{to} which describe the position(s) of the input
\code{data}. All position arguments need not be specified as missing members
will be filled in by their default values. Default: \code{list(from=1, by=1, units=character())}.}

\item{title.data}{a character string representing the name of the input
\code{data}. Default: \code{character()}.}

\item{units}{a string denoting the units of the time series. Default: \code{character()} (no units).}

\item{wavelet}{a character string denoting the filter type.
See \code{wavDaubechies} for details. Default: \code{"s8"}.}
}

\value{
an object of class \code{wavTransform}.
}
\details{
%
This DWT imposes an ad hoc storage sytem for odd length
scaling coefficient crystals: if the length of a scaling
coefficient crystal is odd, the last coefficient is "stored"
in the \emph{extra} crystal. During
reconstruction, any extra scaling coefficients are returned to
their proper location. Such as system imposes no spurious
energy in the transform coefficients at the cost of a little
bookkeeping.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{reconstruct}}, \code{\link{wavDaubechies}}, \code{\link{wavMODWT}}, \code{\link{wavMODWPT}}, \code{\link{wavMRD}}, \code{\link{wavDictionary}}, \code{\link{wavIndex}}, \code{\link{wavTitle}}, \code{\link{wavBoundary}}, \code{\link{wavShrink}}.}

\examples{
## calculate the DWT of linear chirp 
linchirp <- make.signal("linchirp", n=1024)
result   <- wavDWT(linchirp, wavelet="s8", n.levels=5, keep.series=TRUE)

## plot the transform shifted for approximate zero 
## phase alignment 
plot(wavShift(result))

## plot summary 
eda.plot(result)

## summarize the transform 
summary(result)
}
\keyword{univar}

