%%  WARNING: This file was automatically generated from the associated 
%%  wav_xform.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_xform.mid file in the project MID directory. Once the 
%%  wav_xform.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavMODWT function

\name{wavMODWT}
\alias{wavMODWT}
\title{The maximal overlap discrete wavelet transform (MODWT)}
\concept{wavelettransforms, maximal overlap discrete wavelet}
\usage{wavMODWT(x, wavelet="s8", n.levels=ilogb(length(x), base=2),
    position=list(from=1,by=1,units=character()), units=character(),
    title.data=character(), documentation=character(), keep.series=FALSE)}
\description{Let \eqn{j, t} be the decomposition level,
and time index, respectively, and
\eqn{s_{0,t}=X_{t=0}^{N-1}}{s(0,t)=X(t) for t=0,...,N-1} where
\eqn{X_t}{X(t)} is a real-valued uniformly-sampled time series. The
\eqn{j^{th}}{jth} level MODWT wavelet
coefficients \eqn{\tilde{d}_{j,t}}{d(j,t)}
and scaling coefficients \eqn{\tilde{s}_{j,t}}{s(j,t)}
are defined as \eqn{\tilde{d}_{j,t} \equiv \sum_{l=0}^{L-1}\widetilde{h}_l
\tilde{s}_{j-1, t - 2^{j-1}\;l \mbox{ mod }N,}}{d(j,t)=sum(h(l) s(j-1, t - 2^(j-1) l mod N))}
and \eqn{\tilde{s}_{j,t} \equiv \sum_{l=0}^{L-1}\widetilde{g}_l
\tilde{s}_{j-1, t - 2^{j-1}\;l \mbox{ mod }N}}{s(j,t)=sum(g(l) s(j-1, t - 2^(j-1) l mod N))}
The variable \eqn{L} is the length of both the scaling
filter \eqn{(g)} and wavelet filter
\eqn{(h)}. The \eqn{\tilde{d}_{j,t}}{d(j,t)}
and \eqn{\tilde{s}_{j,t}}{s(j,t)} are the wavelet and
scaling coefficients, respectively, at decomposition level
\eqn{j} and time index \eqn{t}.
The MODWT is a collection of all wavelet coefficients and the
scaling coefficients at the last level:
\eqn{\mathbf{\tilde{d}_1,\tilde{d}_2},\ldots,\mathbf{\tilde{d}_J,\tilde{s}_J}}{d(1),d(2),...,d(J),s(J)} where
\eqn{\mathbf{\tilde{d}_j}}{d(j)} and
\eqn{\mathbf{\tilde{s}_j}}{s(j)} denote a collection of wavelet
and scaling coefficients, respectively, at level \eqn{j}.}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series.}
\item{documentation}{a character string used to describe the input
\code{data}. Default: \code{character()}.}

\item{keep.series}{a logical value. If \code{TRUE}, the original series
is preserved in the output object. Default: \code{FALSE}.}

\item{n.levels}{the number of decomposition levels. Default: \code{as.integer(floor(logb(length(x),base=2)))}.}

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
The MODWT is a non-decimated form of the discrete wavelet transform (DWT)
having many advantages over the DWT including the ability
to handle arbitrary length sequences and shift invariance (while the
\code{wavDWT} function can handle arbitrary length
sequences, it does so by means of an ad hoc storage sytem for odd length
scaling coefficient crystals. The MODWT needs no such scheme and is
more robust in this respect). The cost of
the MODWT is in its redundancy. For an \eqn{N}
point input sequence, there are \eqn{N} wavelet
coefficients per scale. However, the number of multiplication operations is
\eqn{O(N \log_2(N))}{O(N log2(N))} which is the same as
the fast Fourier transform, and is acceptably fast for most situations.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{reconstruct}}, \code{\link{wavDaubechies}}, \code{\link{wavDWT}}, \code{\link{wavMODWPT}}, \code{\link{wavDictionary}}, \code{\link{wavIndex}}, \code{\link{wavTitle}}, \code{\link{wavBoundary}}.}

\examples{
## calculate the MODWT of linear chirp 
linchirp <- make.signal("linchirp", n=1024)
result   <- wavMODWT(linchirp, wavelet="s8", n.levels=5, keep.series=TRUE)

## plot the transform shifted for approximate zero 
## phase alignment 
plot(wavShift(result))

## plot summary 
eda.plot(result)

## summarize the transform 
summary(result)
}
\keyword{univar}

