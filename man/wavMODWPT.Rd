%%  WARNING: This file was automatically generated from the associated 
%%  wav_xform.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_xform.mid file in the project MID directory. Once the 
%%  wav_xform.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavMODWPT function

\name{wavMODWPT}
\alias{wavMODWPT}
\title{The maximal overlap discrete wavelet packet transform (MODWPT)}
\concept{wavelettransforms, maximal overlap discrete wavelet packet}
\usage{wavMODWPT(x, wavelet="s8", n.levels=ilogb(length(x), base=2),
    position=list(from=1,by=1,units=character()), units=character(),
    title.data=character(), documentation=character())}
\description{Given \eqn{j, n, t} are the decomposition level,
oscillation index, and time index, respectively, the MODWPT is given by
\deqn{
\widetilde{W}_{j,n,t} \equiv \sum_{l=0}^{L-1}\widetilde{u}_{n,l}
\widetilde{W}_{j-1,\lfloor n/2 \rfloor, t - 2^{j-1}\;l \mbox{ mod }N}
}{W(j,n,t)=sum(u(n,l) * W(j-1, floor(n/2), t - 2^(j-1) * l mod N))}

The variable \eqn{L} is the length of the filters defined by
\deqn{
\widetilde{u}_{n,l} \equiv \left\{
\begin{array}{ll}
\widetilde{g_l}/\surd{2}, & \mbox{ if }n \bmod 4=0 \mbox{ or }3; \\
\widetilde{h_l}/\surd{2}, & \mbox{ if }n \bmod 4=1 \mbox{ or }2,
\end{array}\right.
}{u(n,l)=g(l) / sqrt(2) if n mod 4=0 or 3; u(n,l)=h(l) / sqrt(2) if n mod 4=1 or 2; for l=0, ..., L-1}

where \eqn{g} and \eqn{h} are the scaling filter and wavelet filter, respectively.
By definition, \eqn{\widetilde{W}_{0,0,t} \equiv X_t}{W(0,0,t)=X(t)} where
\eqn{\{X_t\}}{X} is the original time series.}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series.}
\item{documentation}{a character string used to describe the input
\code{data}. Default: \code{character()}.}

\item{n.levels}{the number of decomposition levels.
Default: \code{as.integer(floor(logb(length(x),base=2)))}.}

\item{position}{a \code{list} containing the arguments
\code{from, by} and \code{to} which describe the position(s) of the input
\code{data}. All position arguments need not be specified as missing members
will be filled in by their default values. Default: \code{list(from=1, by=1, units=character())}.}

\item{title.data}{a character string representing the name of the input
\code{data}. Default: \code{character()}.}

\item{units}{a string denoting the units of the time series. Default: \code{character()} (no units).}

\item{wavelet}{a character string denoting the filter type. See \code{wavDaubechies} for details. Default: \code{"s8"}.}
}

\value{
an object of class \code{wavTransform}.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{reconstruct}}, \code{\link{wavMRD}}, \code{\link{wavMODWT}}, \code{\link{wavDWT}}, \code{\link{wavDWPT}}, \code{\link{wavDaubechies}}, \code{\link{wavShift}}, \code{\link{wavZeroPhase}}.}

\examples{
## calculate the MODWPT of sunspots series out to 
## 3 levels using Daubechies least asymmetric 
## 8-tap filter set 
z <- wavMODWPT(sunspots, wavelet="s8", n.levels=3)

## plot the transform 
plot(z)

## summarize the transform 
summary(z)
}
\keyword{univar}

