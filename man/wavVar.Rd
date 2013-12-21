%%  WARNING: This file was automatically generated from the associated 
%%  wav_var.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_var.mid file in the project MID directory. Once the 
%%  wav_var.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavVar, plot.wavVar, print.wavVar, summary.wavVar, print.summary.wavVar functions

\name{wavVar}
\alias{wavVar}
\alias{plot.wavVar}
\alias{print.wavVar}
\alias{summary.wavVar}
\alias{print.summary.wavVar}
\title{Discrete wavelet variance estimation}
\concept{waveletwavelet variance}
\usage{wavVar(x, xform="modwt", wavelet="s8", n.levels=NULL,
    position=list(from=1,by=1,units=character()), units=character(),
    documentation=character(), sdf=NULL, sdfargs=NULL,
    sampling.interval=1, n.fft=1024)}
\description{The discrete wavelet variance is a useful alternative to the
spectral density function (SDF) and is seen as an octave-band
regularization of the SDF. The wavelet variance decomposes the
variance of certain stochastic processes on a scale-by-scale
basis, and thus, is very appealing to the analyst studying
physical phenomena which fluctuate both within and across a wide
range of scale.

By definition, the wavelet variance involves an averaged energy summation of
MODWT wavelet coefficients.
While DWT wavelet coefficients can also be used, the statistical
properties are inferior to those of the MODWT wavelet variance.
See the references for more details.

\bold{The MODWT Wavelet Variance}

Let \eqn{N} be the the number of samples in a
time series \eqn{\{X_t\}}{X(t)},
\eqn{L} be the length of the wavelet filter,
\eqn{L_j \equiv (2^j - 1)(L - 1) + 1}{L(j)=(2^j-1)(L-1)+1}
be the equivalent filter width at level $j$ in
a MODWT, and \eqn{\tau_j \equiv 2^{j-1}}{tau(j)=2^(j-1)}
be the scale of the data at level \eqn{j}
for \eqn{j=1,\ldots,J}{j=1,...,J}. Then the
unbiased wavelet variance is defined as
\deqn{ \hat\nu_X^2(\tau_j ) \equiv{1 \over M_j} \sum_{t=L_j - 1}^{N-1}
\tilde{d}_{j,t}^2,}{1/Mj * sum[t=L(j)-1, ..., N-1],}

where \eqn{\tilde{d}_{j,t}}{d(j,t)} are the
MODWT coefficients at level \eqn{j} and
time \eqn{t}, and
\eqn{M_j \equiv N - L_j + 1}{Mj=N - L(j) + 1}.
The unbiased wavelet variance
estimator avoids so-called boundary coefficients which are those
coefficients subject to circular filter operations in a discrete
wavelet transform. The biased estimator is defined as
\deqn{ \widetilde\nu_X^2(\tau_j ) \equiv{1 \over N} \sum_{t=0}^{N-1}
\tilde{d}_{j,t}^2,}{1/N * sum[t=0,...,N-1] (d(j,t))^2,}
and includes the effects of the boundary coefficients.

\bold{The DWT Wavelet Variance}

The DWT can also be used to calculate wavelet variance estimates,
but is not preferred over the MODWT due to its poor statistical
properties.
Let \eqn{N_j \equiv \lfloor N / 2^j \rfloor}{N(j)= floor(N/2^j)}
be the number of DWT wavelet coefficients at level
\eqn{j}, and \eqn{L'_j \equiv \lceil
(L-2)(1 - 2^{-j}) \rceil}{L'(j)=ceiling((L-2)(1-2^(-j)))} be the number of DWT
boundary coefficients at level \eqn{j} (assuming
\eqn{N_j > L'j}{N(j) > L'(j)}). Then the DWT version of the
unbiased wavelet variance is defined as
\deqn{ \hat{\hat\nu}_X^2(\tau_j ) \equiv{1 \over N - 2^jL'_j}
\sum_{t=L'_j - 1}^{N_j-1} d_{j,t}^2,}{1/(N-2^j * L'(j)) sum[t=L'(j)-1,...,N(j)-1] (d(j,t))^2,}

where \eqn{d_{j,t}}{d(j,t)} are the DWT
coefficients at level \eqn{j} and
time \eqn{t}.
Similarly, the DWT version of the biased wavelet variance is defined as
\deqn{ \widetilde{\widetilde\nu}_X^2(\tau_j ) \equiv{1 \over N}
\sum_{t=0}^{N_j-1} d_{j,t}^2.}{1/N sum[t=0,...,N(j)-1 (d(j,t))^2].}}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series.}
\item{documentation}{a character string used to describe the input
\code{data}. Default: \code{character()}.}

\item{n.fft}{a positive integer (greater than one) defining the number of frequencies to use in evaluating the SDF. The
frequencies are uniformly distributed over the interval [0, Nyquist] ala f=[0, 1/P , 2/P, 3/P, ..., (n.freq-1)/P]
where P=2*(n.freq-1)/sampling.interval. Only used when the input SDF is not \code{NULL}. Default: \code{1024}.}

\item{n.levels}{the number of decomposition levels. Default: the
maximum level at which there exists at least one interior wavelet coefficient.}

\item{position}{a \code{list} containing the arguments
\code{from, by} and \code{to} which describe the position(s) of the input
\code{data}. All position arguments need not be specified as missing members
will be filled in by their default values. Default: \code{list(from=1, by=1, units=character())}.}

\item{sampling.interval}{sampling interval of the time series. Default: \code{1}.}

\item{sdf}{a spectral density function of the process corresponding to the input time series.
This input must be a function whose first argument is \code{f}
(representing frequency). At a minimum, the SDF must be defined over frequencies [0, Nyquist]
where \code{Nyquist=1/(2*sampling.interval)}. Additional arguments that are needed to calculate the SDF should
be passed via the \code{sdfargs} parameter.
This argument is used only for calculating mode 2 EDOF.
If the EDOF mode 2 estimates are not desired, specify this
this argument as \code{NULL} and the EDOF mode 2 and corresponding confidence
intervals will not be calculated. See the \code{\link[ifultools]{mutilsSDF}} function for more details. Default: \code{NULL}.}

\item{sdfargs}{a list of arguments passed directly to the SDF function ala
\code{do.call}. Default: \code{NULL} (no additional arguments).}

\item{units}{a string denoting the units of the time series. Default: \code{character()} (no units).}

\item{wavelet}{a character string denoting the filter type.
See \code{wavDaubechies} for details. Default: \code{"s8"}.}

\item{xform}{a character string denoting the type of wavelet transform:
"modwt" or "dwt". Default: \code{"modwt"}.}
}

\value{
an object of class \code{wavVar}.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavVarTest}}, \code{\link{wavEDOF}}, \code{\link{wavVarConfidence}}, \code{\link[ifultools]{mutilsSDF}}.}

\examples{
## create sequence 
x <- make.signal("doppler")

## perform a time independent wavelet variance 
## analysis 
vmod <- wavVar(x)

## plot the result 
plot(vmod, pch=15, title="Wavelet Variance of Doppler")

## calculate wavelet variance estimaates for the 
## ocean series and calculate EDOF mode 2 
## estimates and corresponding 95 percent 
## confidence intervals 
vocean <- wavVar(ocean, sdf=oceansdf, wavelet="d6")

## summarize the results 
plot(vocean, edof=1:3)

summary(vocean)
}
\keyword{univar}

