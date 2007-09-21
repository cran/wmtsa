%%  WARNING: This file was automatically generated from the associated 
%%  wav_var.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_var.mid file in the project MID directory. Once the 
%%  wav_var.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavEDOF function

\name{wavEDOF}
\alias{wavEDOF}
\title{Equivalent degrees of freedom (EDOF) estimates for a chi-squared distribution}
\concept{waveletwavelet varianceequivalent degrees of freedom}
\usage{wavEDOF(x, wavelet="s8", levels=NULL, sdf=NULL, sdfargs=NULL,
    sampling.interval=1, n.fft=1024)}
\description{Let X be a collection of M uncorrelated zero mean Gaussian random variables (RVs).
The sum of the squares of the RVs in X will obey a scaled chi-square distribution
with M degrees of freedom (DOF). If, however, the original Gaussian RVs are (partially)
correlated, we can approximate the distribution of the sum of the squares of
(correlated Gaussian) RVs using a scaled chi-square distribution with the
DOF adjusted for the correlation in the RVs. These adjusted
DOF estimates are known as the
\emph{equivalent degrees of freedom}
(EDOF). In the context of unbiased wavelet variance analysis, the EDOF can be
used to estimate confidence intervals that are guaranteed to have non-negative
bounds.

This program calculates three estimates of the
EDOF for each level of a discrete wavelet transform.
The three modes are described as follows for the MODWT of an
an input sequence \eqn{\{X_t\}_{t=0}^{N-1}}{X(t)}:

\describe{
\item{EDOF 1}{Large sample approximation that requires an SDF estimation via wavelet coefficients.
\deqn{\hat\eta_1 ={ M_j (\hat{s}_{j,0})^2 \over \hat{A}_j},}{n1 = s(j,0)^2 / Aj,}
where \eqn{\hat{s}_{j,\tau}}{s(j,tau)} is the autocovariance sequence defined by
\deqn{ \hat{s}_{j,\tau} \equiv{1 \over M_j} \sum_{t=0}^{M_j - 1}
\tilde{d}_{j,t}^{(int)} \tilde{d}_{j,t + |\tau|}^{(int)}
\;\; 0\le|\tau|\le M_j-1, }{s(j,tau) = (1 / Mj) * sum[t=0,..., Mj - 1]{d(j,t)}}
and
\eqn{\tilde{d}_{j,t}^{(int)}}{d(j,t)} are the
\eqn{M_j$ $j^{th}}{Mj jth} level interior MODWT wavelet coefficients and
\eqn{\hat{A}_j}{Aj} is defined as
\deqn{ \hat{A}_j \equiv{ ( \hat{s}_{j,0} ) \over 2 } + \sum_{\tau=1}^{M_j-1} ( \hat{s}_{j,\tau} )^2.}{Aj = s(j,0)/2 + sum[tau=1,...,Mj-1]s(j,tau)^2.}}

\item{EDOF 2}{Large sample approximation where the SDF is known a priori.
\deqn{\hat\eta_2 ={ 2{\biggl( \sum_{k=1}^{\lfloor (M_j - 1)/ 2 \rfloor} C_j(f_k)
\biggr)}^2  \over \sum_{k=1}^{\lfloor (M_j - 1) / 2 \rfloor} C_j^2(f_k)},}{n2 = 2 * (sum[k=1,...,floor((Mj-1)/2)] Cj(f(k)))^2 / sum[k=1,...,floor((Mj-1)/2)] (Cj(f(k)))^2,}
where \eqn{f_k \equiv k /M_j}{f(k)=k/Mj} and
\eqn{C_j \equiv \widetilde{\mathcal{H}}_j^{(D)}(f) S_X(f)}{Cj = Hj(f) Sx(f)}
is the product of Daubechies wavelet filter squared
gain function and the spectral density function of \eqn{X_t}{X(t)}.}

\item{EDOF 3}{Large sample approximation using a band-pass approximation for the SDF.
\deqn{\hat\eta_3=\mbox{max}\{M_j/2^j, 1\}}{n3 = max(Mj/2, 1)}.}}


See references for more details.}

\arguments{
\item{x}{an object of class \code{wavTransform} or a vector containing a uniformly-sampled real-valued time series.}
\item{levels}{a vector containing the decomposition levels. Default: when
\code{x} is of class \code{wavTransform} then \code{levels} is set to \code{1:n.level}, otherwise
\code{levels} is set to \code{1:J}, where \code{J} is the
maximum wavelet transform level in which there exists at least one interior wavelet coefficient.}

\item{n.fft}{a positive integer (greater than one) defining the number of frequencies to use in evaluating the SDF for
EDOF 2 calculations. The
frequencies are uniformly distributed over the interval [0, Nyquist] ala f=[0, 1/P , 2/P, 3/P, ..., (n.freq-1)/P]
where P=2*(n.freq-1)/sampling.interval. Only used when the input SDF is not \code{NULL}. Default: \code{1024}.}

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

\item{wavelet}{a character string denoting the filter type.
See \code{wavDaubechies} for details.
Only used if input \code{x} is a time series. Default: \code{"s8"}.}
}

\value{
a list containing the EDOF estimates for
modes 1, 2 and 3 as well as the block-dependent unbiased wavelet variance estimates.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavVar}}, \code{\link{wavVarConfidence}}, \code{\link[ifultools]{mutilsSDF}}.}

\examples{
## initialize variables 
n.level <- 9
wavelet <- "d6"
N       <- 4096
phi     <- 0.9

## define input SDF 
S <- function(f, phi)  1/(1 + phi^2 - 2*phi*cos(2*pi*f))
sdfarg <- list(phi=phi)

## create series and MODWT 
set.seed(100)
x <- rnorm(N)
W <- wavMODWT(x, wavelet=wavelet, n.level=n.level)

## calculate EDOF using the wavTransform object 
z1 <- wavEDOF(W, sdf=S, sdfarg=sdfarg)
print(z1)

## calculate EDOF using original time series 
z2 <- wavEDOF(x, wavelet=wavelet, levels=seq(n.level), sdf=S, sdfarg=sdfarg)
print(z2)

## compare the two approaches 
print(all.equal(z1,z2))
}
\keyword{distribution}

