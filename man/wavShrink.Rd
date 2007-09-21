%%  WARNING: This file was automatically generated from the associated 
%%  wav_shrk.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_shrk.mid file in the project MID directory. Once the 
%%  wav_shrk.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavShrink function

\name{wavShrink}
\alias{wavShrink}
\title{Nonlinear denoising via wavelet shrinkage}
\concept{nonlinear noise reductionwavelet}
\usage{wavShrink(x, wavelet="s8",
    n.level=ilogb(length(x), base=2),
    shrink.fun="hard", thresh.fun="universal", threshold=NULL,
    thresh.scale=1, xform="modwt", noise.variance=-1.0,
    reflect=TRUE)}
\description{Performs a decimated or undecimated discrete wavelet transform
on the input series and "shrinks" (decreases the amplitude towards zero) the
wavelet coefficients based on a calculated noise threshold
and specified shrinkage function. The resulting shrunken set of wavelet
transform coefficients is inverted in a synthesis operation, resulting in
a denoised version of the original series.}

\arguments{
\item{x}{a vector containing a uniformly-sampled real-valued time series.}
\item{n.level}{the number of decomposition levels, limited to
\code{floor(logb(length(x),2))}.
Default: \code{floor(logb(length(x),2))}.}

\item{noise.variance}{a numeric scalar representing (an estimate of)
the additive Gaussian white noise variance. If unknown, setting
this value to 0.0 (or less) will prompt the function to automatically
estimate the noise variance based on the median absolute deviation (MAD)
of the scale one wavelet coefficients. Default: -1.}

\item{reflect}{a logical value. If \code{TRUE}, the
last \eqn{L_J = (2^{\mbox{n.level}} - 1)(L - 1) + 1}{Lj = (2^n.level - 1)(L - 1) + 1}
coefficients of the series are reflected (reversed and appended to the end
of the series) in order to attenuate the adverse effect of circular
filter operations on wavelet transform coefficients for
series whose endpoint levels are (highly) mismatched.
The variable \eqn{L_J}{Lj} represents the effective filter length at
decomposition level \code{n.level}, where \eqn{L}
is the length of the wavelet (or scaling) filter.
After waveshrinking and reconstructing,
the first \eqn{N} points of the result are returned, where \eqn{N} is the length
of the original time series.
Default: \code{TRUE}.}

\item{shrink.fun}{a character string denoting the shrinkage function.
Choices are \code{"hard"}, \code{"soft"}, and \code{"mid"}. Default: \code{"hard"}.}

\item{thresh.fun}{a character string denoting the threshold function to
use in calculating the waveshrink thresholds.
\describe{
\item{character string}{Choices are \code{"universal"},
\code{"minimax"}, and \code{"adaptive"}.}

\item{numeric values}{Either a single threshold value or a
vector of values containing \code{n.levels}
thresholds (one threshold per decomposition level).}}

\bold{Note:} if \code{xform == "modwt"}, then only the \code{"universal"}
threshold function is (currently) supported.
Default: \code{"universal"}.}

\item{thresh.scale}{a positive valued numeric scalar which
is used to amplify or attenuate the threshold values at each
decomposition level. The use of this argument signifies a departure
from a model driven estimate of the thresholds and can be used
to tweak the levels to obtain a smoother or rougher result. Default: 1.}

\item{threshold}{explicit setting of the wavelet shrinkage thresholds,
one for each level of the decomposition. If a single threshold is given, it is
replicated appropriately and (if the chosen transform is additionally a MODWT then) these thresholds
are normalized by dividing the threshold at level \eqn{j} by \eqn{2^((j-1)/2)}{2\eqn{\mbox{\textasciicircum}}{^}((j-1)/2)}.
If the number of thresholds is equal to the number of decomposition levels, the thresholds
are unaltered prior to use.
Default: \code{NULL} (thresholds are calculated based on the model defined by the \code{thresh.fun}
and \code{thresh.scale} input arguments).}

\item{wavelet}{a character string denoting the filter type.
See \code{wavDaubechies} for details. Default: \code{"s8"}.}

\item{xform}{a character string denoting the wavelet transform type.
Choices are \code{"dwt"} and \code{"modwt"} for the discrete wavelet transform (DWT)
and maximal overlap DWT (MODWT), respectively. The DWT is a decimated transform
where (at each level) the number of transform coefficients is halved. Given
\eqn{N} is the length of the original time series, the total
number of DWT transform coefficients is \eqn{N}.
The MODWT is a non-decimated transform where the number of
coefficients at each level is \eqn{N} and the total number of
transform coefficients is \eqn{N}*\code{n.level}. Unlike the DWT, the MODWT
is shift-invariant and is seen as a weighted average of all possible
non-redundant shifts of the DWT. See the references for details.
Default: \code{"modwt"}.}
}

\value{
vector containing the denoised series.
}
\details{
%
Assume that an appropriate model for our time series is
\eqn{\mathbf{X}=\mathbf{D} + \mathbf{\epsilon}}{X=D + e} where
\eqn{\mathbf{D}}{D} represents an unknown
deterministic signal of interest and
\eqn{\mathbf{\epsilon}}{e} is some
undesired stochastic noise that is independent and identically
distributed and has a process mean of zero. Waveshrink seeks to
eliminate the noise component
\eqn{\mathbf{\epsilon}}{e} of \eqn{\mathbf{X}}{X} in hopes of obtaining
(a close approximation to) \eqn{\mathbf{D}}{D}. The basic algorithm works as
follows:

\describe{
\item{1}{Calculate the DWT of \eqn{X}.}

\item{2}{Shrink (reduce towards zero) the wavelet coefficients based on a
selected thresholding scheme.}

\item{3}{Invert the DWT.}}


This function support different shrinkage methods and threshold estimation
schemes. Let \eqn{W} represent an arbitrary DWT coefficient and
\eqn{W^{\mbox{(t)}}}{W'} the correpsonding thresholded coefficient
using a threshold of \eqn{\delta}{delta}.
The supported shrinkage methods are

\describe{
\item{hard thresholding}{\deqn{ W^{\mbox{(t)}}=\left\{
\begin{array}{ll}
0,& \mbox{if $|W| \le \delta$;}\\
W,& \mbox{otherwise}
\end{array}
\right.}{W'=0 if |W| <= delta; W otherwise.}}

\item{soft thresholding}{\deqn{ W^{\mbox{(t)}}=\mbox{sign}(W)\;f( |W| - \delta ) \]
where
\[ \mbox{sign}(W) \equiv \left\{
\begin{array}{ll}
+1,& \mbox{if $W > 0$;}\\
0,& \mbox{if $W=0$;}\\
-1,& \mbox{if $W < 0$.}
\end{array}
\right.\]
and
\[
f(x) \equiv \left\{
\begin{array}{ll}
x,& \mbox{if $x \ge 0$;}\\
0,& \mbox{if $x < 0$.}
\end{array}
\right.}{W'=sign{W}*f(|W| - delta) where f(x)=x if x >= 0; 0 otherwise and sign(x)=+1 if x > 0; 0 if x=0; -1 if x < 0.}}

\item{mid thresholding}{\deqn{ W^{\mbox{(t)}}=\mbox{sign}(W)\;g( |W| - \delta ) \]
where
\[ g( |W| - \delta ) \equiv \left\{
\begin{array}{ll}
2 f(|W| - \delta),& \mbox{if $|W| < 2\delta$;}\\
|W|,& \mbox{otherwise.}
\end{array}
\right.}{W'=sign{W}*g(|W| - delta) where g(W)=2*f(|W|-delta) if |W| < 2*delta; |W| otherwise.}}}


\emph{Hard} thresholding reduces to zero all coefficients that do not
exceed the threshold. \emph{Soft} thresholding pushes toward zero any
coefficient whose magnitude exceeds the threshold, and zeros the
coefficient otherwise. \emph{Mid} thresholding represents a compromise
between hard and soft thresholding such that coefficients whose
magnitude exceeds twice the threshold are not adjusted, those between
the threshold and twice the trhreshold are shrunk, and those
below the threshodl are zeroed.

The threshold is selected based on a model of the noise. The supported
techniques for estimating the noise threshold are

\describe{
\item{universal}{\eqn{\sqrt{2\sigma^2_\varepsilon \log(N)}}{delta=sqrt( 2*var{noise}*log(N) )}
where  is the number of samples in the time series. As the noise variance
is typically unknown, it is estimated based on the median absolute deviation
of the absolute value of the scale one wavelet coefficients (and scaled by
dividing the result by \eqn{0.6745} so that if the series were Gaussian white noise, the
correct variance would be returned). The universal threshold is defined so that
if the original time series was solely comprised of Gaussian noise, then all the
wavelet coefficients would be (correctly) set to zero using a hard thresholding scheme.
Inasmuch, the universal threshold results in highly smoothed output.}

\item{minimax}{These thresholds are used with soft and hard thresholding,
and are precomputed based on a minimization of a theoretical upperbound
on the asymptotic risk. The minimax thresholds are always smaller than the
universal threshold for a given sample size, thus resulting in relatively
less smoothing.}

\item{adaptive}{These are scale-adaptive thresholds, based on the minimization of
Stein's Unbiased Risk Estimator for each level of the DWT. This method is only
available with soft shrinkage. As a caveat, this threshold can produce poor
results if the data is too sparse (see the references for details).}}


Finally, the user has the choice of using either a decimated (standard)
form of the discrete wavelet transform (DWT) or an undecimated version
of the DWT (known as the Maximal Overlap DWT (MODWT)).
Unlike the DWT, the MODWT is a (circular) shift-invariant transform so
that a circular shift in the original time series produces an equivalent shift of
the MODWT coefficients. In addition, the MODWT can be interpreted as
a \emph{cycle-spun} version of the DWT, which is achieved by averaging
over all non-redundant DWTs of shifted versions of the original series. The z
is a smoother version of the DWT at the cost of an increase in computational
complexity (for an N-point series, the DWT requires \eqn{O(N)} multiplications
while the MODWT requires \eqn{O(N\log_2N}{O(N log2(N))} multiplications.
}
\references{
Donoho, D. and Johnstone, I. \emph{Ideal Spatial Adaptation
by Wavelet Shrinkage}. Technical report, Department of Statistics,
Stanford University, 1992.

Donoho, D. and Johnstone, I. \emph{Adapting to Unknown Smoothness
via Wavelet Shrinkage}. Technical report, Department of Statistics,
Stanford University, 1992.

D. B. Percival and A. T. Walden, \emph{Wavelet Methods for
Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavDaubechies}}, \code{\link{wavDWT}}, \code{\link{wavMODWT}}.}

\examples{
## MODWT waveshrinking using various thresh.scale 
## values on sunspots series 
x <- as.vector(sunspots)
tt <- as.numeric(time(sunspots))
thresh <- seq(0.5,2,length=4)
ws <- lapply(thresh, function(k,x)
    wavShrink(x, wavelet="s8",
       shrink.fun="hard", thresh.fun="universal",
       thresh.scale=k, xform="modwt"), x=x)

stackPlot(x=tt, y=data.frame(x, ws),
    ylab=c("sunspots",thresh),
    xlab="Time")

## DWT waveshrinking using various threshold 
## functions 
threshfuns <- c("universal", "minimax", "adaptive")
ws <- lapply(threshfuns, function(k,x)
    wavShrink(x, wavelet="s8",
       thresh.fun=k, xform="dwt"), x=x)

stackPlot(x=tt, y=data.frame(x, ws),
ylab=c("original", threshfuns),
xlab="Normalized Time")
}
\keyword{nonlinear}

