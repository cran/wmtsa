%%  WARNING: This file was automatically generated from the associated 
%%  wav_wtmm.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_wtmm.mid file in the project MID directory. Once the 
%%  wav_wtmm.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavCWTPeaks function

\name{wavCWTPeaks}
\alias{wavCWTPeaks}
\title{Peak detection in a time series via the CWT}
\concept{peak detectionwavelet}
\usage{wavCWTPeaks(x, snr.min=3, scale.range=NULL, length.min=10,
    noise.span=NULL, noise.fun="quantile", noise.min=NULL)}
\description{Finds the local maxima in a time series via a CWT tree.}

\arguments{
\item{x}{an object of class \code{wavCWTTree}.}
\item{length.min}{the minimum number of points along a CWT tree branch and within the specified \code{scale.range} needed in order for that branches peak to be considered a peak candidate.
Default: \code{10}.}

\item{noise.fun}{a character string defining the function to apply to the local noise estimates in order to sumarize
and quantify the local noise level into a scalar value. See the \bold{DETAILS} section for more information. Supported values
are
\describe{
\item{"quantile"}{quantile(x, probs=0.95)}

\item{"sd"}{sd(x)}

\item{"mad"}{mad(x, center=0)}}


where \code{x} is a vector of smallest-scale CWT coefficients whose time indices are near that of the branch termination time.
Default: \code{"quantile"}.}

\item{noise.min}{the minimum allowed estimated local noise level. Default: \code{quantile(attr(x,"noise"), prob=0.05)},
where \code{x} is the input \code{wavCWTTree} object.}

\item{noise.span}{the span in time surrounding each branche's temrination point to use in forming
local noise estimates and (ultimately) peak SNR estimates. Default: \code{NULL,max(0.01 * diff(range(times)), 5*sampling.interval)},
where \code{times} and \code{sampling.interval} are attributes of the input \code{wavCWTTree} object.}

\item{scale.range}{the range of CWT scales that a peak must fall into in order to be considered a peak candidate.
Default: \code{scale[range(which(branch.hist > quantile(branch.hist,prob=0.8)))]}, where \code{branch.hist} is an
attribute of the input \code{wavCWTTree} object. This default results in isolating the bulk of energetic CWT branches, but
the user in encouraged to reduce the scale range in order to attenuate the computational burden of the peak detection scheme.}

\item{snr.min}{the minimum allowed peak signal-to-noise ratio. Default: \code{3}.}
}

\value{
a list of \code{x} and \code{y} vectors identifying the peaks in the original time series.
The pruning criteria (\code{snr.min}, \code{scale.range}, \code{length.min}, \code{noise.span}, \code{noise.fun}, \code{noise.min})
are attached are attached as attributes. In addition, a \code{peaks} attribute is
attached and corresponds to a \code{data.frame} containing the following information for each peak:

\item{branch}{index of the associated branch in the CWT tree}

\item{itime}{index location in time}

\item{iscale}{index location in scale}

\item{time}{location in time}

\item{scale}{location in scale}

\item{extrema}{CWT value}

\item{iendtime}{index location of branch termination time, i.e., the index of the point in the time series corresponding to the current peak}
}
\details{
%
The local maxima of the CWT are linked together
to form so-called branches, where each branch represents one \emph{ridge} of the
CWT time-scale terrain. The collection of branches forms a tree, as output
by the \code{\link{wavCWTTree}} function. The \code{wavCWTpeaks} function prunes the branches of the
input CWT tree and records the termination time (i.e., the time associated with
point of the branch that is closest to scale zero) as the time index associated
with the local peak of the corresponding time series. Information regarding the collection of isolated peaks
is returned as a \code{\link[base]{data.frame}} object.

The tree branches are pruned in the following ways:
\describe{
\item{peak SNR}{an estimate of SNR at peak value is greater than or equal to the specified \code{snr.min}.
A peak SNR estimate is formed as follows:
For each branch of the input CWT tree, a subset of
CWT coefficients is collected such that the CWT coefficients are both local to the branch termination time and correspond
to the smallest analyzed CWT scale. The user specified \code{noise.span} argument is used to define the boundaries of
each subset in time ala \eqn{[B - }\code{noise.span}, \eqn{B + }\code{noise.span}\eqn{]}, where \eqn{B} is the branch termination
time. Each CWT subset is assumed to be representative of the local noise levels near the corresponding branch termination time
and \code{noise.fun} is used to quantify (and summarize) each level resulting in a scalar $z$. The minimum value of $z$ is
specified by the user ala the \code{noise.min} argument. Finally, the ratio \eqn{|P|/|W|} is used to form an estimate of
the local signal-to-noise ration (SNR) for the corresponding branch, where \eqn{P} is the maximum CWT value along the
branch in the CWT time-scale plane.}

\item{scale}{the scale corresponding to the peak is larger than the minimum of the specified \code{scale.range}.}

\item{branch length}{the length of the branch within the specified \code{scale.range} is greater than or equal to the specified minimum \code{length.min}.}

\item{endpoint}{the index of the terminating time of the branch is on the interval \eqn{(W, N-W)}, where \eqn{N} is the length of the
series and \eqn{W} is integer equivalent of \eqn{1/4} the length of the \code{noise.span} or \eqn{3}, whichever is greater.}}


\bold{NOTE:} For peak detection, the wavelet filters used to form the CWT must maintain an (approximate) zero phase
property so that the CWT coefficients can be meaningfully aligned with the events of the original
time series. Currently, only the so-called Mexican hat wavelet maintains this property due to
the even-symmetry of the filter's impulse response. Therefore, only the Mexican hat wavelet ("gaussian2")
is currently supported for CWT-based peak detection. See the \code{\link{wavCWTFilters}} and \code{\link{wavCWT}}
function for more information.
}
\references{
Pan Du, Warren A. Kibbe, and Simon M. Lin,
``Improved peak detection in mass spectrum by incorporating continuous wavelet transform-based pattern matching",
\emph{Bioinformatics}, \bold{22}, 2059--2065 (2006).

J.F. Muzy, E. Bacry, and A. Arneodo., ``The multifractal formalism revisited with wavelets.",
\emph{International Journal of Bifurcation and Chaos}, \bold{4}, 245--302 (1994).
}

\seealso{
\code{\link{wavCWTTree}}, \code{\link{wavCWT}}, \code{\link{wavCWTFilters}}.}

\examples{
## create linchirp series 
linchirp <- make.signal("linchirp")

## calculate the CWT 
W <- wavCWT(linchirp)

## form CWT tree 
z <- wavCWTTree(W)

## estimate the peak locations using default 
## scale.range 
p <- wavCWTPeaks(z)

## plot an overlay of the original series and the 
## peaks 
x <- as(linchirp@positions,"numeric")
y <- linchirp@data
plot(x, y, type="l", xlab="time", ylab="linchirp")
points(p, pch=16, col="red", cex=1.2)
}
\keyword{univar}
\keyword{models}
\keyword{nonlinear}

