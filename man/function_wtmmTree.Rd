%%  WARNING: This file was automatically generated from the associated 
%%  wav_wtmm.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_wtmm.mid file in the project MID directory. Once the 
%%  wav_wtmm.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wtmmTree, [.wtmmTree, plot.wtmmTree functions

\name{wtmmTree}
\alias{wtmmTree}
\alias{[.wtmmTree}
\alias{plot.wtmmTree}
\title{Tree map of continuous wavelet transform modulus maxima}
\concept{singularity detectionwavelet}
\usage{wtmmTree(x, bridge.gaps=FALSE, strength.min=0.,
    n.octave.min=2, tolerance=0.0, border=FALSE)}
\description{This function first finds the modulus maxima locations (in time and in scale)
of the continuous wavelet transform input. The set of modulus maxima are then
segmented into \emph{branches}, where each branch represents a collection of WTMM that
correspond to the same ridge in the WTMM time-scale plane.
A coarse-to-fine scale strategy is used to identify the members
of each branch as follows: (i) a single WTMM at the coarsest
scale is selected as the start of a given branch,
(ii) the closest neighboring WTMM in time at the next finest
scale is then added to the branch, (iii) step ii is repeated
until the smallest scale is reached or an apparent break
occurs in the branch across scale, and (iv) steps i-iii are
repeated until all WTMM have been accounted. Branches are
allowed to wrap off one end of the time-scale plane in time
and back onto the other. A branch is not grown unless
the nearest neighbor candidate at the next finest scale is
close in time to the last recorded branch member, where "close"
is defined as being less than the current scale of the neighbor
candidate. This means that the window in time for admissible
neighbor WTMM candidates (at the next finest scale) shrinks
proportionally with scale.}

\arguments{
\item{x}{an object of class \code{wavCWT} (as produced by the \code{wavCWT} function).}
\item{border}{a logical flag. If \code{border} is \code{TRUE} WTMM branches which fall off one end
of the CWT branch in time and return on the other are allowed. Otherwise,
these border branches are pruned from the returned list. Default: \code{FALSE}.}

\item{bridge.gaps}{a logical flag. Gaps encountered in the search for
additional branch members may be bridged by setting the
\code{bridge.gaps} boolean to \code{TRUE}. Default: \code{FALSE}.}

\item{n.octave.min}{a pruning factor for excluding non-persistent branches. If
a WTMM branch does not span this number of octaves, it is
excluded from the tree. Default: 2.}

\item{strength.min}{a pruning factor for excluding weak branches. A given WTMM
is excluded from the current branch if it less than or equal
to the product of the maximum WTMM at the current scale and
\code{strength.min}. This parameter must be in the range [0,1].
If set to zero, all WTMM are allowed. If set to unity,
only the maximum WTMM at the current scale is allowed.
See the \code{DETAILS} section for more details.
Default: 0.}

\item{tolerance}{a tolerance vector used to find moduls maxima in the given CWT.
This vector must be as long as there are scales in the CWT such that
the \eqn{j^{th}}{jth} element defines the
tolerance to use in finding modulus maxima at
the \eqn{j^{th}}{jth} scale of the CWT.
Default: 0.}
}

\value{
an object of class \code{WaveletWTMMTree}.
}

\section{S3 METHODS}{
\describe{

\item{[}{extracts a subset of branches from the tree. For example, to extract branches 2 through 5,
use the syntax \code{x[2:5]}. To extract branches which terminate near times 0.47, 0.3, and 1.4, use the syntax
\code{x[time=c(0.47, 0.3, 1.4)]}. To extract all branches which terminate between times 1.2 and 1.5, use the syntax
\code{x[range=c(1.2, 1.5)]}.}

\item{plot}{plots the WTMM tree. The plot method also supports the following optional arguments
(assume that the variable \code{x} is an output of the \code{wtmmTree} function):
\describe{
\item{fit}{A logical flag. If \code{fit=TRUE}, a subset of branches (limited to four) are fit with
various linear regression models on a log(|WTMM|) versus log(scale) basis. The models are specified by the
optional \code{models} argument. This scheme illustrates the process by which 
exponents are estimated using the WTMM branches. For example, to see the regressions over chains 10 through 13,
issue \code{plot(x[10:13], fit=TRUE)}. Default: \code{FALSE}.}

\item{models}{A vector of character strings denoting the linear models to use in illustring the calculation
of  exponents. This argumetn is used only if \code{fit=TRUE}. Default: \code{c("lm", "lmsreg", "ltsreg")}.}

\item{labels}{Logical flag. If \code{TRUE}, the branch number is placed at the head of each branch. Default: \code{TRUE}.}

\item{wtmm}{A logical flag. If \code{TRUE}, all of the (non-pruned and unbranched) WTMM are
plotted in the time-scale plane. Default: \code{FALSE}.}

\item{pch}{The marker used in plotting branch points via the \code{par} function. Default: \code{"o"}.}}}

}
}

\details{
%
A point in the CWT \eqn{W(t,j)} is defined
as a maximum if \eqn{|W(t-1,j)| + tol < |W(t,j)|} and \eqn{|W(t+1,j)| + tol < |W(t,j)|}
where \code{tol} is a (scale-dependent) tolerance specified by the user.
The search algorithm is also adpated to identify plateaus in the data,
and will select the the middle of the plateau as a maximum location
when encountered. The data \eqn{|W(t,j)|} is first scaled so that its
maximum value is 1.0, so the tolerances should be adjusted accordingly.
Since the CWT coefficients are in effect a result band-pass filtering operations,
the large scale coefficients form a smoother curve than do the small
scale coefficients. Thus, the tolerance vector allows the user to specify
scale-dependent tolerances, helping to weed out undesirable local maxima.
It is recommended that the tolerance be set proportional to the scale,
e.g., \code{tolerance=C / sqrt(scale)} where \eqn{C} is a constant
\eqn{0 < C < 1}{0 < C < 1}.
}
\references{
J.F. Muzy, E. Bacry, and A. Arneodo., ``The multifractal formalism revisited with wavelets.",
\emph{International Journal of Bifurcation and Chaos}, \bold{4}, 245--302 (1994).
}

\seealso{
\code{\link{wavCWT}},
\code{\link{wavCWTFilters}}.
}

\examples{
## calculate the CWT of a random walk series and 
## the plot the result with an overlay of the 
## original time series above the CWT image 
set.seed(100)
x <- cumsum(rnorm(1024))
x.cwt <- wavCWT(x)
plot(x.cwt, series=TRUE)

## create a WTMM tree and plot the first 50 
## branches found 
x.tree <- wtmmTree(x.cwt)
plot(x.tree[1:50])

## plot an illustration of the Holder exponent 
## estimation process. select branches between 
## times 100 and 200 (only the first four found 
## will be fitted) 
plot(x.tree[range=c(100, 200)], fit=TRUE)
}
\keyword{univar}
\keyword{models}
\keyword{nonlinear}

