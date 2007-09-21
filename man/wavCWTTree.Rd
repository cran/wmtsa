%%  WARNING: This file was automatically generated from the associated 
%%  wav_wtmm.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_wtmm.mid file in the project MID directory. Once the 
%%  wav_wtmm.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavCWTTree, [.wavCWTTree, plot.wavCWTTree, print.summary.wavCWTTree, print.wavCWTTree, summary.wavCWTTree functions

\name{wavCWTTree}
\alias{wavCWTTree}
\alias{[.wavCWTTree}
\alias{plot.wavCWTTree}
\alias{print.summary.wavCWTTree}
\alias{print.wavCWTTree}
\alias{summary.wavCWTTree}
\title{Tree map of continuous wavelet transform extrema}
\concept{singularity detectionwavelet}
\usage{wavCWTTree(x, n.octave.min=1, tolerance=0.0, type="maxima")}
\description{This function first finds the extrema locations (in time and in scale)
of the continuous wavelet transform input. The set of extrema are then
subdivided into sets of \emph{branches}, where each branch represents
a collection of extrema that correspond to the same ridge in the CWT time-scale plane.
A coarse-to-fine scale strategy is used to identify the members
of each branch as follows: (i) a single extremum at the coarsest
scale is selected as the start of a given branch,
(ii) the closest neighboring extremum in time at the next finest
scale is then added to the branch, (iii) step ii is repeated
until the smallest scale is reached or an apparent break
occurs in the branch across scale, and (iv) steps i-iii are
repeated until all extrema have been accounted. A branch is not grown unless
the nearest neighbor candidate at the next finest scale is
close in time to the last recorded branch member, where "close"
is defined as being less than the current scale of the neighbor
candidate. This means that the window in time for admissible
neighbor extrema candidates (at the next finest scale) shrinks
proportionally with scale.}

\arguments{
\item{x}{an object of class \code{wavCWT} (as produced by the \code{wavCWT} function).}
\item{n.octave.min}{a pruning factor for excluding non-persistent branches. If
a branch of connected extrema does not span this number of octaves, it is
excluded from the tree. Default: 1.}

\item{tolerance}{a tolerance vector used to find CWT extrema.
This vector must be as long as there are scales in the CWT such that
the \eqn{j^{th}}{jth} element defines the
tolerance to use in finding modulus maxima at
the \eqn{j^{th}}{jth} scale of the CWT. If not, the last value
is replicated appropriately. Default: 0.}

\item{type}{a character string denoting the type of extrema to seek in the CWT plane.
Supported types are \code{"extrema"}, \code{"maxima"} and \code{"minima"}. Default: \code{"maxima"}.}
}

\value{
an object of class \code{wavCWTTree}. See \bold{DETAILS} section for more information.
}

\section{S3 METHODS}{
\describe{

\item{[}{extracts a subset of branches from the tree. For example, to extract branches 2 through 5,
use the syntax \code{x[2:5]}. To extract branches which terminate near times 0.47, 0.3, and 1.4, use the syntax
\code{x[time=c(0.47, 0.3, 1.4)]}. To extract all branches which terminate between times 1.2 and 1.5, use the syntax
\code{x[range=c(1.2, 1.5)]}.}

\item{plot}{plots the WTMM tree. The plot method also supports the following optional arguments
(assume that the variable \code{x} is an output of the \code{wavCWTTree} function):
\describe{
\item{fit}{A logical flag. If \code{fit=TRUE}, a subset of branches (limited to four) are fit with
various linear regression models on a log(|WTMM|) versus log(scale) basis. The models are specified by the
optional \code{models} argument. This scheme illustrates the process by which 
exponents are estimated using the WTMM branches. For example, to see the regressions over chains 10 through 13,
issue \code{plot(x[10:13], fit=TRUE)}. Default: \code{FALSE}.}

\item{models}{A vector of character strings denoting the linear models to use in illustring the calculation
of  exponents. This argumetn is used only if \code{fit=TRUE}. Default: \code{c("lm", "lmsreg", "ltsreg")}.}

\item{labels}{Logical flag. If \code{TRUE}, the branch number is placed at the head of each branch. Default: \code{TRUE}.}

\item{extrema}{A logical flag. If \code{TRUE}, the locations of the (non-pruned and unbranched) extrema are
marked in the time-scale plane. Default: \code{FALSE}.}

\item{pch}{The marker used in plotting branch points via the \code{par} function. Default: \code{"o"}.}}}

\item{print}{prints a summary of the object.}

\item{summary}{sumamrizes the branches comprising the tree ala a data.frame object}

}
}

\details{
%
A point in the CWT \eqn{W(t,j)} is defined
as an extremum if \eqn{|W(t-1,j)| + tol < |W(t,j)|} and \eqn{|W(t+1,j)| + tol < |W(t,j)|}
where \code{tol} is a (scale-dependent) tolerance specified by the user.
The search algorithm is also adpated to identify plateaus in the data,
and will select the the middle of the plateau as a maximum location
when encountered. The data \eqn{|W(t,j)|} is first scaled so that its
maximum value is 1.0, so the tolerances should be adjusted accordingly.
Since the CWT coefficients are (in effect) a result band-pass filtering operations,
the large scale coefficients form a smoother curve than do the small
scale coefficients. Thus, the tolerance vector allows the user to specify
scale-dependent tolerances, helping to weed out undesirable local maxima.
It is recommended that the tolerance be set proportional to the scale,
e.g., \code{tolerance=C / sqrt(scale)} where \eqn{C} is a constant
\eqn{0 < C < 1}{0 < C < 1}.
The user is also allowed to control the types of peaks to pursue in
the CWT plane: extrema, maxima, or minima. The algorithm (described above) is
adjusted accordingly.

The output object contains a list
of sublists, each sublist corresponds to a single branch in the CWT tree and contains the named vectors:
\describe{
\item{itime}{index locations in time of CWT extrema}

\item{iscale}{index locations in scale of CWT extrema}

\item{time}{times associated with CWT extrema}

\item{scale}{scales associated with CWT extrema}

\item{extrema}{CWT extrema values}}


In addition, the returned object contains the following attributes:

\describe{
\item{iendtime}{integer vector of indices corresponding to the locations in time where the branches terminated as the scale appraoches zero.}

\item{endtime}{numeric vector containing branch termination times}

\item{time}{numeric vector of times corresponding to the original time series}

\item{scale}{numeric vector of scales used to form the CWT}

\item{extrema.mask}{a binary matrix (of the same dimension as the CWT matrix) containing a 1 where there exists a corresponding extremum value in the CWT plane}

\item{noise}{a numeric vector containing the first scale's CWT coefficients. Statistical analysis of these data are often used as a rough estimate of the (local) noise level(s) in the original time series.}

\item{branch.hist}{a numeric vector containing the sum across time of all extrema values. This data can be used to help automate the selection of scales of interest in the CWT plane.}}
}
\references{
J.F. Muzy, E. Bacry, and A. Arneodo., ``The multifractal formalism revisited with wavelets.",
\emph{International Journal of Bifurcation and Chaos}, \bold{4}, 245--302 (1994).
}

\seealso{
\code{\link{wavCWT}}, \code{\link{wavCWTFilters}}.}

\examples{
## create linchirp series 
linchirp <- make.signal("linchirp")

## calculate the CWT 
W <- wavCWT(linchirp)

## form CWT tree 
W.tree <- wavCWTTree(W)

## print the object 
print(W.tree)

## summarize the object 
summary(W.tree)

## plot thea CWT image with a tree overlay 
## (R-only) 
plot(W)
if (is.R()) plot(W.tree, extrema=TRUE, add=TRUE)

## plot all CWT tree branches 
plot(W.tree)

## plot a subset of CWT tree branches 
plot(W.tree[5:10])

## plot an illustration of the Holder exponent 
## estimation process. select branches between 
## times 0.2 and 0.4 (only the first four found 
## will be fitted) 
plot(W.tree[range=c(0.2, 0.4)], fit=TRUE)
}
\keyword{univar}
\keyword{models}
\keyword{nonlinear}

