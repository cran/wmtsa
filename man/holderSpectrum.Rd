%%  WARNING: This file was automatically generated from the associated 
%%  wav_wtmm.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_wtmm.mid file in the project MID directory. Once the 
%%  wav_wtmm.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the holderSpectrum function

\name{holderSpectrum}
\alias{holderSpectrum}
\title{The Holder spectrum of a time series}
\concept{Holder spectrumsingularity detection}
\usage{holderSpectrum(x, n.scale.min=3, fit=lmsreg)}
\description{Using a  tree, this function returns time localized  exponent
estimations for a given time series.}

\arguments{
\item{x}{an object of class \code{wavCWTTree}.}
\item{fit}{a linear regression function to use in fitting the resulting data. Default: \code{lmsreg}.}

\item{n.scale.min}{the minimum number of scales (points) that a given \emph{suitable}
branch segment must have before being considered as an admissible candidate for
 exponent estimation. Default: 3.}
}

\value{
a list containing the estimated  exponents, associated times and corresponding branch number.
}
\details{
%
Many real-world time series contain sharp dicontinuities (cusps) which can be
attributed to rapid changes in the observed system. These cusps are called
\emph{singularities} and their strength can be quantified via localized  exponents
as follows: Let \eqn{f(t)} be a continuous real-valued function containing a singularity at
time \eqn{t_0}{t0}. The   exponent \eqn{h(t_0)}{h(t0)} is defined
as the supremum (least upper bound) of all exponents \eqn{h} which satisfies the condition
\deqn{|f(t) - P_n(t - t_0)| \le C|t - t_0|^{h(t_0)},}{|f(t) - Pn(t - t0)| <= C|t - t0|^h(t0)}
where \eqn{P_n(t - t_0)}{Pn(t - t0)} is a polynomial of degree \eqn{n \le h(t_0)}{n <= h(t0)}
and \eqn{C} is a constant.
The collection of exponents for a given time series denotes the so-called  spectrum.
Mallat demonstrated that a cusp singularity at time \eqn{t_0}{t0} can be estimated via the CWT by noting that
the wavelet transform modulus maxima behave as \eqn{W_{a,t_0}(f) \propto |a|^{h(t_0)}}{W(a,t0) ~ |a|^h(t0)}
as the scale \eqn{a \rightarrow 0}{approaches 0}.

Thus, the strength of cusp singularities in a given time series can be quantified by
\describe{
\item{i}{Calculate the CWT of the time series.}

\item{ii}{Find the modulus maxima of the CWT (WTMM).}

\item{iii}{Link the WTMM into separate branches based (mainly) on their position in time to form a WTMM tree.}

\item{iv}{For each branch in the tree, perform an exponential fit of the WTMM over an
admissible range of scale and as the scale approaches zero. The resulting \emph{scaling exponent} is an estimate of
the local  exponent for the time series. The occurrence of the singularity in time is recorded
as the location in time where the WTMM converges as the scale nears zero.}}


In practice, the above technique can be unstable when applied to observational data due to
negative moment divergences and so-called \emph{outliers} which correspond to the end points
of sample singularities. One must also be very careful in selecting an appropriate scaling
region of a tree branch before fitting the data. We accomplish this by first segmenting a
given tree branch into regions which exhibit approximate linear behavior in the log(scale)-log(WTMM) space,
and subsequently selecting the region corresponding to the smallest scales for  exponent estimation.
Furthermore, through the \code{n.scale.min} argument, the user can control the minimum number of scales (points)
that must exist in the isolated scaling region before a
 exponent estimation is recorded.
}
\references{
S.G. Mallat, \emph{A Wavelet Tour of Signal Processing (2nd Edition)}, Academic Press, Cambridge, 1999.

S.G. Mallat and W.L. Hwang, ``Singularity detection and processing with wavelets",
\emph{IEEE Transactions on Information Theory}, \bold{38}, 617--643 (1992).

S.G. Mallat and S. Zhong, ``Complete signal representation with multiscale edges",
\emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{14}, 710--732 (1992).

J.F. Muzy, E. Bacry, and A. Arneodo, ``The multifractal formalism revisited with wavelets.",
\emph{International Journal of Bifurcation and Chaos}, \bold{4}, 245--302 (1994).
}

\seealso{
\code{\link{wavCWT}}, \code{\link{wavCWTFilters}}, \code{\link{wavCWTTree}}.}

\examples{
## create series with a linear trend and two 
## cusps: h(x = 1) = 0.5 and h(x = 15) = 0.3 
cusps <- function(x) -0.2 * abs(x-1)^0.5 - 0.5* abs(x-15)^0.3 + 0.00346 * x + 1.34
x <- seq(-5, 20, length=1000)
y <- signalSeries(cusps(x), x)

## calculate CWT using Mexican hat filter 
W <- wavCWT(y, wavelet="gaussian2")

## calculate WTMM and extract first two branches 
## in tree corresponding to the cusps 
W.tree <- wavCWTTree(W)[1:2]

## plot the CWT tree overlaid with a scaled 
## version of the time series to illustrate 
## alignment of branches with cusps 
yshift <- y@data - min(y@data)
yshift <- yshift / max(yshift) * 4 - 4.5
plot(W.tree, xlab="x")
lines(x, yshift, lwd=2)
text(6.5, -1, "f(x) = -0.2|x-1|^0.5 - 0.5|x-15|^0.3 + 0.00346x + 1.34", cex=0.8)

## estimate Holder exponents 
holder <- holderSpectrum(W.tree)
print(holder)
}
\keyword{univar}
\keyword{models}
\keyword{nonlinear}

