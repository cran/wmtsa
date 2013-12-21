%%  WARNING: This file was automatically generated from the associated 
%%  wav_var.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_var.mid file in the project MID directory. Once the 
%%  wav_var.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the wavVarTest, plot.wavVarTest, print.wavVarTest functions

\name{wavVarTest}
\alias{wavVarTest}
\alias{plot.wavVarTest}
\alias{print.wavVarTest}
\title{Homogeneity test for discrete wavelet transform crystals}
\concept{waveletwavelet variance}
\usage{wavVarTest(x, wavelet="s8", n.levels=NULL,
    significance=c(0.1,0.05,0.01), lookup=TRUE, n.realization=10000,
    n.repetition=3, tolerance=1e-6)}
\description{Tests for homogeneity of variance for each
scale of a discrete wavelet transform (DWT) decomposition. Based on the
assumption that the DWT decorrelates colored noise processes, the
interior wavelet coefficients in a given scale (\eqn{{\mathbf d}_j^{int}}{dj}) can be regarded as a
zero mean Gaussian white noise process. For a homogeneous distribution
of \eqn{{\mathbf d}_j^{int}}{dj}, there is an expected linear increase in the cumulative energy
as a function of time. The so called D-statistic denotes the maximum
deviation of the \eqn{{\mathbf d}_j^{int}}{dj} from a hypothetical linear cumulative energy
trend. This D-statistic is then compared to a table of critical
D-statistics that defines the distribution of D for various sample
sizes. Comparing the D-statistic of \eqn{{\mathbf d}_j^{int}}{dj} to the corresponding critical
values provides a means of quantitatively rejecting or accepting the
linear cumulative energy hypothesis. This function performs
this test for an ensemble of distribution probabilities.}

\arguments{
\item{x}{an object of class \code{wavTransform} as output by the \code{wavDWT} function, a
corresponding \code{wavBoundary} object, or a
numeric vector. In the latter case, \code{wavDWT} parameters can be passed to specify the
type of wavelet to use and the number of decomposition levels to perform.}
\item{lookup}{a logical flag for accessing precalculated critical D-statistics. The
critical D-statistics are calculated for a variety of sample sizes
and significances.  If \code{lookup} is \code{TRUE}, this table is
accessed. The table is stored as the matrix object
\code{D.table.critical} and is loaded with the package.
Missing table values are calculated using the input arguments:
\code{n.realization, n.repetition}
and \code{tolerance}.
Default: \code{TRUE}.}

\item{n.levels}{the number of decomposition levels. Valid only for input not of class
\code{wavTransform} or \code{wavBoundary}.
Default: the maximum decomposition level that contains at least one interior wavelet coefficient.}

\item{n.realization}{an integer specifying the number of realizations to generate in a
Monte Carlo simulation for calculating the D-statistic(s). This
parameter is used either when \code{lookup}
is \code{FALSE},
or when \code{lookup}
is \code{TRUE} and the table is missing
values corresponding to the specified significances. Default: \code{10000}.}

\item{n.repetition}{an integer specifying the number of Monte Carlo simulations to
perform. This parameter is coordinated with the \code{n.realization}
parameter. Default: \code{3}.}

\item{significance}{a numeric vector of real values on the interval (0,1).
Qualitatively the significance is the fraction of times that the
linear cumulative energy hypothesis is incorrectly rejected. It is
equal to the difference of the distribution probability (p) and unity. Default: \code{c(0.1, 0.05, 0.01)}.}

\item{tolerance}{a numeric real scalar that specifies the amplitude threshold to use in
estimating critical D-statistic(s) via the Inclan-Tiao approximation.
Setting this parameter to a higher value
results in a lesser number of summation terms at the expense of obtaining
a less accurate approximation. Default: \code{1e-6}.}

\item{wavelet}{a character string denoting the filter type. Valid only for input not of class
\code{wavTransform} or \code{wavBoundary}.
Default: \code{"s8"}.}
}

\value{
an object of class \code{wavVarTest}.
}
\details{
%
An Inclan-Tiao approximation of critical D-statistics is used for sample
sizes \eqn{N \ge 128}{N >= 128} while a
Monte Carlo technique is used for
\eqn{N < 128}{N < 128}.
For the Monte Carlo technique, the D-statistic for a
Gaussian white noise sequence of length N is calculated. This
process is repeated \code{n.realization} times,
forming a distribution of the D-statistic.
The critical values corresponding to the significances
are calculated a total of \code{n.repetition}
times, and averaged to form
an approximation to the D-statistic(s).
Because the Monte Carlo study can be both computationally and memory
intensive, it is highly recommended that lookup be set to
\code{TRUE}, its default value.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavVar}}, \code{\link{wavDWT}}, \code{\link{D.table}}.}

\examples{
## perform a homogeneity of variance test for a 
## DWT decomposition of a long memory process 
## realization 
homogeneity <- wavVarTest(fdp045)
}
\keyword{univar}

