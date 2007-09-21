%%  WARNING: This file was automatically generated from the associated 
%%  wav_var.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_var.mid file in the project MID directory. Once the 
%%  wav_var.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the D.table function

\name{D.table}
\alias{D.table}
\title{Critical D-statistic table generation}
\concept{wavelet}
\usage{D.table(n.sample=c(127, 130), significance=c(0.1, 0.05, 0.01),
    lookup=TRUE, n.realization=10000, n.repetition=3,
    tolerance=1e-6)}
\description{The \code{D-statistic} denotes the maximum
deviation of sequence from a hypothetical linear cumulative energy
trend. The critical \code{D-statistics} define the distribution of D for a
zero mean Gaussian white noise process. Comparing the sequence
\code{D-statistic} to the corresponding critical
values provides a means of quantitatively rejecting or accepting the
linear cumulative energy hypothesis. The table is generated for an
ensemble of distribution probabilities and sample sizes.}

\arguments{

\item{lookup}{a logical flag for accessing precalculated
critical \code{D-statistics}. The
critical \code{D-statistics} are calculated for a variety of sample sizes
and significances.  If lookup is \code{TRUE} (recommended), this table is
accessed. The table is stored as the matrix object \code{D.table.critical}. Missing table values are calculated
using the input arguments: \code{n.realization}, \code{n.repetition}, and \code{tolerance}. Default: \code{TRUE}.}

\item{n.realization}{an integer specifying the number of realizations to generate in a
Monte Carlo simulation for calculating the \code{D-statistic(s)}. This
parameter is used either when lookup is \code{FALSE},
or when lookup is \code{TRUE} and the table is missing
values corresponding to the specified significances. Default: \code{10000}.}

\item{n.repetition}{an integer specifying the number of Monte Carlo simulations to
perform. This parameter coordinates with the n.realization
parameter. Default: \code{3}.}

\item{n.sample}{a vector of integers denoting the sample sizes for which critical
\code{D-statistics} are created. Default: \code{c(127,130)}.}

\item{significance}{a numeric vector of real values in the interval (0,1).
The significance is the fraction of times that the
linear cumulative energy hypothesis is incorrectly rejected. It is
equal to the difference of the distribution probability (p) and unity.
Default: \code{c(0.1, 0.05, 0.01)}.}

\item{tolerance}{a numeric real scalar that specifies the amplitude threshold to use in
estimating critical \code{D-statistic(s)} via the Inclan-Tiao approximation.
Setting this parameter to a higher value
results in a lesser number of summation terms at the expense of obtaining
a less accurate approximation. Default: \code{1e-6}.}
}

\value{
a matrix containing the critical
\code{D-statistics} corresponding to the supplied sample sizes and
significances.
}
\details{
%
A precalculated critical \code{D-statistics} object
(\code{D.table.critical}) exists
on the package workspace and was built for a variety of sample sizes and
significances using 3 repetitions and \code{10000}
realizations/repetition. This \code{D.table} function should be used in
cases where specific \code{D-statistic}s are missing from
\code{D.table.critical}.
Note: the results of the \code{D.table} value should not be returned to a
variable named \code{D.table.critical} as it will override the
precalculated table available in the package.

An Inclan-Tiao approximation of critical \code{D-statistics} is used for sample
sizes \code{n.sample} \eqn{\ge 128}{>= 128} while a Monte Carlo technique is used for
\code{n.sample} \eqn{< 128}{< 128}. For the
Monte Carlo technique, the \code{D-statistic} for a
Gaussian white noise sequence of length \code{n.sample} is calculated. This
process is repeated n.realization times, forming a distribution of the
\code{D-statistic}.  The critical values corresponding to the significances
are calculated a total of n.repetition times, and averaged to form
an approximation to the \code{D-statistic(s)}.
}
\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{D.table.critical}}.}

\examples{
D.lookup <- D.table(significance=c(10,5,1)/100,
    n.realization=100, n.sample=125:130, lookup=FALSE)
}
\keyword{distribution}

