%%  WARNING: This file was automatically generated from the associated 
%%  wav_var.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the wav_var.mid file in the project MID directory. Once the 
%%  wav_var.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the D.statistic function

\name{D.statistic}
\alias{D.statistic}
\title{D-statistic}
\concept{waveletwavelet variance}
\usage{D.statistic(x)}
\description{Given a numeric vector \code{x},
this function calculates the
the maximum departure of \code{x} from an expected
linear increase in cumulative energy based on
a white noise model.}

\arguments{
\item{x}{a numeric vector.}

}
\seealso{
\code{\link{wavVar}}, \code{\link{wavVarTest}}, \code{\link{D.table}}, \code{\link{D.table.critical}}.}

\examples{
## compare the D-statistic for a white noise 
## realization and a random walk. the random 
## walk D-statistic will be relatively large in 
## comparison to that of the white noise 
## realization, signifying a stronger departure 
## from an expected increase in cumulative 
## energy. 
D.statistic(rnorm(1024))
D.statistic(cumsum(rnorm(1024)))
}
\keyword{univar}

