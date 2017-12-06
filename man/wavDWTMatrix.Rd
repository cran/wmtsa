\name{wavDWTMatrix}
\alias{wavDWTMatrix}
\title{
Discrete wavelet transform matrix.
}

\description{
Returns the discrete wavelet transform matrix for the specified wavelet and 
decomposition level.
}

\usage{
wavDWTMatrix(wavelet="d4", J=4, J0=J)
}

\arguments{
  \item{wavelet}{character string defining Daubechies wavelet filter.}
  \item{J}{integer scalar defining maximum decomposition level.}
  \item{J0}{integer scaling defining the partial decomposition level. Note
  that necessarily J0 <= J, otherwise an error will be thrown.}
}

\details{
The discrete wavelet transform of a real-valued uniformly-sampled time series
\eqn{X_t}{X(t)} for t=0,...,N-1 can be written in matrix form as 
\eqn{\mathbf{W} = \mathcal{W}\mathbf{X}}{W = BX}, where \eqn{\mathbf{W}}{W} are
the N x 1 transform coefficients and \eqn{\mathcal{W}}{B} is an orthonormal 
matrix of size N x N known as the discrete wavelet transform matrix. 
Generally, the discrete wavelet transform is calculated
via a computationally efficient pyramid algorithm and not via the matrix
form above. However, it serves as a good illustration and is primarily used 
for pedagogical interest. We highly recommend using the equivalent and highly
efficient alternative functions in this package for computing the discrete
wavelet transform such as \code{\link{wavDWT}}.

Note that the length of the time series is taken to be N=2^J. You can obtain the
discrete wavelet transform matrix for a partial transform by setting J0 < J. By
default, J0=J.
}

\value{
discrete wavelet transform matrix
}

\references{
D. B. Percival and  A. T. Walden,
\emph{Wavelet Methods for Time Series Analysis}, Cambridge University Press, 2000.
}

\seealso{
\code{\link{wavDaubechies}}, \code{\link{wavDWT}}.
}

\examples{
# Plot some rows of the wavelet transform matrix for the Haar filter
B <- wavDWTMatrix(wavelet = 'haar', J = 4)
wavStemPlot(B[1,])
wavStemPlot(B[15,])
wavStemPlot(B[16,])

# Compare the wavelet transform coefficients via the efficient
# wavDWT and non-efficient matrix multiplaction methods
N <- 256
X <- rnorm(N)
J0 <- 3
J <- as.integer(log2(N))
wavelet <- 's8'
W <- wavDWT(x = X, wavelet = wavelet, n.levels = J0)
B <- wavDWTMatrix(wavelet = wavelet, J = J, J0 = J0)
all.equal.numeric(as.vector(unlist(W$data)), as.vector(B \%*\% X))
}