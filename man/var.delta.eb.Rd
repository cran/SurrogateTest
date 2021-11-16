\name{var.delta.eb}
\alias{var.delta.eb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Variance estimation, generally not called directly by the user
}
\usage{
var.delta.eb(Axzero, Adeltazero, Aszero, Bxone, Bdeltaone, 
Bsone, Bxzero, Bdeltazero, Bszero, t, landmark = landmark, 
extrapolate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Axzero}{
%%     ~~Describe \code{Axzero} here~~
observed event times in the control group in  Study A
}
  \item{Adeltazero}{
%%     ~~Describe \code{Adeltazero} here~~
event/censoring indicators in the control group in Study A
}
  \item{Aszero}{
%%     ~~Describe \code{Aszero} here~~
surrogate marker values in the control group in Study A, NA for individuals not observable at the time the surrogate marker was measured
}
  \item{Bxone}{
%%     ~~Describe \code{Bxone} here~~
observed event times in the treatment group in  Study B
}
  \item{Bdeltaone}{
%%     ~~Describe \code{Bdeltaone} here~~
event/censoring indicators in the treatment group in  Study B
}
  \item{Bsone}{
%%     ~~Describe \code{Bsone} here~~
surrogate marker values in the treatment group in Study B, NA for individuals not observable at the time the surrogate marker was measured
}
  \item{Bxzero}{
%%     ~~Describe \code{Bxzero} here~~
observed event times in the control group in  Study B
}
  \item{Bdeltazero}{
%%     ~~Describe \code{Bdeltazero} here~~
event/censoring indicators in the control group in Study B
}
  \item{Bszero}{
%%     ~~Describe \code{Bszero} here~~
surrogate marker values in the control group in Study B, NA for individuals not observable at the time the surrogate marker was measured
}
   \item{t}{
%%     ~~Describe \code{t} here~~
time of interest
}
  \item{landmark}{
%%     ~~Describe \code{landmark} here~~
landmark time of interest, t0
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether local constant extrapolation should be used, default is TRUE
}
}
\details{
Variance estimation using the closed form expression under the null hypothesis of no treatment effect; more details are included in the documentation for early.delta.test.
}
\value{
Variance estimate for \eqn{\sqrt{n_{B}}\hat{\Delta}_{EB}(t,t_0)}
}
\references{
Parast L, Cai T, Tian L (2019). Using a Surrogate Marker for Early Testing of a Treatment Effect. Biometrics, 75(4):1253-1263.
}
\author{
Layla Parast
}
\examples{
data(dataA)
data(dataB)

var.delta.eb(Axzero = dataA$x0, Adeltazero = dataA$delta0, 
Aszero = dataA$s0, Bxone = dataB$x1, Bdeltaone = dataB$delta1, 
Bsone = dataB$s1, Bxzero = dataB$x0, Bdeltazero = dataB$delta0, 
Bszero = dataB$s0, t=1, landmark=0.5,  extrapolate = TRUE)

}

\keyword{internal}
