\name{recover.B}
\alias{recover.B}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Recover an estimate of the treatment effect at time t in Study B
}
\description{
Recover an estimate of the treatment effect at time t in Study B
}
\usage{
recover.B(Axzero, Adeltazero, Aszero, Axone, Adeltaone, Asone, Bxzero, Bdeltazero, 
Bszero, Bxone, Bdeltaone, Bsone, t, landmark, extrapolate = T, transform = F)
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
  \item{Axone}{
%%     ~~Describe \code{Axone} here~~
observed event times in the treatment group in  Study A
}
  \item{Adeltaone}{
%%     ~~Describe \code{Adeltaone} here~~
event/censoring indicators in the treatment group in Study A
}
  \item{Asone}{
%%     ~~Describe \code{Asone} here~~
surrogate marker values in the treatment group in Study A, NA for individuals not observable at the time the surrogate marker was measured
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
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether a transformation should be used, default is FALSE
}
}
\details{
Assume there are two randomized studies of a treatment effect, a prior study (Study A) and a current study (Study B). Study A was completed up to some time t, while Study B was stopped at time \eqn{t_0<t}. In both studies, a surrogate marker was measured at time \eqn{t_0} for individuals still observable at \eqn{t_0}. Let \eqn{G} be the binary treatment indicator with \eqn{G=1} for treatment and \eqn{G=0} for control and we assume throughout that subjects are randomly assigned to a treatment group at baseline. Let \eqn{T_K^{(1)}} and \eqn{T_K^{(0)}} denote the time of the primary outcome of interest, death for example, under the treatment and under the control, respectively, in Study K. Let \eqn{S_K^{(1)}} and \eqn{S_K^{(0)}} denote the surrogate marker measured at time \eqn{t_0} under the treatment and the control, respectively, in Study K.

The treatment effect quantity of interest, \eqn{\Delta_K(t)}, is the difference in survival rates by time \eqn{t} under treatment versus control,
  \deqn{ \Delta_K(t)=E\{ I(T_K^{(1)}>t)\} - E\{I(T_K^{(0)}>t)\} = P(T_K^{(1)}>t) - P(T_K^{(0)}>t)}
where \eqn{t>t_0}. Here, we recover an estimate of \eqn{\Delta_B(t)} using Study B information (which stopped follow-up at time \eqn{t_0<t}) and Study A information (which has follow-up information through time t). The estimate is obtained as \deqn{\hat{\Delta}_{EB}(t,t_0)/ \hat{R}_{SA}(t,t_0)}
where \eqn{\hat{\Delta}_{EB}(t,t_0)} is the early treatment effect estimate in Study B, described in the early.delta.test documention, and \eqn{\hat{R}_{SA}(t,t_0)} is the proportion of treatment effect explained by the surrogate marker information at \eqn{t_0} in Study A. This proportion is calculated as \eqn{\hat{R}_{SA}(t,t_0) =\hat{\Delta}_{EA}(t,t_0)/\hat{\Delta}_A(t)}
where \deqn{\hat{\Delta}_A(t)=n_{A1}^{-1}\sum_{i=1}^{n_{A1}}\frac{I(X_{Ai}^{(1)}>t)}{\hat{W}_{A1}^C(t)}-n_{A0}^{-1}\sum_{i=1}^{n_{A0}}\frac{I(X_{Ai}^{(0)}>t)}{\hat{W}_{A0}^C(t)},}
and \eqn{\hat{\Delta}_{EA}(t,t_0)} is parallel to \eqn{\hat{\Delta}_{EB}(t,t_0)} except replacing
\eqn{n_{A0}^{-1} \sum_{i=1}^{n_{A0}} \hat{r}_A^{(0)}(t|S_{Ai}^{(0)}, t_0) \frac{I(X_{Ai}^{(0)} > t_0)}{\hat{W}_{A0}^C(t_0)}}
by \eqn{n_{A0}^{-1}\sum_{i=1}^{n_{A0}}\hat{W}_{A0}^C(t)^{-1}I(X_{Ai}^{(0)}>t),} and \eqn{\hat{W}^C_{Ag}(\cdot)} is the Kaplan-Meier estimator of the survival function for \eqn{C_{A}^{(g)}} for \eqn{g=0,1}.  

Perturbation resampling is used to provide a standard error estimate for the estimate of \eqn{\Delta_B(t)} and a confidence interval.

}
\value{
\item{recovered.deltaB}{The recovered estimate of \eqn{\Delta_B(t)}.}
\item{sd.recovered.deltaB}{The standard error estimate of the recovered estimate of \eqn{\Delta_B(t)}.}
\item{conf.quantile.recovered.deltaB}{A confidence interval for the recovered estimate of \eqn{\Delta_B(t)}.}
}
\references{
Parast L, Cai T, Tian L (2019). Using a Surrogate Marker for Early Testing of a Treatment Effect. Biometrics, In press.

Parast L, Cai T and Tian L (2017). Evaluating Surrogate Marker Information using Censored Data. Statistics in Medicine, 36(11): 1767-1782.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
\donttest{
data(dataA)
data(dataB)
recover.B(Axzero = dataA$x0, Adeltazero = dataA$delta0, Aszero = dataA$s0, Axone 
= dataA$x1, Adeltaone = dataA$delta1, Asone = dataA$s1, Bxzero = dataB$x0, Bdeltazero
= dataB$delta0, Bszero = dataB$s0, Bxone = dataB$x1, Bdeltaone = dataB$delta1, Bsone 
= dataB$s1, t=1, landmark=0.5,  extrapolate = T)

recover.B(Axzero = dataA$x0, Adeltazero = dataA$delta0, Aszero = dataA$s0, Axone 
= dataA$x1, Adeltaone = dataA$delta1, Asone = dataA$s1, Bxzero = dataB$x0, Bdeltazero
= dataB$delta0, Bszero = dataB$s0, Bxone = dataB$x1, Bdeltaone = dataB$delta1, Bsone 
= dataB$s1, t=0.75, landmark=0.5,  extrapolate = T)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% use one of  RShowDoc("KEYWORDS")
\keyword{nonparametric}% __ONLY ONE__ keyword per line
