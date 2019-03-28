\name{design.study}
\alias{design.study}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power and sample size calculation for designing a future study
}
\description{
Power and sample size calculation for designing a future study
}
\usage{
design.study(Axzero, Adeltazero, Aszero, Axone = NULL, Adeltaone = NULL, Asone = 
NULL, delta.ea = NULL, psi = NULL, R.A.given = NULL, t, landmark, extrapolate = T, 
adjustment = F, n = NULL, power = NULL, pi.1 = 0.5, pi.0 = 0.5, cens.rate, transform = F)
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
observed event times in the treatment group in  Study A; optional (user must provide either (1) data from treatment arm in Study A or (2) hypothesized values for delta.ea (or R.A.given)and psi or (3) data from treatment arm in Study A and hypothesized psi (if different from observed treatment effect at t in Study A))
}
  \item{Adeltaone}{
%%     ~~Describe \code{Adeltaone} here~~
event/censoring indicators in the treatment group in Study A; optional (user must provide either (1) data from treatment arm in Study A or (2) hypothesized values for delta.ea (or R.A.given)and psi or (3) data from treatment arm in Study A and hypothesized psi (if different from observed treatment effect at t in Study A))
}
  \item{Asone}{
%%     ~~Describe \code{Asone} here~~
surrogate marker values in the treatment group in Study A, NA for individuals not observable at the time the surrogate marker was measured; optional (user must provide either (1) data from treatment arm in Study A or (2) hypothesized values for delta.ea (or R.A.given) and psi or (3) data from treatment arm in Study A and hypothesized psi (if different from observed treatment effect at t in Study A))
}
  \item{delta.ea}{
%%     ~~Describe \code{delta.ea} here~~
hypothesized value for the early treatment effect at time t0; optional (user must provide either (1) data from treatment arm in Study A or (2) hypothesized values for delta.ea (or R.A.given) and psi or (3) data from treatment arm in Study A and hypothesized psi (if different from observed treatment effect at t in Study A)), if not given then it is assumed that this quantity equals the osberved early treatment effect at time t0 in Study A
}
  \item{psi}{
%%     ~~Describe \code{psi} here~~
hypothesized value for the treatment effect at time t; optional (user must provide either (1) data from treatment arm in Study A or (2) hypothesized values for delta.ea (or R.A.given) and psi or (3) data from treatment arm in Study A and hypothesized psi (if different from observed treatment effect at t in Study A)), if not given then it is assumed that this quantity equals the osberved treatment effect at time t in Study A
}
  \item{R.A.given}{
%%     ~~Describe \code{psi} here~~
hypothesized value for the proportion of treatment effect on the primary outcome explained by surrogate information at t0 in Study A; optional (user must provide either (1) data from treatment arm in Study A or (2) hypothesized values for delta.ea (or R.A.given) and psi or (3) data from treatment arm in Study A and hypothesized psi (if different from observed treatment effect at t in Study A))
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
  \item{adjustment}{
%%     ~~Describe \code{adjustment} here~~
TRUE or FALSE; indicates whether adjustment that is needed when survival past time t is high should be used, default is FALSE if survival past t0 is < 0.90 in both arms arm of Study A, otherwise default is true if survival past t0 is >= 0.90 in either arm of Study A
}
  \item{n}{
%%     ~~Describe \code{n} here~~
total sample size for future study (Study B); optional (user needs to provide either n or power)
}
  \item{power}{
%%     ~~Describe \code{power} here~~
desired power for testing at time t0 for future study (Study B); optional (user needs to provide either n or power)
}
  \item{pi.1}{
%%     ~~Describe \code{pi.1} here~~
proportion of total sample size in future study (Study B) that would be assigned to the treatment group, default is 0.5
}
  \item{pi.0}{
%%     ~~Describe \code{pi.0} here~~
proportion of total sample size in future study (Study B) that would be assigned to the treatment group, default is 0.5
}
  \item{cens.rate}{
%%     ~~Describe \code{cens.rate} here~~
censoring in the future study (Study B) is assumed to follow an exponential distribution with censoring rate equal to this specificed value
}
  \item{transform}{
%%     ~~Describe \code{transform} here~~
TRUE or FALSE; indicates whether a transformation should be used, default is FALSE.
}
}
\details{
Assume information is available on a prior study, Study A, examining the effectiveness of a treatment up to some time of interest, t. The aim is to plan a future study, Study B, that would be conducted only up to time \eqn{t_0<t} and a test for a treatment effect would occur at \eqn{t_0}. In both studies, we assume a surrogate marker is/will be measured at time \eqn{t_0} for individuals still observable at \eqn{t_0}. Let \eqn{G} be the binary treatment indicator with \eqn{G=1} for treatment and \eqn{G=0} for control and we assume throughout that subjects are randomly assigned to a treatment group at baseline. Let \eqn{T_K^{(1)}} and \eqn{T_K^{(0)}} denote the time of the primary outcome of interest, death for example, under the treatment and under the control, respectively, in Study K. Let \eqn{S_K^{(1)}} and \eqn{S_K^{(0)}} denote the surrogate marker measured at time \eqn{t_0} under the treatment and the control, respectively, in Study K.

The null and alternative hypotheses of interest are:
\deqn{H_0: \Delta_B(t) \equiv P(T_B^{(1)}>t) - P(T_B^{(0)}>t) = 0}
\deqn{H_1: \Delta_B(t) = \psi >0}
Here, we plan to test \eqn{H_0} in Study B using the test statistic \deqn{Z_{EB}(t,t_0) = \sqrt{n_B}\frac{\hat{\Delta}_{EB}(t,t_0)}{\hat{\sigma}_{EB}(t,t_0)}} (see early.delta.test documentation). The estimated power at a type I error rate of 0.05 is thus
\deqn{ 1 - \Phi \left\{1.96 - \frac{\sqrt{n_B}\hat{R}_{SA}(t, t_0)\psi  }{ \hat{\sigma}_{EB0}(t,t_0\mid \hat{r}_A^{(0)}, W_{B}^{C})} \right \}}

where \eqn{\hat{R}_{SA}(t,t_0) =\hat{\Delta}_{EA}(t,t_0)/\hat{\Delta}_A(t)}, and \deqn{\hat{\Delta}_A(t)=n_{A1}^{-1}\sum_{i=1}^{n_{A1}}\frac{I(X_{Ai}^{(1)}>t)}{\hat{W}_{A1}^C(t)}-n_{A0}^{-1}\sum_{i=1}^{n_{A0}}\frac{I(X_{Ai}^{(0)}>t)}{\hat{W}_{A0}^C(t)},}
and \eqn{\hat{\Delta}_{EA}(t,t_0)} is parallel to \eqn{\hat{\Delta}_{EB}(t,t_0)} except replacing
\eqn{n_{A0}^{-1} \sum_{i=1}^{n_{A0}} \hat{r}_A^{(0)}(t|S_{Ai}^{(0)}, t_0) \frac{I(X_{Ai}^{(0)} > t_0)}{\hat{W}_{A0}^C(t_0)}}
by \eqn{n_{A0}^{-1}\sum_{i=1}^{n_{A0}}\hat{W}_{A0}^C(t)^{-1}I(X_{Ai}^{(0)}>t),} and \eqn{\hat{W}^C_{Ag}(\cdot)} is the Kaplan-Meier estimator of the survival function for \eqn{C_{A}^{(g)}} for \eqn{g=0,1}.  In addition, \eqn{\hat{\sigma}_{EB0}(t, t_0| \hat{r}_A^{(0)}, W_{B}^{C})^2 =} \deqn{\frac{1}{\pi_{B0}\pi_{B1}}\left[ \frac{\hat\mu_{AB2}^{(0)}(t, t_0, \mid \hat r_A^{(0)})}{W_{B}^{C}(t_0)}-\hat\mu_{AB1}^{(0)}(t, t_0, \mid \hat r_A^{(0)})^2\left\{1+\int_0^{t_0}\frac{\lambda_{B}^{C}(u)du}{\hat{W}_{A0}^{T}(u)W_{B}^{C}(u)}\right\}\right]} assuming that the survival function of the censoring distribution is \eqn{W_{B}^{C}(t)} in both arms, where \eqn{\pi_{Bg}=n_{Bg}/n_B} and \eqn{\hat{W}_{A0}^{T}(\cdot)} is the Kaplan-Meier estimator of  the survival function of \eqn{T_A^{(0)}} based on the observations from Study A, and \deqn{\hat\mu_{ABm}^{(0)}(t, t_0, \mid \hat r_A^{(0)})=n_{A0}^{-1}\sum_{i=1}^{n_{A0}}\frac{\hat{r}_A^{(0)}(t|S_{Ai}^{(0)}, t_0)^mI(X_{Ai}^{(0)}>t_0)}{\hat{W}_{A0}^{C}(t_0)}} where \eqn{\hat{r}_A^{(0)}(t|s, t_0)} is provided in the early.delta.test documentation.


This can be re-arranged to calculate the sample size needed in Study B to achieve a power of \eqn{100(1-\beta)\%}:
\deqn{n_B=\left \{ \hat{\sigma}_{EB0}(t,t_0\mid \hat{r}_A^{(0)},W_{B}^{C}) \left (\frac{1.96 - \Phi^{-1}(\beta)}{\hat{R}_{SA}(t,t_0)\psi } \right ) \right \}^2.}

When the outcome rate is low (i.e., survival rate at t is high), an adjustment to the variance calculation is needed. This is automatically implemented if the survival rate at t in either arm is 0.90 or higher. 

}
\value{
\item{n}{Total sample size needed for Study B at the given power (if power is provided by user).}
\item{power}{Estimated power for Study B at the given sample size (if sample size is provided by user).}
}
\references{
Parast L, Cai T, Tian L (2019). Using a Surrogate Marker for Early Testing of a Treatment Effect. Biometrics, In press.
}
\author{
Layla Parast
}
\examples{
data(dataA)
design.study(Axzero = dataA$x0, Adeltazero = dataA$delta0, Aszero = dataA$s0, 
Axone = dataA$x1, Adeltaone = dataA$delta1, Asone = dataA$s1, t=1, landmark=0.5, 
power = 0.80, cens.rate=0.5)

design.study(Axzero = dataA$x0, Adeltazero = dataA$delta0, Aszero = dataA$s0, 
Axone = dataA$x1, Adeltaone = dataA$delta1, Asone = dataA$s1, t=1, landmark=0.5, 
n=2500, cens.rate=0.5)

design.study(Axzero = dataA$x0, Adeltazero = dataA$delta0, Aszero = dataA$s0, 
Axone = dataA$x1, Adeltaone = dataA$delta1, Asone = dataA$s1, t=1, landmark=0.5, 
power = 0.80, cens.rate=0.5, psi = 0.05)
}
\keyword{survival}
\keyword{nonparametric}
\keyword{htest}