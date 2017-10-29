\name{early.delta.test}
\alias{early.delta.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate and test the early treatment effect
}
\description{
Estimates the early treatment effect estimate and provides two versions of the standard error; tests the null hypothesis that this treatment effect is equal to 0 }
\usage{
early.delta.test(Axzero, Adeltazero, Aszero, Bxzero, Bdeltazero, Bszero, Bxone, 
Bdeltaone, Bsone, t, landmark, perturb = T, extrapolate = T)
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
  \item{perturb}{
%%     ~~Describe \code{perturb} here~~
TRUE or FALSE; indicates whether the standard error estimate obtained using perturbation resampling should be calculated
}
  \item{extrapolate}{
%%     ~~Describe \code{extrapolate} here~~
TRUE or FALSE; indicates whether local constant extrapolation should be used, default is TRUE
}
}
\details{
Assume there are two randomized studies of a treatment effect, a prior study (Study A) and a current study (Study B). Study A was completed up to some time t, while Study B was stopped at time \eqn{t_0<t}. In both studies, a surrogate marker was measured at time \eqn{t_0} for individuals still observable at \eqn{t_0}. Let \eqn{G} be the binary treatment indicator with \eqn{G=1} for treatment and \eqn{G=0} for control and we assume throughout that subjects are randomly assigned to a treatment group at baseline. Let \eqn{T_K^{(1)}} and \eqn{T_K^{(0)}} denote the time of the primary outcome of interest, death for example, under the treatment and under the control, respectively, in Study K. Let \eqn{S_K^{(1)}} and \eqn{S_K^{(0)}} denote the surrogate marker measured at time \eqn{t_0} under the treatment and the control, respectively, in Study K.

The treatment effect quantity of interest, \eqn{\Delta_K(t)}, is the difference in survival rates by time \eqn{t} under treatment versus control,
  \deqn{ \Delta_K(t)=E\{ I(T_K^{(1)}>t)\} - E\{I(T_K^{(0)}>t)\} = P(T_K^{(1)}>t) - P(T_K^{(0)}>t)}
where \eqn{t>t_0}. Here, we estimate an early treatment effect quantity using surrogate marker information defined as, 
\deqn{\Delta_{EB}(t,t_0) = P( T_B^{(1)} > t_0) \int r(t|s,t_0) dF_B^{(1)} (s|t_0) - P( T_B^{(0)} > t_0) \int r(t|s,t_0) dF_B^{(0)} (s|t_0)}
where \eqn{r(t|s,t_0) = P(T_{A}^{(0)} > t | T_{A}^{(0)} > t_0, S_{A}^{(0)}=s)} and \eqn{F_B^{(g)}(s|t_0) = P(S_B^{(g)} \le s \mid T_B^{(g)} > t_0)}.

To test the null hypothesis that \eqn{\Delta_B(t) = 0}, we test the null hypothesis \eqn{\Delta_{EB}(t,t_0) = 0} using the test statistic \deqn{Z_{EB}(t,t_0) = \sqrt{n_B}\frac{\hat{\Delta}_{EB}(t,t_0)}{\hat{\sigma}_{EB}(t,t_0)}}
where \eqn{\hat{\Delta}_{EB}(t,t_0)} is a consistent estimate of \eqn{\Delta_{EB}(t,t_0)} and \eqn{\hat{\sigma}_{EB}(t,t_0)} is the estimated standard error of \eqn{\sqrt{n_B}\{\hat{\Delta}_{EB}(t,t_0)-\Delta_{EB}(t, t_0)\}.}  We reject the null hypothesis when \eqn{|Z_{EB}(t,t_0) | > \Phi^{-1}(1-\alpha/2)} where \eqn{\alpha} is the Type 1 error rate. 

To obtain \eqn{\hat{\Delta}_{EB}(t,t_0)}, we use 
\deqn{  \hat{\Delta}_{EB}(t,t_0) = n_{B1}^{-1} \sum_{i=1}^{n_{B1}} \hat{r}_A^{(0)}(t|S_{Bi}^{(1)}, t_0) \frac{I(X_{Bi}^{(1)} > t_0)}{\hat{W}_{B1}^C(t_0)} -  n_{B0}^{-1} \sum_{i=1}^{n_{B0}} \hat{r}_A^{(0)}(t|S_{Bi}^{(0)}, t_0) \frac{I(X_{Bi}^{(0)} > t_0)}{\hat{W}_{B0}^C(t_0)}}
where \eqn{\hat{W}^C_{k g}(u)} is the Kaplan-Meier estimator of \eqn{W_{k g}^{C}(u)=P(C_{k}^{(g)} > u)}  and 
 \eqn{\hat{r}_A^{(0)}(t|s,t_0) = \exp\{-\hat{\Lambda}_A^{(0)}(t\mid s,t_0) \}}, where
\deqn{\hat{\Lambda}_A^{(0)}(t \mid t_0,s) = \int_{t_0}^t \frac{\sum_{i=1}^{n_{A0}} I(X_{Ai}^{(0)}>t_0) K_h\{\gamma(S_{Ai}^{(0)}) - \gamma(s)\}dN_{Ai}^{(0)} (z)}{\sum_{i=1}^{n_{A0}}  K_h\{\gamma(S_{Ai}^{(0)}) - \gamma(s)\} Y_{Ai}^{(0)}(z)}}
is a consistent estimate of \eqn{\Lambda_A^{(0)}(t\mid t_0,s ) =  -\log [r_A^{(0)}(t\mid t_0,s)],} \eqn{Y_{Ai}^{(0)}(t) = I(X_{Ai}^{(0)} \geq t)},  \eqn{N_{Ai}^{(0)}(t) = I(X_{Ai}^{(0)} \leq t) \delta_{Ai}^{(0)},  K(\cdot)} is a smooth symmetric density function, \eqn{K_h(x) = K(x/h)/h} and \eqn{\gamma(\cdot)} is a given monotone transformation function.  For the bandwidth \eqn{h}, we require the standard undersmoothing assumption of \eqn{h=O(n_g^{-\gamma})} with \eqn{\gamma \in (1/4,1/2)} in order to eliminate the impact of the bias of the conditional survival function on the resulting estimator.

The quantity \eqn{\hat{\sigma}_{EB}(t,t_0)} is obtained using either a closed form expression under the null or a perturbation resampling approach. If a confidence interval is desired, perturbation resampling is required.
} 

\value{
\item{delta.eb}{The estimate early treatment effect, \eqn{\hat{\Delta}_{EB}(t,t_0)}.}
\item{se.closed}{The standard error estimate of the early treatment effect using the closed form expression under the null.}
\item{Z.closed}{The test statistic using the closed form standard error expression.}
\item{p.value.closed}{The p-value using the closed form standard error expression.}
\item{se.perturb}{The standard error estimate of the early treatment effect using perturbation resampling, if perturb = T.} 
\item{Z.perturb}{The test statistic using the perturbed standard error estimate, if perturb = T.}
\item{p.value.perturb}{The p-value using the perturbed standard error estimate, if perturb = T}
\item{delta.eb.CI}{The confidence interval for the early treatment effect, using perturbation resampling, if perturb = T}
}
\references{
%% ~put references to the literature/web site here ~
Parast L, Cai T, Tian L (2017). Using a Surrogate Marker for Early Testing of a Treatment Effect, Under Review.
}
\author{
%%  ~~who you are~~
Layla Parast
}
\examples{
data(dataA)
data(dataB)
early.delta.test(Axzero = dataA$x0, Adeltazero = dataA$delta0, Aszero = dataA$s0, 
Bxzero = dataB$x0, Bdeltazero = dataB$delta0, Bszero = dataB$s0, Bxone = dataB$x1, 
Bdeltaone = dataB$delta1, Bsone = dataB$s1, t=1, landmark=0.5, perturb = FALSE, 
extrapolate = TRUE)

early.delta.test(Axzero = dataA$x0, Adeltazero = dataA$delta0, Aszero = dataA$s0, 
Bxzero = dataB$x0, Bdeltazero = dataB$delta0, Bszero = dataB$s0, Bxone = dataB$x1, 
Bdeltaone = dataB$delta1, Bsone = dataB$s1, t=0.75, landmark=0.5, perturb = FALSE, 
extrapolate = TRUE)

\donttest{
early.delta.test(Axzero = dataA$x0, Adeltazero = dataA$delta0, Aszero = dataA$s0, 
Bxzero = dataB$x0, Bdeltazero = dataB$delta0, Bszero = dataB$s0, Bxone = dataB$x1, 
Bdeltaone = dataB$delta1, Bsone = dataB$s1, t=1, landmark=0.5, perturb = TRUE, 
extrapolate = TRUE)	
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% use one of  RShowDoc("KEYWORDS")
\keyword{nonparametric}% __ONLY ONE__ keyword per line
\keyword{htest}
