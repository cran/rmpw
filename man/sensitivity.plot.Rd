% Please edit documentation in R/rmpw.R
\name{sensitivity.plot}
\alias{sensitivity.plot}
\title{Sensitivity Analysis Plot for Causal Mediation Analysis Using Weighting Approach}
\usage{
sensitivity.plot(est.ie, est.de, est.se.ie, est.se.de, X.omit.bias, effect, type)
}
\arguments{
\item{est.ie}{The effect size of the original natural indirect effect (NIE) estimate obtained from an RMPW analysis.}

\item{est.de}{The effect size of the original natural direct effect (NDE) estimate obtained from an RMPW analysis.}

\item{est.se.ie}{The estimated SE of the effect size of the NIE estimate.}

\item{est.se.de}{The estimated SE of the effect size of the NDE estimate.}

\item{X.omit.bias}{The bias table returned by sensitivity()}

\item{effect}{If effect = "NIE", returns sensitivity plot for the natural indirect effect; If effect = "NDE", returns sensitivity plot for the natural direct effect}

\item{type}{If type = "unmeasured", returns sensitivity plot for unmeasured covariates that are comparable to observed pretreatment covariates; if type = "omitted", returns sensitivity plot for omitted covariates. If in sensitivity(), X.unmeasure.pre = NULL, then here type can only be "omitted". If in sensitivity(), X.omit.pre = NULL and X.omit.post = NULL, then here type can only be "unmeasured".}
}
\value{
A graphical display of the sensitivity analysis results returned by sensitivity().
}
\description{
This function generates sensitivity analysis plots, one for natural direct effect (NDE) and one for natural indirect effect (NIE). The solid curves indicate amounts of bias that may lead to a qualitative change in the conclusion of statistical inference. The number at the end of each curve denotes the bias value represented by that curve. The bold dashed curve on each side corresponds to a threshold value of bias that is just great enough to start changing the analytic conclusion. Each star corresponds to an observed covariate or a set of covariates that, if omitted, would contribute a bias great enough to change the analytic conclusion; each triangle corresponds to a covariate (or a set of covariates), the omission of which would be inconsequential. These covariates are listed in the sensitivity analysis table in the output of sensitivity(). This function is not applicable when there are omitted or unmeasured confounders for the treatment assignment. This is because, in such cases, there are as many as four sensitivity parameters, difficult to represent in a single plot.
}
\examples{
data(Riverside)
omit.bias = sensitivity(est.ie = -0.111, est.de = 0.158, est.se.ie = 0.059, 
    est.se.de = 0.108, outcome = "trunc_dep12sm2", mediator = "emp", 
    treatment = "treat", X = c("emp_prior", "pqtrunc50", "pqtrunc51", 
        "pqtrunc52", "pqtrunc53", "pqtrunc30", "hispanic", "pqtrunc49", 
        "nevmar"), X.omit.pre = c("AFDC3660", "pqtrunc25", "nohsdip"), 
    X.omit.post = "AFDC0_Y1", m.scale = "discrete", t.rand = TRUE, t.confound = FALSE, 
    data = Riverside)
sensitivity.plot(est.ie = -0.111, est.de = 0.158, est.se.ie = 0.059, 
    est.se.de = 0.108, X.omit.bias = omit.bias, effect = "NIE", type = "omitted")
sensitivity.plot(est.ie = -0.111, est.de = 0.158, est.se.ie = 0.059, 
    est.se.de = 0.108, X.omit.bias = omit.bias, effect = "NDE", type = "omitted")
}
\references{
Hong, G., Qin, X., & Yang, F. (in press). Weighting-based sensitivity analysis in causal mediation studies. Journal of Educational and Behavioral Statistics. \doi{10.3102/1076998617749561}
}
\author{
Xu Qin, Guanglei Hong, and Fan Yang
}
