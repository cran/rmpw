% Please edit documentation in R/rmpw.R
\name{sensitivity}
\alias{sensitivity}
\title{Sensitivity Analysis for Causal Mediation Analysis Using Weighting Approach}
\usage{
sensitivity(est.ie, est.de, est.se.ie, est.se.de, outcome, mediator, treatment, 
    X, X.omit.pre = NULL, X.omit.post = NULL, X.unmeasure.pre = NULL, m.scale, 
    t.rand = TRUE, t.confound = FALSE, data)
}
\arguments{
\item{est.ie}{The effect size of the original natural indirect effect (NIE) estimate obtained from an RMPW analysis.}

\item{est.de}{The effect size of the original natural direct effect (NDE) estimate obtained from an RMPW analysis.}

\item{est.se.ie}{The estimated SE of the effect size of the NIE estimate.}

\item{est.se.de}{The estimated SE of the effect size of the NDE estimate.}

\item{outcome}{The name of the outcome variable (string).}

\item{mediator}{The name of the mediator variable (string).}

\item{treatment}{The name of the treatment variable (string).}

\item{X}{A vector of names of the observed pretreatment covariates already adjusted for in the original analysis (string).}

\item{X.omit.pre}{An optional vector of names of the observed pretreatment covariates that are omitted from the original analysis and may confound the mediator-outcome relationship (string). The default is NULL.}

\item{X.omit.post}{An optional vector of names of the observed pretreatment covariates that are omitted from the original analysis and may confound the mediator-outcome relationship and preceding the focal mediator (string). The default is NULL. X, X.omit.pre, and X.omit.post are mutually exclusive.}

\item{X.unmeasure.pre}{An optional vector of names of the observed pretreatment confounders that are already adjusted for in the original analysis and are comparable to some potential unmeasured confounders (String). X.unmeasure.pre is a subset of X.}

\item{m.scale}{Scale of the mediator ("discrete" or "continuous").}

\item{t.rand}{A logical value. If TRUE, treatment is randomized. If FALSE, treatment is not randomized. The default is TRUE.}

\item{t.confound}{A logical value. If TRUE, X.omit.pre may also confound the treatment-mediator or treatment-outcome relationships in addition to confounding the mediator-outcome relationship. The default is FALSE. If X.omit.pre = NULL or z.rand = TRUE, then z.confound = FALSE.}

\item{data}{A data frame containing the variables in the model.}
}
\value{
A sensitivity analysis table that contains rho, sigma, the effect size of the actual bias, the effect size of the modified estimate, and the effect size of the modified confidence interval, associated with each omitted covariate or each combination of omitted covariates, for both NIE and NDE.
}
\description{
This function generates a sensitivity analysis table. When only the mediator-outcome relationship is possibly confounded, the function computes the effect size of actual bias associated with each omitted pretreatment or posttreatment covariate or their combinations; it also computes the effect size of potential bias associated with an unmeasured confounder comparable to an observed pretreatment confounder that was already adjusted for. When the treatment assignment is also subjected to hidden selection bias, the function additionally assesses potential confounding of each omitted or unmeasured pretreatment covariate that may confound both the mediator-outcome relationship and the treatment assignment.
}
\examples{
data(Riverside)
omit.bias = sensitivity(est.ie = -0.111, est.de = 0.158, est.se.ie = 0.059,
    est.se.de = 0.108, outcome = "trunc_dep12sm2", mediator = "emp", treatment = "treat",
    X = c("emp_prior", "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30",
        "hispanic", "pqtrunc49", "nevmar"), X.omit.pre = c("AFDC3660", "pqtrunc25",
        "nohsdip"), X.omit.post = "AFDC0_Y1", m.scale = "discrete", t.rand = TRUE,
    t.confound = FALSE, data = Riverside)
}
\references{
Hong, G., Qin, X., & Yang, F. (in press). Weighting-based sensitivity analysis in causal mediation studies. Journal of Educational and Behavioral Statistics. \doi{10.3102/1076998617749561}
}
\author{
Xu Qin, Guanglei Hong, and Fan Yang
}
