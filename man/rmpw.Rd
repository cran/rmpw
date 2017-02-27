% Please edit documentation in R/rmpw.R
\name{rmpw}
\alias{rmpw}
\title{Causal Mediation Analysis Using Weighting Approach}
\usage{
rmpw(data, treatment, mediator, outcome, propensity_x, outcome_x, decomposition)
}
\arguments{
\item{data}{The data set for analysis.}

\item{treatment}{The name of the treatment variable (string).}

\item{mediator}{The name of the mediator variable (string).}

\item{outcome}{The name of the outcome variable (string).}

\item{propensity_x}{A vector of variable names (string) of pretreatment confounders, which will be included in the propensity score model.}

\item{outcome_x}{A vector of variable names (string) of pretreatment confounders, which will be included in the outcome model.}

\item{decomposition}{Type of decomposition. When decomposition = 1, the total treatment effect will be decomposed into pure direct effect (DE.0), total and pure indirect effect (IE.1 and IE.0), and natural treatment-by-mediator interaction effect (IE.1 - IE.0). When decomposition = 2, the total treatment effect will be decomposed into pure indirect effect (IE.0), total and pure direct effect (DE.1 and DE.0), and natural treatment-by-mediator interaction effect (DE.1 - DE.0).}
}
\value{
A list contains the estimates of the causal effects and the coefficients of the pretreatment covariates.
}
\description{
Causal Mediation Analysis Using Weighting Approach
}
\examples{
data(Riverside)
rmpw(data = Riverside, treatment = "treat", mediator = "emp", 
    outcome = "trunc_dep12sm2", propensity_x = c("emp_prior", 
        "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", 
        "hispanic", "pqtrunc49", "nevmar"), outcome_x = c("emp_prior", 
        "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", 
        "hispanic", "pqtrunc49", "nevmar"), decomposition = 0)
rmpw(data = Riverside, treatment = "treat", mediator = "emp", 
    outcome = "trunc_dep12sm2", propensity_x = c("emp_prior", 
        "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", 
        "hispanic", "pqtrunc49", "nevmar"), outcome_x = c("emp_prior", 
        "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", 
        "hispanic", "pqtrunc49", "nevmar"), decomposition = 1)
rmpw(data = Riverside, treatment = "treat", mediator = "emp", 
    outcome = "trunc_dep12sm2", propensity_x = c("emp_prior", 
        "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", 
        "hispanic", "pqtrunc49", "nevmar"), outcome_x = c("emp_prior", 
        "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", 
        "hispanic", "pqtrunc49", "nevmar"), decomposition = 2)
}
\references{
Hong, G., Deutsch, J., & Hill, H. D. (2015). Ratio-of-mediator-probability weighting for causal mediation analysis in the presence of treatment-by-mediator interaction. Journal of Educational and Behavioral Statistics, 40 (3), 307-340. \doi{10.3102/1076998615583902}
}
\author{
Xu Qin and Guanglei Hong
}
