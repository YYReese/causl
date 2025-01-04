##' Checking the ATE given a causl model
##'
##' @param n number of samples required
##' @param causl_model object of class `causl_model`
##' @param formulas list of lists of formulas
##' @param family families for variables and copula
##' @param pars list of lists of parameters
##' @param link list of link functions
##' @param estimand either `"ate"` (the default), `"att"`, `"att"`, or `"ato"`
##' @param method either `"inversion"` (the default), `"inversion_mv"`, or `"rejection"`
##' @param control list of options for the algorithm
##' 
##' @return a list containing correct and biased ATE with estimated se 
##' using sandwich method
##' 
##' @export
causl_checking <- function(n, formulas = forms, family = fams, 
                           pars=pars, link=NULL, estimand="ate",
                           method="inversion",control=list()){
  # Sample from the model (default ATE)
  #samples <- rfrugal(n, model, control=control)
  samples <- rfrugalParam(n, formulas,
                          family, pars, estimand=estimand, method="inversion")
  dato <- samples[!is.infinite(samples$Y),]
  
  # Regression checking
  # simplest case here, X is 1-dim treatment set 
  X_form <- formulas[[2]]
  X_fam <- family_vals$family[family_vals$val == family[[2]]]
  modX <- glm(X_form, family=X_fam, data=dato) 
  ps <- predict(modX, type = "response")
  
  # Trim extreme propensity scores to avoid instability
  ps <- pmax(pmin(ps, 0.99), 0.01) 
  if (estimand == "ate"){
    # weights for ATE (1/ps for the treated, 1/(1-ps) for the controlled)
    wt <- dato$X/ps + (1-dato$X)/(1-ps) # weights for ATE
  }
  else if (estimand == "att"){
    # weights for ATT (1 for the treated, ps/(1-ps) for the controlled)
    wt <- dato$X + (1-dato$X)*ps/(1-ps) 
  }
  else if (estimand == "atc"){
    # weights for ATC ((1-ps)/ps for the treated, 1 for the controlled)
    wt <- dato$X*(1-ps)/ps + (1-dato$X)
  }
  else { 
    # weights for ATO ((1-ps)*ps for both)
    wt <- (1-ps)*ps
    wt <- wt/sum(wt)
  }
  # Correct model for ATE
  # simplest case here, Y is 1-dim outcome set 
  Y_form <- formulas[[3]]
  Y_fam <- family_vals$family[family_vals$val == family[[3]]]
  dato$wt <- wt
  modY <- glm(Y_form, family=Y_fam, weights=wt, data=dato) 
  
  est_correct <- modY$coefficients[2] # estimand estimate
  seY <- sqrt(diag(sandwich::sandwich(modY)))[2] # standard error estimate
  
  # Unweighted (biased) model for ATE
  modYw <- glm(Y_form, family=Y_fam, data=dato) 
  est_biased <- modYw$coefficients[2] # estimand biased estimate
  seY_biased <- sqrt(diag(sandwich::sandwich(modYw)))[2] 
  
  return(list(correct = list(estimand=est_correct, se=seY),
              biased = list(estimand = est_biased, se=seY_biased)))
}

# Example
#n <- 1e5
#forms <- list(list(U ~ 1,Z~1), 
#              X ~ U+Z,
#              Y ~ X, 
#              ~ 1)

#pars <- list(Z = list(beta = 0, phi=1), #Z~Gamma(shape=1/phi,scale=phi*inv_link(mu)), link(mu)=eta, eta=X*beta
#             U = list(beta = 0, phi=1),
#             X = list(beta = c(0,0.2,1)),
#             Y = list(beta=c(0, 0.5), phi=1),
 #            cop = list(beta = c(0.5)))

#fams <- list(list(1,3),5,1,1) # 1 for Gaussian, 3 for Gamma, 5 for binomial


#rfrugalParam(n, formulas = forms,
 #             family = fams, pars=pars, estimand="ate",
#              method="inversion")


# Define the model
#model <- causl_model(formulas=forms, family=fams, pars=pars, dat=NULL, method="inversion",
#                     kwd="cop")

#debugonce(causl_checking)
#out <- causl_checking(n, model)
