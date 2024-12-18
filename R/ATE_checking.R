##' Checking the ATE given a causl model
##'
##' @param n number of samples to approximate
##' @param model object of class `causl_model`
##' @param control list of options for the algorithm
##' 
##' @return a list containing correct and biased ATE with estimated se 
##' using sandwich method
causl_checking <- function(n, model, control=list()){
  # Sample from the model (default ATE)
  samples <- rfrugal(n, model, control=control)
  dato <- samples[!is.infinite(samples$Y),]
  
  # Regression checking
  # simplest case here, X is 1-dim treatment set 
  X_form <- model$formulas[[2]][[1]]
  X_fam <- family_vals$family[family_vals$val == model$family[[2]]]
  modX <- glm(X_form, family=X_fam, data=dato) 
  ps <- predict(modX, type = "response") 
  wt <- dato$X/ps + (1-dato$X)/(1-ps) # weights for ATE
  
  # Correct model for ATE
  # simplest case here, Y is 1-dim outcome set 
  Y_form <- model$formulas[[3]][[1]]
  Y_fam <- family_vals$family[family_vals$val == model$family[[3]]]
  dato$wt <- wt
  modY <- glm(Y_form, family=Y_fam, weights=wt, data=dato) 
  
  ate_correct <- modY$coefficients[2] # ate estimate
  seY <- sqrt(diag(sandwich::sandwich(modY)))[2] # ate standard error estimate
  
  # Unweighted (biased) model for ATE
  modYw <- glm(Y_form, family=Y_fam, data=dato) 
  ate_biased <- modYw$coefficients[2] # ate biased estimate
  seY_biased <- sqrt(diag(sandwich::sandwich(modYw)))[2] 
  
  return(list(correct = list(ate=ate_correct, se=seY),
              biased = list(ate = ate_biased, se=seY_biased)))
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
