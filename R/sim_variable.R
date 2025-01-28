##' Simulate a single variable using the inversion method
##'
##' @param n sample size
##' @param formulas list consisting of a formula for the output variables and a list of formulae for the pair-copula
##' @param family list containing family variable
##' @param pars list with two entries, first a list of parameters for response, and second a further list of parameters for pair-copula
##' @param link list of same form as `family`
##' @param dat data frame of current variables
##' @param quantiles data frame of quantiles
##'
##' @return The data frame `dat` with an additional column given by the left-hand side of `formula[[1]]`.
##'
##' @description
##' Each entry `formulas`, `family`, `pars`, `link` is a list
##' with two entries, the first referring to the variable being simulated and the
##' second to the pair-copulas being used.
##'
##' @export
sim_variable <- function (n, formulas, family, pars, link, dat, quantiles) {
  qY <- runif(n)

  ## get variables
  vnm <- lhs(formulas[[1]])
  if (length(vnm) != 1L) stop("Unable to extract variable name")

  quantiles[[vnm]] <- qY

  LHS_cop <- lhs(formulas[[2]])

  for (i in rev(seq_along(formulas[[2]]))) {
    X <- model.matrix(formulas[[2]][[i]], data=dat)
    # eta <- X %*% pars[[2]][[i]]$beta

    ## rescale quantiles for pair-copula
    qs <- cbind(quantiles[[LHS_cop[[i]]]], qY)
    qY <- rescale_cop(qs, X=X, beta=pars[[2]][[i]]$beta, family=family[[2]][[i]],
                     par2=pars[[2]][[i]]$par2) #, link=link[[2]][j])
  }

  ## now rescale to correct margin
  X <- model.matrix(delete.response(terms(formulas[[1]])), data=dat)

  Y <- rescale_var(qY, X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  # Y <- rescale_var(runif(n), X=X, family=family[[1]], pars=pars[[1]], link=link[[1]])
  dat[[vnm]] <- Y
  # quantiles[[vnm]] <- qY
  attr(dat, "quantiles") <- quantiles

  return(dat)
}

##' @param n number of samples
##' @param family list containing family variable for X and pair-copulas
##' @param pars list with two entries, first a list of parameters for X, and second a further list of parameters for pair-copulas
##' @param dat data frame of current variables Z1,Z2|Z1,...ZK|Zk-1:1
##' @param quantiles data frame of quantiles
##' 
##' @return a matrix with each column being F_{X|Zk:1} for all k<=K, 
##' and T|ZK:1 samples 
sim_cX <- function(n, family, pars, dat, quantiles){
  # sample from X marginal 
  if (family[[1]] == "logistic"){
    beta <- pars[[1]][[1]]
    phi <- pars[[1]][[2]]
    t_star <- rlogis(n, beta,phi)
    ut_star <- plogis(t_star,beta,phi)
    threshold <- 0 # T*=0 as threshold. T=1 if T*>0; T=0 if T*<=0
    u_threshold <- plogis(threshold,beta,phi) 
  }
  dimZ <- dim(quantiles)[2]
  
  # Initialise matrices for conditional probabilities
  p <- matrix(ncol=dimZ+1, nrow=n) # P(T*<t|Zk:1)
  ps <- matrix(ncol=dimZ+1,nrow=n) # P(T=1|Zk:1) = P(T*>0|Zk:1)
  scaled_qZs <- matrix(ncol=dimZ,nrow=n)
  p[,1] <- ut_star
  ps[,1] <- rep(1-u_threshold,n)
  
  params <- pars[[2]]
  for (i in 1:dimZ){
    if (family[[2]] == 1){
      #qnorm(U[,2])*sqrt(1-param^2)+param*qnorm(U[,1])
      #p[,i+1] <- pnorm((qnorm(p[,i])  - params[[i]] *qnorm( quantiles[,i])) / sqrt(1 - params[[i]]^2))
      p[,i+1] <- pnorm(qnorm(p[,i])*sqrt(1-params[[i]]^2)+params[[i]]*qnorm(quantiles[,i]))
      #ps[,i+1] <- 1-pnorm((qnorm((1-ps[,i]))  - params[[i]] *qnorm( quantiles[,i])) / sqrt(1 - params[[i]]^2))
      ps[,i+1] <- 1-pnorm(qnorm(1-ps[,i])*sqrt(1-params[[i]]^2)+params[[i]]*qnorm(quantiles[,i]))
      
      #integrand <- function(uz) {
      #  1-pnorm((qnorm((1-ps[,i]))  - params[[i]]*qnorm(uz)) / sqrt(1 - params[[i]]^2))
      #}
      input_data <- data.frame(ps[,i],quantiles[,i])
      scaled_qZs[,i] <- apply(input_data, 1, function(row) {
        integrand <- function(uz) {
          1 - pnorm((qnorm((1 - row[1])) - params[[i]] * qnorm(uz)) / sqrt(1 - params[[i]]^2))
        }
        integrate(integrand, lower = 0, upper = row[2])$value / row[1]
      })
    }
    ## (to be completed)
    else if (family == 2) {
      # Y <- sqrt(phi)*qt(U, df=pars$par2) + eta
      p[,i+1] <- cVCopula(data.frame(quantiles[,i],p[,i]), copula = tCopula, param = params[[i]][[1]], par2=params[[i]][[2]], inverse=FALSE)
      ps[,i+1] <- 1-cVCopula(data.frame(quantiles[,i],1-ps[,i]), copula = tCopula, param = params[[i]][[1]], par2=params[[i]][[2]], inverse=FALSE)
    }
    else if (family == 3) {
      p[,i+1] <- cVCopula(data.frame(quantiles[,i],p[,i]), copula = claytonCopula, param = param[[i]], inverse=FALSE)
      ps[,i+1] <- 1-cVCopula(data.frame(quantiles[,i],1-ps[,i]), copula = tCopula, param = params[[i]], par2=params[[i]][[2]], inverse=FALSE)
    }
    else if (family == 4) {
      p[,i+1] <- cVCopula(data.frame(quantiles[,i],p[,i]), copula = gumbelCopula, param = params[[i]], inverse=FALSE)
      ps[,i+1] <- 1-cVCopula(data.frame(quantiles[,i],1-ps[,i]), copula = tCopula, param = params[[i]], par2=params[[i]][[2]], inverse=FALSE)
      
    }
    else if (family == 5) {
      p[,i+1] <- cVCopula(data.frame(quantiles[,i],p[,i]), copula = frankCopula, param = params[[i]], inverse=FALSE)
      ps[,i+1] <- 1-cVCopula(data.frame(quantiles[,i],1-ps[,i]), copula = tCopula, param = params[[i]], par2=params[[i]][[2]], inverse=FALSE)
    }
    else if (family == 6) {
      p[,i+1] <- cVCopula(data.frame(quantiles[,i],p[,i]), copula = joeCopula, param = param[[i]], inverse=FALSE)
      ps[,i+1] <- 1-cVCopula(data.frame(quantiles[,i],1-ps[,i]), copula = tCopula, param = params[[i]], par2=params[[i]][[2]], inverse=FALSE)
    }
    else stop("family must be between 0 and 5")
  }
  # simulate discrete X
  X <- ifelse(p[,dimZ+1] > u_threshold, 1, 0)
  
  return(list(ps = ps, X_samples = X, scaled_qZs=scaled_qZs ))
}

