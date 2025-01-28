##' Simulate for single time-step
##'
##' @param out data frame for output
##' @param proc_inputs output of `process_inputs()`
## @param control list of control parameters
##'
##' @details `sim_inversion` and `sim_rejection` correspond to
##' performing the sampling by inversion or using rejection sampling.
##'
##' @export
sim_inversion <- function (out, proc_inputs) {
  ## unpack proc_inputs
  formulas <- proc_inputs$formulas
  pars <- proc_inputs$pars
  family <- proc_inputs$family
  link <- proc_inputs$link
  estimand <- proc_inputs$estimand
  dZ <- proc_inputs$dim[1]; dX <- proc_inputs$dim[2]; dY <- proc_inputs$dim[3]
  LHS_Z <- proc_inputs$LHSs$LHS_Z; LHS_X <- proc_inputs$LHSs$LHS_X; LHS_Y <- proc_inputs$LHSs$LHS_Y; kwd <- proc_inputs$kwd
  famZ <- proc_inputs$family[[1]]; famX <- proc_inputs$family[[2]]; famY <- proc_inputs$family[[3]]; famCop <- proc_inputs$family[[4]]
  # dC <-
  
  vars <- proc_inputs$vars
  
  ## get quantiles, if any available
  if (is.null(proc_inputs$quantiles)) {
    quantiles <- out
  }
  else {
    quantiles <- cbind(proc_inputs$quantiles, out[vars])
  }
  
  
  ## code to get causal order
  order <- proc_inputs$order
  
  ## sample size
  n <- nrow(out)
  
  for (i in seq_along(order)) {
    vnm <- vars[order[i]]
    
    ## code to get Y quantiles conditional on different Zs
    if (order[i] > dZ+dX) {
      ## simulate Y variable
      # qY <- runif(n)
      wh <- order[i] - dZ - dX
      # print(wh)
      
      ## code to use sim_variable
      forms <- list(formulas[[3]][[wh]], formulas[[4]][[wh]])
      fams <- list(family[[3]][[wh]], family[[4]][[wh]])
      prs <- list(pars[[vnm]], pars[[kwd]][[vnm]])
      lnk <- list(link[[3]][wh], list()) # link[[4]][[wh]])
      
      if (any(is.na(lhs(forms[[2]])))) {
        forms[[2]] <- `lhs<-`(forms[[2]], c(LHS_Z[rank(order[seq_len(dZ)])],
                                            LHS_Y[rank(order[dZ+dX+seq_len(i-1 - dZ-dX)])]))
      }
      
      ## code to rescale qU for different estimands
      #U <- as.numeric(out[[LHS_Z[1]]])
      if (estimand == "att"){
        #qU <- ecdf(U[out[LHS_X]==1])(U)
        #quantiles[LHS_Z[1]] <- qU
        quantiles[LHS_Z] <- scaled_qZs
      }
      else if (estimand == "atc"){
        qU <- ecdf(U[out[LHS_X]==0])(U)
        quantiles[LHS_Z[1]] <- qU
      }
      else if (estimand == "ato"){
        # to do
      }
      
      out <- sim_variable(n=nrow(out), formulas=forms, family=fams, pars=prs,
                          link=lnk, dat=out, quantiles=quantiles)
      quantiles <- attr(out, "quantiles")
      attr(out, "quantiles") <- NULL
      
      # out[[vars[order[i]]]] <- sim_Y(n, formulas=formulas[[4]][[wh]],
      #                                family=family[[4]][[wh]],
      #                                pars=pars[[kwd]][[wh]],
      #                                formY = formulas[[3]][[wh]],
      #                                famY=family[[3]][wh],
      #                                parsY=pars[[LHS_Y[wh]]],
      #                                linkY=link[[3]][wh], qZ=quantiles, vars=vars,
      #                                dat=out)
      
      # for (j in seq_len(dZ)) {
      #   curr_qZ <- qZs[[vars[j]]]
      #   X <- model.matrix(formulas[[4]][[wh]][[j]], data=out)
      #   curr_fam <- family[[4]][wh,j]
      #   curr_par <- pars[[kwd]]$beta[[wh]][[j]]
      #   # eta <- X %*% curr_par
      #   qY <- rescale_cop(cbind(curr_qZ,qY), X=X, pars=curr_par, family=curr_fam) #, link=link[[4]][i,j])
      # }
      # ##
      # X <- model.matrix(formulas[[3]][[wh]], data=out)
      # qY <- rescale_var(qY, X=X, family=famY[[wh]], pars=pars[[LHS_Y[wh]]],
      #                  link=link[[3]][wh])
      #
      # out[[vars[order[i]]]] <- qY
    }
    else {
      ## code to simulate Z and X variables in causal order
      curr_link <- unlist(link)[order[i]]
      
      if (vnm %in% LHS_Z) {
        print("simulating Zs")
        curr_form <- formulas[[1]][[order[i]]]
        curr_fam <- famZ[[order[i]]]
        
        trm <- terms(curr_form)
        # curr_form2 <- delete.response(terms(curr_form))
        MM <- model.matrix(delete.response(trm), data=out)
        if (nrow(MM) != nrow(out)) {
          if (length(attr(trm, "factors")) == 0) {
            if (attr(trm, "intercept") == 1) MM <- matrix(1, nrow=nrow(out), ncol=1)
            else MM <- matrix(0, nrow=nrow(out), ncol=0)
          }
          else warning(paste0("Missing entries for ", vnm))
        }
        eta <- MM %*% pars[[vnm]]$beta
        oth_pars <- pars[[vnm]]
        curr_phi <- pars[[vnm]]$phi
        tmp <- glm_sim(family=curr_fam, eta=eta, phi=curr_phi, other_pars=pars[[vnm]], link=curr_link)
        if (vnm %in% LHS_Z) quantiles[[vnm]] <- attr(tmp, "quantile")
        attr(tmp, "quantile") <- NULL
        out[[vnm]] <- tmp
      }
      else{
        print("Simulating X")
        # X_curr_form <- formulas[[2]][[order[i]-dZ]]
        # X_fam <- famX[[order[i]-dZ]][[1]] 
        X_fam <- "logistic" # now simply set as logistic distribution
        # copX_fam <- famX[[order[i]-dZ]][[2]]
        # copX_fam <- famCop[[1]]
        copX_fam <- 1
        X_res <- sim_cX(n, family=list(X_fam,copX_fam), 
                        pars= list(list(beta=0,phi=1), list(0.5,0.2)),
                        dat=out[LHS_Z], quantiles=quantiles[LHS_Z])
        #e.g.pars=list(list(beta=0, phi=1), list(0.5,0.2,0.4))
        out[[vnm]] <- X_res$X_samples
        scaled_qZs <- X_res$scaled_qZs
        
      }
    }
  }
  
  attr(out, "qZ") <- quantiles
  
  return(out)
}

