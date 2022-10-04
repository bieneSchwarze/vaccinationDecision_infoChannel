fitModel <- function(imp, eq, D_cc, we=TRUE){
  
  # prestuff
  # -----------------------------------------------------------------
  m_cc <- glm(formula=eq, family=binomial(link="logit"), data=D_cc)
  #summary(m_cc)
  #m_cc_cl <- coeftest(m_cc, vcov. = vcovCL(m_cc, cluster = D_cc$idh, type = "HC0"))
  #dfobj <- svydesign(id = ~id, weights = ~phrf_full, data = D_cc)
  #m_ww_cc <- svyglm(formula=eq,design=dfobj,family = quasibinomial)
  #summary(m_ww_cc)
  #m_cc_ww_cl <- coeftest(m_ww_cc, vcov. = vcovCL(m_ww_cc, cluster = D_cc$idh, type = "HC0"))
  #stargazer(m_cc_cl, m_cc_ww_cl, type="text")
  # -----------------------------------------------------------------
  
  # Do this looping over all imputed data sets, then pool
  modList <- vector(length=imp$m, mode="list") 
  nbetas <- nrow(summary(m_cc)$coef)
  matrBeta <- matrix(NA, ncol=imp$m, nrow=nbetas)
  matrSE <- matrix(NA, ncol=imp$m, nrow=nbetas)
  for(i in 1:imp$m){
    #cat("Imp.It: ",i,"\n")
    #i <- 1
    D_i <- complete(imp, action=i)
    if(we){
      dfobj <- svydesign(id = ~id, weights = ~phrf_full, data = D_i)
    } else {
      dfobj <- svydesign(id = ~id, weights = ~1, data = D_i)
    }
    mod_i <- svyglm(formula=eq,design=dfobj,family = quasibinomial)
    mod_i_cl <- coeftest(mod_i, vcov. = vcovCL(mod_i, cluster = D_i$idh, type = "HC0"))
    modList[[i]] <- mod_i_cl
    matrBeta[,i] <- mod_i_cl[, "Estimate"]
    matrSE[,i] <- mod_i_cl[, "Std. Error"]  
  }
  betas <- apply(matrBeta, 1, mean)
  var_within <- apply(matrSE^2, 1, mean)
  var_between <- apply((matrSE-var_within)^2,1,mean)/(imp$m-1)
  var_total <- (1+(1/imp$m))*var_between + var_within   
  wald_pooled <- betas^2/var_total
  lambda <- (var_between + var_between/imp$m)/var_total #  fraction of missing information
  df_old <- (imp$m-1)/lambda^2
  df_obs <- ((nrow(D_cc)-nbetas) + 1)/((nrow(D_cc)-nbetas) + 3) * (nrow(D_cc)-nbetas)*(1-lambda)
  df_adj <- df_old*df_obs/(df_old + df_obs)
  alpha <- 0.05
  cis <- cbind(betas-qt(1-alpha/2, df=df_adj)*sqrt(var_total),betas+qt(1-alpha/2, df=df_adj)*sqrt(var_total)) # 95% CI
  pvalue <- pt(-abs(wald_pooled),df=df_adj)             
  res <- cbind(betas, pvalue, cis)
  rownames(res) <- rownames(mod_i_cl)
  colnames(res) <- c("Estimate","Pr(>|t|)", "CI_low", "CI_up")
  return(res)
}
