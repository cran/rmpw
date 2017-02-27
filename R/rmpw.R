#' Causal Mediation Analysis Using Weighting Approach
#'
#'
#' @param data The data set for analysis.
#' @param treatment The name of the treatment variable (string).
#' @param mediator The name of the mediator variable (string).
#' @param outcome The name of the outcome variable (string).
#' @param propensity_x A vector of variable names (string) of pretreatment confounders, which will be included in the propensity score model. 
#' @param outcome_x A vector of variable names (string) of pretreatment confounders, which will be included in the outcome model. 
#' @param decomposition Type of decomposition. When decomposition = 1, the total treatment effect will be decomposed into pure direct effect (DE.0), total and pure indirect effect (IE.1 and IE.0), and natural treatment-by-mediator interaction effect (IE.1 - IE.0). When decomposition = 2, the total treatment effect will be decomposed into pure indirect effect (IE.0), total and pure direct effect (DE.1 and DE.0), and natural treatment-by-mediator interaction effect (DE.1 - DE.0).
#' @return A list contains the estimates of the causal effects and the coefficients of the pretreatment covariates.
#' @author Xu Qin and Guanglei Hong
#' @references Hong, G., Deutsch, J., & Hill, H. D. (2015). Ratio-of-mediator-probability weighting for causal mediation analysis in the presence of treatment-by-mediator interaction. Journal of Educational and Behavioral Statistics, 40 (3), 307-340. \doi{10.3102/1076998615583902}
#' @export
#' @importFrom stats as.formula binomial coef fitted glm lm pnorm predict
#' @examples 
#' data(Riverside)
#' rmpw(data = Riverside, treatment = "treat", mediator = "emp", outcome = "trunc_dep12sm2", propensity_x = c("emp_prior", "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", "hispanic", "pqtrunc49", "nevmar"), outcome_x = c("emp_prior", "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", "hispanic", "pqtrunc49", "nevmar"), decomposition = 0)
#' rmpw(data = Riverside, treatment = "treat", mediator = "emp", outcome = "trunc_dep12sm2", propensity_x = c("emp_prior", "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", "hispanic", "pqtrunc49", "nevmar"), outcome_x = c("emp_prior", "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", "hispanic", "pqtrunc49", "nevmar"), decomposition = 1)
#' rmpw(data = Riverside, treatment = "treat", mediator = "emp", outcome = "trunc_dep12sm2", propensity_x = c("emp_prior", "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", "hispanic", "pqtrunc49", "nevmar"), outcome_x = c("emp_prior", "pqtrunc50", "pqtrunc51", "pqtrunc52", "pqtrunc53", "pqtrunc30", "hispanic", "pqtrunc49", "nevmar"), decomposition = 2)

rmpw = function(data, treatment, mediator, outcome, propensity_x, outcome_x, decomposition){
  if(decomposition == 0){
    weight0 = function(m, t, X, data){
      data1 = data[which(data[, t] ==1), ]
      data0 = data[which(data[, t] ==0), ]
      formula = as.formula(paste(m, "~", paste(X, collapse="+")))
      formula_x = as.formula(paste(t, "~", paste(X, collapse="+")))
      l_x = glm(formula_x, data = data, family = binomial)
      p1_x = fitted(l_x)
      p0_x = 1 - p1_x
      iptw1 = (sum(data[, t] == 1)/nrow(data))/p1_x
      iptw0 = (sum(data[, t] == 0)/nrow(data))/p0_x
      data$iptw[data[, t] == 1] = iptw1[data[, t] == 1] # For the use in the bias calculation when T is not randomized and X.omit is pretreatment
      data$iptw[data[, t] == 0] = iptw0[data[, t] == 0]
      l1 = glm(formula, data = data1, family = binomial)
      l0 = glm(formula, data = data0, family = binomial)
      p1 = predict(l1, data, type = "response")
      p0 = predict(l0, data, type = "response")
      data$rmpw[data[, t] == 1 & data[, m] == 1] = (p0/p1)[data[, t] == 1 & data[, m] == 1]
      data$rmpw[data[, t] == 1 & data[, m] == 0] = ((1 - p0)/(1 - p1))[data[, t] == 1 & data[, m] == 0]
      data$rmpw[data[, t] == 0] = 1
      data$p1 = p1 # For the use in the function of "est"
      data$p0 = p0

      return(data)
    }
    
    #### Direct and Indirect Effect Estimation
    est = function(y, m, t, X, data) {
      data = weight0(m, t, X, data)
      n = nrow(data)
      data = data[order(data[, t], decreasing = T), ]
      muAlpha = data$p0
      muBeta = data$p1
      #for anything related to propensity model estimation, will need to use the orginal muAlpha, muBeta
      muAlphaOrig = muAlpha
      muBetaOrig = muBeta
      x = matrix(1, n, 1) 
      x = cbind(x, as.matrix(data[, X]))
      nP = ncol(x)
      dimnames(x)[[2]] <- c( "1", X)
      alphaVar = (muAlphaOrig * (1 - muAlphaOrig))
      betaVar =  (muBetaOrig * (1 - muBetaOrig))
      alphaResidue = (data[, m] - muAlphaOrig) * (data[, t] == 0)
      betaResidue = (data[, m] - muBetaOrig) * (data[, t] == 1)
      
      G11 = matrix(0, 2 * nP, 2 * nP)
      G11_Alpha = t(x) %*% diag(alphaVar * (data[, t] == 0)) %*% x
      G11_Beta = t(x) %*% diag(betaVar * (data[, t]==1)) %*% x
      G11[1:nP, 1:nP] = G11_Alpha
      G11[(nP+1):(2*nP), (nP+1):(2*nP)] = G11_Beta
      
      w = matrix(0, n, 3)
      w[,1] = (data[, t] == 0)
      w[,3] = (data[, t] == 1)
      w[,2] = (data[, t] == 1) * data$rmpw
      deltaEstimate = NULL
      wSum = NULL
      hDelta = NULL 
      for (j in 1:3){
        delta = sum(data[, y] * w[,j])/sum(w[, j])    
        hDelta = cbind(hDelta, (data[, y] - delta) * w[, j])
        deltaEstimate = cbind(deltaEstimate, delta)
        wSum = cbind(wSum, sum(w[, j]))
      }
      hAlphaBeta = cbind(x * (alphaResidue), x * (betaResidue))
      hCombined = cbind(hAlphaBeta, hDelta);
      B0 = t(hCombined) %*% hCombined
      G22 = diag(as.numeric(wSum))
      
      dWi_dAlpha_noX = ((data[, t] == 1) * muAlpha * (1-muAlpha) * (data[, m]/muBeta - (1-data[, m])/(1-muBeta))) 
      dWi_dBeta_noX = ((data[, t] == 1) * (-data[, m] * muAlpha * (1 - muBeta)/muBeta + (1 - data[, m]) * (1 - muAlpha) * muBeta/(1 - muBeta))) 
      dhDelta_dAlpha = (-dWi_dAlpha_noX * (data[, y] - deltaEstimate[2])) * x
      dhDelta_dBeta = (-dWi_dBeta_noX * (data[, y] - deltaEstimate[2])) * x
      G12 = matrix(0, 3, 2 * nP)
      for (j in 1:nP) {
        G12[2, j] = sum(dhDelta_dAlpha[, j])
        G12[2, j + nP] = sum(dhDelta_dBeta[, j])
      }
      
      A0 = matrix(0, (nP * 2 + 3), (nP * 2 + 3))
      A0[(2 * nP + 1):(2 * nP + 3), (2 * nP + 1):(2 * nP + 3)] = G22
      A0[1:(2 * nP), 1:(2 * nP)] = G11
      A0[(2 * nP + 1):(2 * nP + 3), 1:(2 * nP)] = G12
      
      v_hw = solve(A0) %*% B0 %*% t(solve(A0))
      v_hw_ctrl_counterfacal_tr = v_hw[(2 * nP + 1):(2 * nP + 3), (2 * nP + 1):(2 * nP + 3)]  
      de = deltaEstimate[2] - deltaEstimate[1]
      ie = deltaEstimate[3] - deltaEstimate[2]
      se_de = sqrt(v_hw_ctrl_counterfacal_tr[2, 2] + v_hw_ctrl_counterfacal_tr[1, 1] - 2 * v_hw_ctrl_counterfacal_tr[1, 2])
      se_ie = sqrt(v_hw_ctrl_counterfacal_tr[3, 3] + v_hw_ctrl_counterfacal_tr[2, 2] - 2 * v_hw_ctrl_counterfacal_tr[2, 3])
      
      results = c(de = de, ie = ie, se_de = se_de, se_ie = se_ie, CIL_de = de - 1.96 * se_de, CIU_de = de + 1.96 * se_de, CIL_ie = ie - 1.96 * se_ie, CIU_ie = ie + 1.96 * se_ie)
      
      return(results)
    }
    
    result_now = est(outcome, mediator, treatment, propensity_x, data)
    result = cbind(c(result_now["de"], result_now["ie"]),c(result_now["se_de"], result_now["se_ie"]))
    z = result[, 1]/result[, 2]
    p = NULL
    for(i in 1:nrow(result)){
      p = c(p, (1 - pnorm(abs(z[i]))) * 2)
    }
    result = round(cbind(result, z, p), 4)
    result[p < 0.001, 4] = "<0.001"
    sig = NULL
    sig[p <= 0.001] = "**"
    sig[p > 0.001 & p <= 0.01] = "*"
    sig[p > 0.01 & p <= 0.05] = "."
    sig[p > 0.05] = ""
    result = cbind(result, sig)
    result = as.data.frame(result)
    colnames(result) <- c("Estimate", "Std.Error", "t value", "Pr(>|t|)", "")
    rownames(result) <- c("Natural Direct Effect", "Natural Indirect Effect")
  }
  
  if(decomposition > 0){
    weight1 = function(input_data, treatment, propensity_yx, decomposition = 1) {
      data_tr = input_data[input_data[, treatment] == 1, ]
      data_ctrl = input_data[input_data[, treatment] == 0, ]
      
      #create the model formula based on propensity_yx
      nM = length(propensity_yx);
      fmlaString = paste(propensity_yx[1], "~", paste(propensity_yx[2:nM], collapse="+"))
      fmla <- as.formula(fmlaString);
      
      ##Estimated logit of propensity score
      {
        l_tr = glm(fmla, data = data_tr, family = binomial)
        propensity_l_tr_data_tr = fitted(l_tr)
        logit_l_tr_data_tr = log(propensity_l_tr_data_tr/(1 - propensity_l_tr_data_tr))
        l_ctrl = glm(fmla, data = data_ctrl, family = binomial)
        propensity_l_ctrl_data_ctrl = fitted(l_ctrl)
        logit_l_ctrl_data_ctrl = log(propensity_l_ctrl_data_ctrl/(1 - propensity_l_ctrl_data_ctrl))
        
        logit_l_tr_data_ctrl = predict(l_tr, data_ctrl)
        logit_l_ctrl_data_tr = predict(l_ctrl, data_tr)
      }
      
      logit_l_tr = c(logit_l_tr_data_tr, logit_l_tr_data_ctrl)
      logit_l_ctrl = c(logit_l_ctrl_data_tr, logit_l_ctrl_data_ctrl)
      propensity_l_tr = exp(logit_l_tr)/(exp(logit_l_tr) + 1)#This is namely phi_1  in the memo.
      propensity_l_ctrl = exp(logit_l_ctrl)/(exp(logit_l_ctrl) + 1)#This is namely phi_0 in the memo.
      
      data = rbind(data_tr, data_ctrl)
      data = cbind(data, propensity_l_tr, propensity_l_ctrl)
      
      #### RMPW weight
      w_tr_me1 = (data$propensity_l_ctrl/data$propensity_l_tr)[data[, treatment] == 1 & data[, mediator] == 1]
      w_tr_me0 = ((1 - data$propensity_l_ctrl)/(1 - data$propensity_l_tr))[data[, treatment] == 1 & data[, mediator] == 0]
      w_ctrl_me1 = (data$propensity_l_tr/data$propensity_l_ctrl)[data[, treatment] == 0 & data[, mediator] == 1]
      w_ctrl_me0 = ((1 - data$propensity_l_tr)/(1 - data$propensity_l_ctrl))[data[, treatment] == 0 & data[, mediator] == 0]
      
      rmpw = c(w_tr_me1, w_tr_me0, w_ctrl_me1, w_ctrl_me0)
      
      tr_me1 = data[data[, treatment] == 1 & data[, mediator] == 1, ]
      tr_me0 = data[data[, treatment] == 1 & data[, mediator] == 0, ]
      ctrl_me1 = data[data[, treatment] == 0 & data[, mediator] == 1, ]
      ctrl_me0 = data[data[, treatment] == 0 & data[, mediator] == 0, ]
      data = rbind(tr_me1, tr_me0, ctrl_me1, ctrl_me0)
      data = cbind(rmpw, data)  
      data_nodup = data 
      
      #then create data set with duplications
      newdata = data
      newdata.ctrl=newdata[which(newdata[, treatment]==0),]
      newdata.tr=newdata[which(newdata[, treatment]==1),]
      
      if(decomposition == 1){
        
        d1_rmpw = c(rep(c(0,1),dim(newdata.tr)[1]),  rep(c(0, 0),dim(newdata.ctrl)[1]))
        d0_rmpw = c(rep(c(0,0),dim(newdata.tr)[1]),  rep(c(0, 1),dim(newdata.ctrl)[1]))
        
        newdata.tr.dup=NULL
        for(j in 1:dim(newdata.tr)[2]){
          newdata.tr.dup=cbind(newdata.tr.dup,rep(newdata.tr[,j],rep(2,dim(newdata.tr)[1])))
        }
        colnames(newdata.tr.dup)=colnames(newdata.ctrl)
        newdata.tr.dup=as.data.frame(newdata.tr.dup)
        newdata.tr.dup$rmpw[seq(2,length(newdata.tr.dup[,1]),2)]=1
        
        newdata.ctrl.dup=NULL
        for(j in 1:dim(newdata.ctrl)[2]){
          newdata.ctrl.dup=cbind(newdata.ctrl.dup,rep(newdata.ctrl[,j],rep(2,dim(newdata.ctrl)[1])))
        }
        colnames(newdata.ctrl.dup)=colnames(newdata.ctrl)
        newdata.ctrl.dup=as.data.frame(newdata.ctrl.dup)
        newdata.ctrl.dup$rmpw[seq(1,length(newdata.ctrl.dup[,1]),2)]=1
        
        data_dup=cbind(d1_rmpw,d0_rmpw, rbind(newdata.tr.dup,newdata.ctrl.dup))  
        result = list(data_dup= data_dup, data_nodup = data_nodup)
      } else {
        
        d1_rmpw = c(rep(c(0,1),dim(newdata.tr)[1]),  rep(c(0, 0),dim(newdata.ctrl)[1]))
        d0_rmpw = c(rep(c(1,0),dim(newdata.tr)[1]),  rep(c(0, 1),dim(newdata.ctrl)[1]))
        
        newdata.tr.dup=NULL
        for(j in 1:dim(newdata.tr)[2]){
          newdata.tr.dup=cbind(newdata.tr.dup,rep(newdata.tr[,j],rep(2,dim(newdata.tr)[1])))
        }
        colnames(newdata.tr.dup)=colnames(newdata.ctrl)
        newdata.tr.dup=as.data.frame(newdata.tr.dup)
        newdata.tr.dup$rmpw[seq(1,length(newdata.tr.dup[,1]),2)]=1
        
        newdata.ctrl.dup=NULL
        for(j in 1:dim(newdata.ctrl)[2]){
          newdata.ctrl.dup=cbind(newdata.ctrl.dup,rep(newdata.ctrl[,j],rep(2,dim(newdata.ctrl)[1])))
        }
        colnames(newdata.ctrl.dup)=colnames(newdata.ctrl)
        newdata.ctrl.dup=as.data.frame(newdata.ctrl.dup)
        newdata.ctrl.dup$rmpw[seq(1,length(newdata.ctrl.dup[,1]),2)]=1
        
        data_dup=cbind(d1_rmpw,d0_rmpw, rbind(newdata.tr.dup,newdata.ctrl.dup))
        
        #add two intermediate variable
        data_dup$trd0_rmpw = data_dup[, treatment]*data_dup$d0_rmpw
        data_dup$trd1_rmpw = data_dup[, treatment]*data_dup$d1_rmpw
        
        result = list(data_dup= data_dup, data_nodup = data_nodup)    
      }
      
      return(result)
    }
    
    
    est_rmpw = function(data_dup, data_nodup, treatment, outcome_yx, propensity_yx, decomposition = 1) {
      
      
      #at first we estimate the parameters
      
      fmlaString = paste(outcome_yx[1], "~", paste(outcome_yx[2:length(outcome_yx)], collapse="+"))
      outcomeFmla <- as.formula(fmlaString);
      
      l = lm(outcomeFmla, data=data_dup, weights=rmpw)
      beta = as.numeric(coef(l))
      
      nG = 4 #experiment, control, and their counterfactual groups, totally 4
      nP = length(outcome_yx)-1 - (nG-1)
      
      m = length(propensity_yx)-1
      nM = m+1
      
      propensity_x = propensity_yx[2:nM]
      outcome = outcome_yx[1]
      
      #then let's estimate the covariance matrix between those coefficients
      data = data_nodup
      n = length(data[, 1])
      n_tr = sum(data[, treatment] == 1)
      n_ctrl = sum(data[, treatment] == 0)
      
      phi_0 = data$propensity_l_ctrl
      phi_1 = data$propensity_l_tr
      
      #predictors for propensity model
      x = matrix(1,n,1) #constant term, 
      x = cbind(x, as.matrix(data[,propensity_x]))
      
      #either Y or Y minus the part contributed by the covariates in the outcome model
      useY = data[,outcome] 
      
      if(nP > 0)
      {
        outcome_xo = outcome_yx[(nG+1):length(outcome_yx)]
        xo = as.matrix(data[, outcome_xo])
        useY =  data[,outcome] - xo%*%beta[(nG+1):length(outcome_yx)]
      }
      
      #one by one, I follow the memo to implement it.
      
      s_0_noX = (data[, mediator] - phi_0)*(1-data[, treatment])
      s_1_noX = (data[, mediator] - phi_1)*data[, treatment]
      
      h1 = cbind(x*s_0_noX, x*s_1_noX)
      
      #elements in beta are gamma_0, gamma_DE.0, gamma_
      gamma_0 = beta[1]
      
      if(decomposition == 1) {
        gamma_DE.0 = beta[2]
        gamma_IE.1 = beta[3]
        gamma_IE.0 = beta[4]
        
        #will define mean of the four groups under frist decomposition method
        gamma_star.0 = gamma_0 + gamma_IE.0
        gamma_star.1 = gamma_0 + gamma_DE.0
        gamma_1 = gamma_0 + gamma_DE.0 + gamma_IE.1
      } else {
        
        gamma_IE.0 = beta[2]
        gamma_DE.1 = beta[3]
        gamma_DE.0 = beta[4]
        
        #will define mean of the four groups under frist decomposition method
        gamma_star.0 = gamma_0 + gamma_IE.0
        gamma_star.1 = gamma_0 + gamma_DE.0
        gamma_1 = gamma_0 + gamma_IE.0 + gamma_DE.1
      }
      
      gammaEstimate = cbind(gamma_0, gamma_star.0, gamma_star.1, gamma_1)
      dimnames(gammaEstimate)[[2]] <- c( "gamma_0","gamma_*0", "gamma_*1","gamma_1")
      
      m_weight = matrix(1, n, nG)
      m_weight[,1] = 1-data[, treatment]
      m_weight[,2] = (1-data[, treatment])*data$rmpw
      m_weight[,3] = data[, treatment]*data$rmpw
      m_weight[,4] = data[, treatment]
      
      
      h2 = matrix(0, n, nG)
      for (i in 1:nG) {
        h2[,i] = (useY-gammaEstimate[i])*m_weight[,i]
      }
      
      h = cbind(h1, h2)
      
      
      if(nP > 0)
      {
        i=1
        h3_noX = h2[,i]
        for (i in 2:nG) {
          h3_noX = h3_noX + h2[,i]
        }
        h3 = xo*h3_noX    
        h = cbind(h, h3)
      }
      
      B0 = t(h)%*%h
      
      #then we calculate A0 matrix
      G11 = matrix(0, 2*nM, 2*nM)
      G11_0 = t(x)%*%diag(-phi_0*(1-phi_0)*(1-data[, treatment]))%*%x
      G11_1 = t(x)%*%diag(-phi_1*(1-phi_1)*data[, treatment])%*%x
      G11[1:nM, 1:nM] = G11_0
      G11[(nM+1):(2*nM), (nM+1):(2*nM)] = G11_1
      
      diag_G22 = NULL
      for (i in 1:nG) {
        diag_G22 = cbind(diag_G22, -sum(m_weight[,i]))    
      }  
      G22 = diag(as.vector(diag_G22))
      
      G21 = matrix(0, nG, 2*nM)
      G21_row2_noWeight_noX = (useY-gammaEstimate[2])*(1-data[, treatment])
      G21_row2_part1_noX =  G21_row2_noWeight_noX*(-phi_1*(1-phi_0)/phi_0*data[, mediator] + (1-phi_1)*phi_0/(1-phi_0)*(1-data[, mediator]))
      G21_row2_part2_noX =  G21_row2_noWeight_noX*(1/phi_0*data[, mediator] - 1/(1-phi_0)*(1-data[, mediator]))*phi_1*(1-phi_1)
      
      G21_row3_noWeight_noX = (useY-gammaEstimate[3])*data[, treatment]  
      G21_row3_part1_noX = G21_row3_noWeight_noX*(1/phi_1*data[, mediator] - 1/(1-phi_1)*(1-data[, mediator]))*phi_0*(1-phi_0)
      G21_row3_part2_noX = G21_row3_noWeight_noX*(-phi_0*(1-phi_1)/phi_1*data[, mediator] + (1-phi_0)*phi_1/(1-phi_1)*(1-data[, mediator]))
      
      for (j in 1:nM) {
        G21[2, j] = sum(G21_row2_part1_noX*x[,j])
        G21[2, j+nM] =  sum(G21_row2_part2_noX*x[,j])
        G21[3, j] = sum(G21_row3_part1_noX*x[,j])
        G21[3, j+nM] =  sum(G21_row3_part2_noX*x[,j])
      }
      
      A0 = matrix(0, (nM*2+nG+nP), (nM*2+nG+nP))
      A0[1:(2*nM), 1:(2*nM)] = G11
      A0[(2*nM+1):(2*nM+nG), (2*nM+1):(2*nM+nG)] = G22
      A0[(2*nM+1):(2*nM+nG), 1:(2*nM)] = G21
      
      if(nP > 0)
      {
        
        i=1
        #pay attention to the negative sign here
        G33_noX = - m_weight[,i]
        #pay attention the index i is in 1:nG
        for (i in 2:nG) {
          #pay attention to the negative sign here
          G33_noX = G33_noX - m_weight[,i]
        }
        G33 = t(xo)%*%diag(as.vector(G33_noX))%*%xo
        
        
        G31 = matrix(0, nP, 2*nM)
        
        G31_part1_noX = ( (useY-gammaEstimate[2])*(1-data[, treatment])*(-phi_1*(1-phi_0)/phi_0*data[, mediator] + (1-phi_1)*phi_0/(1-phi_0)*(1-data[, mediator])) +
                            (useY-gammaEstimate[3])*data[, treatment]*(1/phi_1*data[, mediator] - 1/(1-phi_1)*(1-data[, mediator]))*phi_0*(1-phi_0))
        
        G31_part2_noX = ( (useY-gammaEstimate[2])*(1-data[, treatment])*(1/phi_0*data[, mediator] - 1/(1-phi_0)*(1-data[, mediator]))*phi_1*(1-phi_1) +
                            (useY-gammaEstimate[3])*data[, treatment]*(-phi_0*(1-phi_1)/phi_1*data[, mediator] + (1-phi_0)*phi_1/(1-phi_1)*(1-data[, mediator])))
        
        G31[,1:nM] =  t(xo)%*%diag(as.vector(G31_part1_noX))%*%x  
        G31[,(nM+1):(2*nM)] = t(xo)%*%diag(as.vector(G31_part2_noX))%*%x  
        
        
        #G32
        G23 = matrix(0, nG, nP)
        
        for(i in 1:nG){
          #pay attention to the negative sign here
          for(j in 1:nP) G23[i,j] = - sum(m_weight[,i]*xo[,j])      
        }
        
        #then we fill them in A0
        A0[(2*nM+nG+1):(2*nM+nG+nP), 1:(2*nM)] = G31
        A0[(2*nM+nG+1):(2*nM+nG+nP), (2*nM+1):(2*nM+nG)] = t(G23)
        A0[(2*nM+1):(2*nM+nG), (2*nM+nG+1):(2*nM+nG+nP)] = G23
        A0[(2*nM+nG+1):(2*nM+nG+nP), (2*nM+nG+1):(2*nM+nG+nP)] = G33
      }
      
      #now we calculate the covariance matrix
      CMatrix = solve(A0) %*% B0 %*% t(solve(A0))
      
      #covariance matrix for ( "gamma_0","gamma_*0", "gamma_*1","gamma_1")
      C22 = CMatrix[(2*nM+1):(2*nM+nG), (2*nM+1):(2*nM+nG)]
      
      if(decomposition == 1) {
        #from C22 we will calculate the covariance matrix between (gamma_0, DE.0, IE.1, IE0, and IE.1-IE.0)
        convertMatrix = matrix(0, nG, nG+1)  
        convertMatrix[,1] = c(1, 0, 0, 0)
        convertMatrix[,2] = c(-1, 0, 1, 0)
        convertMatrix[,3] = c(0, 0, -1, 1)
        convertMatrix[,4] = c(-1, 1, 0, 0)
        convertMatrix[,5] = c(1, -1, -1, 1)
      } else {
        #from C22 we will calculate the covariance matrix between (gamma_0, IE.0, DE.1, DE0, and DE.1-DE.0)
        convertMatrix = matrix(0, nG, nG+1)  
        convertMatrix[,1] = c(1, 0, 0, 0)
        convertMatrix[,2] = c(-1, 1, 0, 0)
        convertMatrix[,3] = c(0, -1, 0, 1)
        convertMatrix[,4] = c(-1, 0, 1, 0)
        convertMatrix[,5] = c(1, -1, -1, 1)
      }
      
      gammaCov = t(convertMatrix)%*%C22%*%convertMatrix
      gammaSE = sqrt(diag(gammaCov))
      
      #covariance matrix for the coefficients of the covariate in the outcome model
      if(nP > 0){
        C33 = CMatrix[(2*nM+nG+1):(2*nM+nG+nP), (2*nM+nG+1):(2*nM+nG+nP)]
        lamdaCov = C33
        if(nP > 1) lamdaSE = sqrt(diag(lamdaCov)) else lamdaSE = sqrt(lamdaCov)
        coefValue = c(beta[1:nG], beta[3]-beta[4], beta[(nG+1):(nG+nP)])  
        coefSE = c(gammaSE, lamdaSE)
        if(decomposition == 1){
          resultLabel = c("Gamma.0", "Natural Direct Effect", "Natural Indirect Effect", "Pure Indirect Effect", "T-by-M Interaction Effect", outcome_x) 
        } else{
          resultLabel = c("Gamma.0", "Pure Indirect Effect", "Total Direct Effect", "Natural Direct Effect", "T-by-M Interaction Effect", outcome_x)
        }
      } else {    
        coefValue = c(beta[1:nG], beta[3]-beta[4])
        coefSE = gammaSE
        if(decomposition == 1){
          resultLabel = c("Gamma.0", "Natural Direct Effect", "Natural Indirect Effect", "Pure Indirect Effect", "T-by-M Interaction Effect") 
        } else{
          resultLabel = c("Gamma.0", "Pure Indirect Effect", "Total Direct Effect", "Natural Direct Effect", "T-by-M Interaction Effect")
        }
      }
      
      result = list(coefValue=coefValue, coefSE=coefSE, resultLabel=resultLabel)
      return(result)
    }
    
    propensity_yx = c(mediator, propensity_x) #response and covariates for propensity score models
    rmpw_list = weight1(data, treatment, propensity_yx, decomposition)
    data_nodup = rmpw_list$data_nodup
    data_dup = rmpw_list$data_dup
    
    outcome_yx = c(outcome, treatment, "d1_rmpw", "d0_rmpw");
    outcome_yx = c(outcome_yx, outcome_x)
    result_now = est_rmpw(data_dup, data_nodup, treatment, outcome_yx, propensity_yx, decomposition)
    result = t(rbind(result_now$coefValue, result_now$coefSE))
    z = result[, 1]/result[, 2]
    p = NULL
    for(i in 1:nrow(result)){
      p = c(p, (1 - pnorm(abs(z[i]))) * 2)
    }
    result = round(cbind(result, z, p), 4)
    result[p < 0.001, 4] = "<0.001"
    sig = NULL
    sig[p <= 0.001] = "**"
    sig[p > 0.001 & p <= 0.01] = "*"
    sig[p > 0.01 & p <= 0.05] = "."
    sig[p > 0.05] = ""
    result = cbind(result, sig)
    result = as.data.frame(result)
    rownames(result) <- result_now$resultLabel
    colnames(result) <- c("Estimate", "Std.Error", "t value", "Pr(>|t|)", "")
  }
  
  return(result)
}
