#One stage competing risks AFT model
library(lme4)
library(tidyverse)
library(ggpubr)
library(geepack)
library(simstudy)

expit = function(x) exp(x)/(1 + exp(x))

aft = function(data, treat.mod, cens.mod, outcome.mod, cause){
  
  #Treatment model and censoring model 
  
  treat_mod = glm(treat.mod, data = data, family = binomial)
  cens_mod = glm(cens.mod, data = data, family = binomial)
  
  treat_prob = predict(treat_mod, type = "response")
  cens_prob = predict(cens_mod, type = "response")

  #Weights are wrong for censored observations but doesent matter since they are excluded
  weights = abs(data$a - treat_prob)/cens_prob
  
  #Estimation of blip parameters via random effects AFT model
  index = with(data, epsilon == cause & delta == 1)
  data = data[index,]
  weights = weights[index]
  
  #Fitting linear mixed model
  
  model = lmer(formula = outcome.mod, data = data, weights = weights)
  fixef(model)
  
}

aft2 = function(data, treat.mod, cens.mod, outcome.mod, cause){
  
  #Treatment model and censoring model 
  
  treat_mod = glm(treat.mod, data = data, family = binomial)
  cens_mod = glm(cens.mod, data = data, family = binomial)
  
  treat_prob = predict(treat_mod, type = "response")
  cens_prob = predict(cens_mod, type = "response")
  
  #Weights are wrong for censored observations but doesent matter since they are excluded
  weights = abs(data$a - treat_prob)/cens_prob
  
  #Estimation of blip parameters via random effects AFT model
  index = with(data, epsilon == cause & delta == 1)
  data = data[index,]
  w = weights[index]
  
  data = data %>% mutate(w = w)
  
  #Fitting GEE with exchangeable correlation structure
  
  data = data %>% arrange(group)
  model = geeglm(formula = outcome.mod, data = data, weights = w, id = group, corstr = "exchangeable")
  coef(model)
}

sim_data = function(n, cens = T, psi1_fixed, psi2_fixed){
  #Direct approach
  x <- rnorm(n, 0, 2)
  z = rnorm(n, 0, 2)
  a <- rbinom(n, 1, expit(0.5 + x +z))
  
  #Generate epsilon
  epsilon = rbinom(n,1,expit(x+0.5))+1
  
  
  #Generate groups and random intercepts for each group
  group = sample(50, n, replace = T)
  rand_effects = rnorm(50, mean = 0, sd = 1)
  ind_effects = map_dbl(group, function(x) rand_effects[x])
  
  beta1 <- c(1, 1.5, 2)
  psi1 <- psi1_fixed
  h1beta <- model.matrix(~x + z)
  h1psi <- model.matrix(~x)
  T1 <- exp(h1beta[epsilon== 1, ] %*% beta1+ a[epsilon == 1] *h1psi[epsilon ==1, ] %*% psi1 + ind_effects[epsilon ==1]+ rnorm(sum(epsilon == 1),sd = 0.3))
  
  beta2 <- c(1, 3,0.2)
  psi2<- psi2_fixed
  h2beta <- model.matrix(~x + z)
  h2psi <- model.matrix(~x)
  T2 <- exp(h2beta[epsilon == 2, ] %*% beta2+ a[epsilon == 2] *h2psi[epsilon ==2, ] %*% psi2 + ind_effects[epsilon ==2] +rnorm(sum(epsilon == 2),sd = 0.3))
  
  
  #Construct dataframe
  Y = rep(0,n)
  Y[epsilon == 1] = T1
  Y[epsilon == 2] = T2
  
  group = as.factor(group)
  data = tibble(x,z,a,group,Y,epsilon)
  
  #Replace times by censoring time for those who were censored
  if(cens){  
    delta <- rbinom(n, 1, expit(-x + 1.9-0.3*z))
    C <- rexp(sum(delta == 0), 1/300)
    Y[delta == 0] = C
    data = data %>% mutate(delta = delta)
  }
  
  data
}

sim_stochastic_regime = function(train, test, psi1_fixed, psi2_fixed){
  
  #Generate training data and obtain blip estimates for both causes 
  psi1_hat = aft(train, a~x+z, delta~x+z, log(Y)~ x+z+ a+a:x + (1|group), 1)[c(4,5)]
  psi2_hat = aft(train, a~x+z, delta~x+z, log(Y)~ x+z+ a+a:x + (1|group), 2)[c(4,5)]
  
  #Model for probability of failing from cause 1
  cause_data = train %>% filter(delta == 1) %>% mutate(ind = as.numeric(epsilon ==1))
  cause_mod = glm(ind~ x+z, data = cause_data, family = binomial)
  
  #Get the optimal regime depending on true underlying cause
  psi1 = psi1_fixed
  psi2 = psi2_fixed
  #Note : have to make sure this works the way i think it does
  n_test = nrow(test)
  test = test %>% mutate(blip = ifelse(epsilon ==1, cbind(1,x) %*% psi1_fixed, cbind(1,x) %*% psi2_fixed),opt = as.numeric(blip > 0))
  
  #Allocate treatment according to weighted rule
  test  = test %>% mutate(cause_prob = predict(cause_mod,newdata = test, type = "response"), blip_weighted = (cbind(1,x) %*% psi1_hat)*cause_prob + (cbind(1,x) %*% psi2_hat)*(1-cause_prob),
                          opt_weighted = as.numeric(blip_weighted > 0))
  
  #Allocate treatment according to greedy rule
  
  test = test %>% mutate(blip_greedy = ifelse(cause_prob > 0.5,cbind(1,x) %*% psi1_hat,cbind(1,x) %*% psi2_hat), opt_greedy = as.numeric(blip_greedy > 0))
  
  #Identify proportion of optimal treatment (POT)
  pot_w = sum(test$opt == test$opt_weighted)/n_test
  pot_g = sum(test$opt == test$opt_greedy)/n_test
  
  #Compute mean blip evaluated at true cause (can split this into two different metrics if not interpretable)
  mean_blip_w = mean(test$blip * test$opt_weighted - test$blip* (1-test$opt_weighted))
  mean_blip_g = mean(test$blip * test$opt_greedy - test$blip* (1-test$opt_weighted))
  
  c(pot_w,mean_blip_w,test$blip_weighted, pot_g,mean_blip_g, test$blip_greedy)
  
}

plot_results = function(res, test){
  
  #Plot results for all 3 regimes
  
  test = test %>% mutate(blip = ifelse(epsilon ==1, cbind(1,x) %*% psi1_fixed, cbind(1,x) %*% psi2_fixed))
  #Weighted regime
  
  weighted_data = res[,1:(n_test+2)]
  
  weighted_pot = weighted_data[,1]
  weighted_blip = weighted_data[,2]
  
  blip_quantiles = t(apply(weighted_data[,3:(n_test+2)], 2, function(x) quantile(x, probs = c(0.025,0.5, 0.975)))) %>% as_tibble() 
  colnames(blip_quantiles) <- c("lower", "mid", "upper")
  
  
  dat = bind_cols(test, blip_quantiles) %>% arrange(blip)
  
  dat1 = dat %>% filter(epsilon == 1)
  dat2 = dat %>% filter(epsilon == 2)
  
  #IDEA : FOR EACH CAUSE, CAN PLOT BOTH STRATEGIES AT ONCE IN DIFFERENT COLOURS (ONLY PROBLEM IS 95\% CI)
  # Plotting
  
  colors <- c("Chosen Strategy" = "black", "True Blip" = "red")
  
  w.plot1 =ggplot(dat1, aes(x = 1:nrow(dat1))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = mid, color = "Chosen Strategy"), linewidth= 1.2)  + geom_line(aes(y = blip, color = "True Blip"), linewidth = 1.2)+ labs(x="Observation number (Ordered)", y = "Median Difference in Log Survival Time",title ="Cause 1", caption = "(a) Weighted Strategy", color = "Legend") + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  w.plot2 =ggplot(dat2, aes(x = 1:nrow(dat2))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = mid,color = "Chosen Strategy"), linewidth= 1.2)  + geom_line(aes(y = blip, color = "True Blip"), linewidth = 1.2)+ labs(x="Observation number (Ordered)", y = "Median Difference in Log Survival Time",title ="Cause 2", caption = "(c) Weighted Strategy",color = "Legend") + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  
  #Greedy regime
  greedy_data = res[,(n_test+3):(2*n_test+4)]
  
  greedy_pot = greedy_data[,1]
  greedy_blip = greedy_data[,2]
  
  blip_quantiles = t(apply(greedy_data[,3:(n_test+2)], 2, function(x) quantile(x, probs = c(0.025,0.5, 0.975)))) %>% as_tibble() 
  colnames(blip_quantiles) <- c("lower", "mid", "upper")
  
  dat = bind_cols(test, blip_quantiles) %>% arrange(blip)
  greedy_blip = mean(dat$blip)
  
  dat1 = dat %>% filter(epsilon == 1)
  dat2 = dat %>% filter(epsilon == 2)
  
  #Plotting
  
  g.plot1 =ggplot(dat1, aes(x = 1:nrow(dat1))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = mid, color = "Chosen Strategy"), linewidth= 1.2)  + geom_line(aes(y = blip, color = "True Blip"), linewidth = 1.2)+labs(x="Observation number (Ordered)", y = "Median Difference in Log Survival Time", title ="Cause 1", caption = "(b) Greedy Strategy", color = "Legend") + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  g.plot2 =ggplot(dat2, aes(x = 1:nrow(dat2))) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70")+ geom_line(aes(y = mid,color = "Chosen Strategy"), linewidth= 1.2)  + geom_line(aes(y = blip, color = "True Blip"), linewidth = 1.2)+ labs(x="Observation number (Ordered)", y = "Median Difference in Log Survival Time", title ="Cause 2", caption = "(d) Greedy Strategy",color = "Legend") + theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0.5))+scale_color_manual(values = colors)
  
  #POT 
  pot = cbind(weighted_pot, greedy_pot)
  
  #Mean blip
  mean_blip = cbind(weighted_blip, greedy_blip)
  means =colMeans(mean_blip)
  
  #Arrange plots
  blip_plot=ggarrange(w.plot1, g.plot1,w.plot2, g.plot2, ncol = 2, nrow = 2, common.legend = T, legend = "right")
  
  #Return all plots
  list(blip_plot, pot, means)
}