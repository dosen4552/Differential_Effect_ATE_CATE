library('MASS')
library('nnet')
library('EnvStats')
library('matrixStats')
library('parallel')
# New bound ATE functions
bound_ATE_new = function(dataset, alpha = 0.05){
  Z1 = data.matrix(dataset[,'Treatment']) 
  Z2 = data.matrix(dataset[,'Z']) 
  Y = data.matrix(dataset[,'Y'])
  d = dim(dataset)[2] - 5
  X = data.matrix(dataset[,5:(5+d-1)]) 
  # define factor 00 mapped to 0, 01 mapped to 1, 10 mapped to 2, 11 mapped to 3
  Z1Z2 = factor(Z1 * 2 + Z2)
  multinomial_model <- multinom(Z1Z2 ~ X)
  propens <- fitted(multinomial_model)
  propens01 <- propens[,"1"]
  propens10 <- propens[,"2"]
 
  mu_p = mean((Z1 == 1 & Z2 == 0) * Y/propens10) - mean((Z1 == 0 & Z2 == 1) * Y / propens01  )
  propens = glm(Z1 ~ X + Z2, family = 'binomial')$fitted.value
  mu_m =  mean(Z1 * Y/propens) - mean((1-Z1) * Y / (1 - propens)  )
  sigma_p = sqrt(sum( ((Z1 == 1 & Z2 == 0) * Y/propens10 - (Z1 == 0 & Z2 == 1) * Y / propens01 - mu_p)^2 )/(dim(dataset)[1]))
  sigma_m = sqrt(sum( (Z1 * Y/propens - (1-Z1) * Y / (1 - propens) - mu_m)^2 )/(dim(dataset)[1]))
  if(mu_p > mu_m){
    tau_u = mu_p
    sigma_u = sigma_p
    tau_l = mu_m
    sigma_l = sigma_m
  }else{
    tau_l = mu_p
    sigma_l = sigma_p
    tau_u = mu_m
    sigma_u = sigma_m
  }
  lower = c(tau_l, tau_l - qnorm(1-alpha/2)*sigma_l, tau_l + qnorm(1-alpha/2)*sigma_l  )
  upper = c(tau_u, tau_u - qnorm(1-alpha/2)*sigma_u, tau_u + qnorm(1-alpha/2)*sigma_u  )
  result = rbind(lower, upper)
  colnames(result) = c('estimate',paste0((1-alpha/2)*100, '% lower'),paste0((1-alpha/2)*100, '% upper'))
  rownames(result) = c('lower bound','upper bound')
  return(result)
}

bound_ATE_new_bootstrap = function(dataset, alpha = 0.05, num_bootstrap = 1000, beta = 0.005){
  mu_1 = rep(0, num_bootstrap)
  mu_2 = rep(0, num_bootstrap)
  
  Z1 = data.matrix(dataset[,'Treatment']) 
  Z2 = data.matrix(dataset[,'Z']) 
  Y = data.matrix(dataset[,'Y'])
  d = dim(dataset)[2] - 5
  X = data.matrix(dataset[,5:(5+d-1)]) 
  # define factor 00 mapped to 0, 01 mapped to 1, 10 mapped to 2, 11 mapped to 3
  Z1Z2 = factor(Z1 * 2 + Z2)
  invisible(capture.output(multinomial_model <- nnet::multinom(Z1Z2 ~ X), file = NULL))
  propens <- fitted(multinomial_model)
  propens01 <- propens[,"1"]
  propens10 <- propens[,"2"]
  
  mu_p_s = mean((Z1 == 1 & Z2 == 0) * Y/propens10) - mean((Z1 == 0 & Z2 == 1) * Y / propens01  )
  propens = glm(Z1 ~ X + Z2, family = 'binomial')$fitted.value
  mu_m_s =  mean(Z1 * Y/propens) - mean((1-Z1) * Y / (1 - propens)  )
  
  # Set up the cluster
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  
  # Define the bootstrap function
  bootstrap_function <- function(i) {
    valid_sample <- FALSE
    while (!valid_sample) {
      # Generate bootstrap sample
      temp <- dataset[sample(1:dim(dataset)[1], dim(dataset)[1], replace = TRUE),]
      Z1 <- data.matrix(temp[,'Treatment'])
      Z2 <- data.matrix(temp[,'Z'])
      Y <- data.matrix(temp[,'Y'])
      d <- dim(temp)[2] - 5
      X <- data.matrix(temp[,5:(5+d-1)])
      
      # Create the Z1Z2 factor
      Z1Z2 <- factor(Z1 * 2 + Z2, levels = c(0, 1, 2, 3))
      
      # Check if each stratum has sufficient observations
      stratum_counts <- table(Z1Z2)
      if (all(stratum_counts > 1)) {
        Z1Z2 <- factor(Z1 * 2 + Z2)
        invisible(capture.output(multinomial_model <- nnet::multinom(Z1Z2 ~ X), file = NULL))
        propens <- fitted(multinomial_model)
        propens01 <- propens[,"1"]
        propens10 <- propens[,"2"]
        
        mu_p <- mean((Z1 == 1 & Z2 == 0) * Y / propens10) - mean((Z1 == 0 & Z2 == 1) * Y / propens01)
        propens <- glm(Z1 ~ X + Z2, family = 'binomial')$fitted.value
        mu_m <- mean(Z1 * Y / propens) - mean((1 - Z1) * Y / (1 - propens))
        if (is.finite(mu_p) && is.finite(mu_m)) {
          valid_sample <- TRUE
        } else{
          cat("non-finite estimates.\n")
        }
        
      } else{
        cat("Resampling due to missing observations in one or more strata.\n")
      }
    }
    
    
    if (mu_p >= mu_m) {
      return(c(mu_1 = mu_m, mu_2 = -mu_p))
    } else {
      return(c(mu_1 = mu_p, mu_2 = -mu_m))
    }
  }
  
  # Run bootstrap in parallel
  results <- parLapply(cl, 1:num_bootstrap, bootstrap_function)
  results <- do.call(rbind, results)
  mu_1 <- results[, "mu_1"]
  mu_2 <- results[, "mu_2"]
  
  # Stop the cluster
  stopCluster(cl)
  
  if(mu_p_s >= mu_m_s){
    sigma_p = sd(mu_1)
    sigma_m = sd(mu_2)
  }else{
    sigma_p = sd(mu_2)
    sigma_m = sd(mu_1)
  }
  
  Bn = pmax(sqrt(dim(dataset)[1]) * ( (mu_1) - mean(mu_1) ) / sigma_p, sqrt(dim(dataset)[1]) * (  (mu_2) - mean(mu_2) ) / sigma_m   )
  invBn = quantile(Bn, 1 - beta)
  lambda1 = min(c(mean(mu_1) + sigma_p *  invBn/sqrt(dim(dataset)[1] ), 0) )
  lambda2 = min(c(mean(mu_2) + sigma_m *  invBn/sqrt(dim(dataset)[1] ), 0) )
  
  
  cn = quantile(pmax(sqrt(dim(dataset)[1]) * ( mean(mu_1) - mu_1 + lambda1 ) / 
                       sigma_p, sqrt(dim(dataset)[1]) * ( mean(mu_2) - mu_2 + lambda2 ) / sigma_m)  , 1 - alpha + beta  )
  
  l1 = c()
  for (t in seq(-100,100,0.01)) {
    
    if(max(c(sqrt(dim(dataset)[1]) * (mean(mu_1) - t)/ sigma_p, sqrt(dim(dataset)[1]) * (mean(mu_2) + t)/ sigma_m ) ) < cn ){
      l1 = append(l1, t)
    }
  }
  
  u = max(l1)
  l = min(l1)
  
  
  
  
  lower = c(min(c(mu_m_s,mu_p_s)), l  )
  upper = c(max(c(mu_m_s,mu_p_s)), u  )
  result = rbind(lower, upper)
  colnames(result) = c('estimate',paste0((1-alpha)*100, '% CI'))
  rownames(result) = c('lower bound','upper bound')
  return(result)
}

get_ATE = function(dataset){
  Z1 = dataset[,'Treatment'] 
  Y = dataset[,'Y']
  Y_c = dataset[,'Y_c']
  n = dim(dataset)[1]
  Y1 = rep(0,n)
  Y0 = rep(0,n)
  for (i in 1:n) {
    if(Z1[i] == 1){
      Y1[i] = Y[i]
      Y0[i] = Y_c[i]
    }else{
      Y1[i] = Y_c[i]
      Y0[i] = Y[i]
    }
  }
  return(mean(Y1 - Y0))
}

main_estimation = function(x = 1, n_sim = 500){
  
  
  for (n in c(1000,2000,5000)) {
    for (d in c(5,10,20)) {
      for (m in c('homo','heter')) {
        for (delta in c(0,1)) {
          for (p in c(0.7,0.8,0.9)) {
            
            ate_list = rep(0,n_sim)
            cate_list = rep(0,n_sim)
            ci_ate_list = rep(0, n_sim)
            ci_cate_list = rep(0, n_sim)
            
            for (folder in seq(1,n_sim,1)) {
              if (folder %% 10 == 0) {
                print(paste("Iteration:", folder))
              }
              path = 'Differential_Effect/Simulation/'
              path_result = 'Differential_Effect/Simulation_new/'
              file_dic = paste0(path,'n_', n, '/d_', d, '/', m, 
                                '/delta_', delta,'/p_0', p*10, '/')
              file_dic_result = paste0(path_result,'n_', n, '/d_', d, '/', m, 
                                '/delta_', delta,'/p_0', p*10, '/')
              dataset = read.csv(paste0(file_dic, folder,'.csv'))
              b_ate = bound_ATE_new_bootstrap(dataset)

              ate = get_ATE(dataset)

              if(ate >= b_ate[1,1] & ate <= b_ate[2,1] ){
                ate_list[folder] = 1
              }
              if(ate >= b_ate[1,2] & ate <= b_ate[2,2]){
                ci_ate_list[folder] = 1
              }
              
              print(paste0('Current n is ', n, '. Current d is ', d, '. Current model is '
                           , m, '. Current delta is ', delta,
                           '. Current p is ', p, '. Current folder is ', folder))
              
            }
            result1 = t(rbind(ate_list,  ci_ate_list))
            colnames(result1) = c('ATE','CI ATE')
            result2 = t(apply(result1,2,mean))
            
            # Check if the directory exists, and create it if it doesn't
            if (!dir.exists(file_dic_result)) {
              dir.create(file_dic_result, recursive = TRUE)
            }
            
            write.csv(result1, file = paste0(file_dic_result,'result1.csv'),row.names = FALSE)
            write.csv(result2, file = paste0(file_dic_result,'result2.csv'),row.names = FALSE)
            
          }
        }
      }
    }
  }
}

# main_estimation()

generate_table = function(n = 1000, target = 'ATE'){
  s = 0
  for (d in c(5,10,20)) {
    for (m in c('homo','heter')) {
      
      subtable = matrix(rep(0,6*1), 1,6)
      k = 1
      
      for (delta in c(0,1)) {
        
        for (p in c(0.7,0.8,0.9)) {
          path = 'Differential_Effect/Simulation_new/'
          file_dic = paste0(path,'n_', n, '/d_', d, '/', m, 
                            '/delta_', delta,'/p_0', p*10, '/')
          dataset = read.csv(paste0(file_dic,'result2.csv'))
          result = dataset[,target]
          subtable[,k] = result
          k = k + 1
          
        }
      }
      s = s + 1
      if(s == 1){wholetable = subtable}
      if(s >= 2){wholetable = rbind(wholetable, subtable)}
    }
  }
  
  
  wholetable = format(wholetable, digits = 3)
  wholetable = cbind((c('Homogeneous','Heterogeneous','Homogeneous','Heterogeneous',
                        'Homogeneous','Heterogeneous') ), wholetable)
  colnames(wholetable) = c('Outcome model','p = 0.7', 'p = 0.8', 'p = 0.9', 'p = 0.7', 
                           'p = 0.8', 'p = 0.9')
  return(wholetable)
  
}

# print(xtable(generate_table(n=1000,target = 'CI.ATE')))  
# print(xtable(generate_table(n=2000,target = 'CI.ATE')))  
# print(xtable(generate_table(n=5000,target = 'CI.ATE')))  

# missing 1000, 5, heter, 1, 0.9

main_estimation_subset = function(x = 1, n_sim = 500){
  
  
  for (n in c(1000,2000,5000)) {
    for (d in c(5,10,20)) {
      for (m in c('homo','heter')) {
        for (delta in c(0,1)) {
          for (p in c(0.7,0.8,0.9)) {
            
            if (n == 1000) {
              next
            } else if (n == 2000 && d == 5 && m == 'homo') {
              next
            } else if (n == 2000 && d == 5 && m == 'heter' && delta == 0) {
              next
            } else if (n == 2000 && d == 5 && m == 'heter' && delta == 1 && p %in% c(0.7)) {
              next
            }else {
              
            }
            
            ate_list = rep(0,n_sim)
            cate_list = rep(0,n_sim)
            ci_ate_list = rep(0, n_sim)
            ci_cate_list = rep(0, n_sim)
            
            relative_length_ate_list = rep(0, n_sim)
            
            for (folder in seq(1,n_sim,1)) {
              if (folder %% 10 == 0) {
                print(paste("Iteration:", folder))
              }
              path = 'C:/Users/Kan Chen/Dropbox/kan_bingkai_jeff_dylan/Differential_Effect/Simulation/'
              path_result = 'C:/Users/Kan Chen/Dropbox/kan_bingkai_jeff_dylan/Differential_Effect/Simulation_new/'
              file_dic = paste0(path,'n_', n, '/d_', d, '/', m, 
                                '/delta_', delta,'/p_0', p*10, '/')
              file_dic_result = paste0(path_result,'n_', n, '/d_', d, '/', m, 
                                       '/delta_', delta,'/p_0', p*10, '/')
              dataset = read.csv(paste0(file_dic, folder,'.csv'))
              b_ate = bound_ATE_new_bootstrap(dataset)
              
              ate = get_ATE(dataset)
              
              relative_length_ate_list[folder] = (b_ate[2,2] - b_ate[1,2])/(b_ate[2,1] - b_ate[1,1])
              print((b_ate[2,2] - b_ate[1,2])/(b_ate[2,1] - b_ate[1,1]))
              
              # Check for finite values (neither NA nor NaN) before evaluating the conditions
              if (is.finite(ate) && is.finite(b_ate[1,1]) && is.finite(b_ate[2,1])) {
                if (ate >= b_ate[1,1] & ate <= b_ate[2,1]) {
                  ate_list[folder] = 1
                }
              }
              
              if (is.finite(ate) && is.finite(b_ate[1,2]) && is.finite(b_ate[2,2])) {
                if (ate >= b_ate[1,2] & ate <= b_ate[2,2]) {
                  ci_ate_list[folder] = 1
                }
              }
              
              print(paste0('Current n is ', n, '. Current d is ', d, '. Current model is '
                           , m, '. Current delta is ', delta,
                           '. Current p is ', p, '. Current folder is ', folder))
              
            }
            result1 = t(rbind(ate_list,  ci_ate_list))
            colnames(result1) = c('ATE','CI ATE')
            result2 = t(apply(result1,2,mean))
            
            result3 = t(relative_length_ate_list)
            colnames(result3) = c('relative_length_ate')
            result3 = t(apply(result3,2,mean))
            
            # Check if the directory exists, and create it if it doesn't
            if (!dir.exists(file_dic_result)) {
              dir.create(file_dic_result, recursive = TRUE)
            }
            
            write.csv(result1, file = paste0(file_dic_result,'result1.csv'),row.names = FALSE)
            write.csv(result2, file = paste0(file_dic_result,'result2.csv'),row.names = FALSE)
            write.csv(result3, file = paste0(file_dic_result,'result3.csv'),row.names = FALSE)
          }
        }
      }
    }
  }
}

 main_estimation_subset()
