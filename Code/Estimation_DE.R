library('MASS')


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

get_CATE = function(x,dataset){
  Z1 = dataset[,'Treatment'] 
  Y = dataset[,'Y']
  Y_c = dataset[,'Y_c']
  n = dim(dataset)[1]
  d = dim(dataset)[2] - 5
  X = dataset[,5:(5+d-1)]
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
  ite = Y1 - Y0
  cate = kernel_regression(x,  X[,1], ite, b =bw.ucv(X[,1]))
  return(cate)
}

kde_f_hat = function(x, X){
  h = bw.ucv(X)
  K = gausinKernel(x-X,h)
  de = sum(K)/(length(X)[1]*h)
  return(de)
}


bound_ATE = function(dataset, alpha = 0.05){
  Z1 = data.matrix(dataset[,'Treatment']) 
  Z2 = data.matrix(dataset[,'Z']) 
  Y = data.matrix(dataset[,'Y'])
  d = dim(dataset)[2] - 5
  X = data.matrix(dataset[,5:(5+d-1)]) 
  mu_p = mean(Y[which(Z1 == 1 & Z2 == 0)]) - mean(Y[which(Z1 == 0 & Z2 == 1)]) 
  propens = glm(Z1 ~ X + Z2, family = 'binomial')$fitted.value
  mu_m =  mean(Z1 * Y/propens) - mean((1-Z1) * Y / (1 - propens)  )
  sigma_p = sqrt(sum((Y[which(Z1 == 1 & Z2 == 0)] - mean(Y[which(Z1 == 1 & Z2 == 0)]))^2)/
    length(which(Z1 == 1 & Z2 == 0))^2 + sum((Y[which(Z1 == 0 & Z2 == 1)] - 
      mean(Y[which(Z1 == 0 & Z2 == 1)]))^2)/length(which(Z1 == 0 & Z2 == 1))^2)
  sigma_m = sqrt(sum( (Z1 * Y/propens - (1-Z1) * Y / (1 - propens) - mean(Z1 * Y/propens - 
            (1-Z1) * Y / (1 - propens)))^2 )/(dim(dataset)[1]^2))
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


#function to calculate Gaussian kernel
gausinKernel = function(x,b){
  K = (1/((sqrt(2*pi))))*exp(-0.5 *(x/b)^2)
  return(K)
}


kernel_regression = function(x,X,Y,b = 0.5){
  K = gausinKernel(x-X,b)
  Ksum = sum(K)
  weight = K/Ksum
  y =  sum(weight*Y)
  return(y)
}

phi = function(y,z1,p){
  return(z1*y/p - (1 - z1) * y / (1 - p) )
}


bound_CATE = function(x, dataset, alpha = 0.05){
  Z1 = data.matrix(dataset[,'Treatment']) 
  Z2 = data.matrix(dataset[,'Z']) 
  Y = data.matrix(dataset[,'Y'])
  d = dim(dataset)[2] - 5
  X = data.matrix(dataset[,5:(5+d-1)]) 
  h1 = bw.ucv(X[which(Z1 == 1 & Z2 == 0),1])
  h2 = bw.ucv(X[which(Z1 == 0 & Z2 == 1),1])
  h3 = bw.ucv(X[,1])
  mu_p_x = kernel_regression(x, X = X[which(Z1 == 1 & Z2 == 0),1], 
   Y[which(Z1 == 1 & Z2 == 0)], b = h1) - 
  kernel_regression(x, X = X[which(Z1 == 0 & Z2 == 1),1], 
   Y[which(Z1 == 0 & Z2 == 1)], b = h2)
  
  sigma_10_x = sqrt(mean((Y[which(Z1 == 1 & Z2 == 0)] - kernel_regression(x, X = X[which(Z1 == 1 & Z2 == 0),1], 
                      Y[which(Z1 == 1 & Z2 == 0)], b = h1))^2 * h1^(-1) * 
  gausinKernel(x-X[which(Z1 == 1 & Z2 == 0),1],h1) / kde_f_hat(x,X[which(Z1 == 1 & Z2 == 0),1])))
  
  sigma_01_x = sqrt(mean((Y[which(Z1 == 0 & Z2 == 1)] - kernel_regression(x, X = X[which(Z1 == 0 & Z2 == 1),1], 
                      Y[which(Z1 == 0 & Z2 == 1)], b = h2))^2 * h2^(-1) * 
  gausinKernel(x-X[which(Z1 == 0 & Z2 == 1),1],h2) / kde_f_hat(x,X[which(Z1 == 0 & Z2 == 1),1])))
  
  V_10_x = sigma_10_x^2 * 1/(2 * sqrt(pi)) / kde_f_hat(x,X[which(Z1 == 1 & Z2 == 0),1])
  V_01_x = sigma_01_x^2 * 1/(2 * sqrt(pi)) / kde_f_hat(x,X[which(Z1 == 0 & Z2 == 1),1])
  V_1_x = V_10_x + V_01_x
  c_x = V_10_x/V_1_x
  p10 = length(which(Z1 == 0 & Z2 == 1))/dim(X)[1]
  p01 = length(which(Z1 == 1 & Z2 == 0))/dim(X)[1]
  sigma_mu_p_x = sqrt(V_1_x * (p10*c_x/h1 + p01*(1-c_x)/h2 )) 
  
  propens = glm(Z1 ~ X + Z2, family = 'binomial')$fitted.value
  mu_m_x = kernel_regression(x, X = X[,1], 
      Z1 * Y/propens - (1-Z1) * Y / (1 - propens), b = h3)
  
  
  sigma_x = sqrt(sum( (phi(Y,Z1,propens) - mu_m_x)^2 )/length(Y)^2)
  
  sigma_mu_m_x = sqrt(sigma_x^2 * 1/(2 * sqrt(pi)) / kde_f_hat(x,X[,1])) / sqrt(h3)
  
  if(mu_p_x > mu_m_x){
    tau_u_x = mu_p_x
    sigma_u_x = sigma_mu_p_x
    tau_l_x = mu_m_x
    sigma_l_x = sigma_mu_m_x
  }else{
    tau_l_x = mu_p_x
    sigma_l_x = sigma_mu_p_x
    tau_u_x = mu_m_x
    sigma_u_x = sigma_mu_m_x
  }
  lower = c(tau_l_x, tau_l_x - qnorm(1-alpha/2)*sigma_l_x, tau_l_x + qnorm(1-alpha/2)*sigma_l_x  )
  upper = c(tau_u_x, tau_u_x - qnorm(1-alpha/2)*sigma_u_x, tau_u_x + qnorm(1-alpha/2)*sigma_u_x  )
  result = rbind(lower, upper)
  colnames(result) = c('estimate',paste0((1-alpha/2)*100, '% lower'),paste0((1-alpha/2)*100, '% upper'))
  rownames(result) = c('lower bound','upper bound')
  return(result)
  
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
              path = 'C:/Users/chenk/Dropbox/kan_bingkai_dylan/Differential_Effect/Simulation/'
              file_dic = paste0(path,'n_', n, '/d_', d, '/', m, 
                                '/delta_', delta,'/p_0', p*10, '/')
              dataset = read.csv(paste0(file_dic, folder,'.csv'))
              b_ate = bound_ATE(dataset)
              b_cate = bound_CATE(x,dataset)
              ate = get_ATE(dataset)
              cate = get_CATE(x,dataset)
              if(cate >= b_cate[1,1] & cate <= b_cate[2,1] ){
                cate_list[folder] = 1
              }
              if(cate >= b_cate[1,2] & cate <= b_cate[2,3]){
                ci_cate_list[folder] = 1
              }
              if(ate >= b_ate[1,1] & ate <= b_ate[2,1] ){
                ate_list[folder] = 1
              }
              if(ate >= b_ate[1,2] & ate <= b_ate[2,3]){
                ci_ate_list[folder] = 1
              }
              
              print(paste0('Current n is ', n, '. Current d is ', d, '. Current model is '
                           , m, '. Current delta is ', delta,
                           '. Current p is ', p, '. Current folder is ', folder))
              
            }
            result1 = t(rbind(ate_list, cate_list, ci_ate_list, ci_cate_list))
            colnames(result1) = c('ATE','CATE','CI ATE','CI CATE')
            result2 = t(apply(result1,2,mean))
            write.csv(result1, file = paste0(file_dic,'result1.csv'),row.names = FALSE)
            write.csv(result2, file = paste0(file_dic,'result2.csv'),row.names = FALSE)

          }
        }
      }
    }
  }
}


#main_estimation()







