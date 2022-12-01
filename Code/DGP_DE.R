# This is DGP
library('MASS')


treatment_received = function(Z, co_status){
  # co_status: complier status; 0: complier 1: always takers 2: never takers
  Treatment = Z
  Treatment[which(co_status == 1)] = 1
  Treatment[which(co_status == 2)] = 0
  return(Treatment)
}



DPG = function(p_at = 0.1, p_nt = 0.1 , p_co = 0.8, d = 10, n = 1000, theta = 2.0, delta = 1.0, model = 'homo'){
  # p_at: proportion of alway takers
  # p_nt: proportion of never takers
  # p_co: proportion of compliers (should be 1 - p_at - p_nt when monotonicity assumption holds)
  # d: dimension of observed covariates X
  # n: number of observations
  # theta: coefficient of T
  X = matrix(rep(0,n*d), n, d) 
  Z = rbinom(n,1,1/2)
  co_status = sample(c(0,1,2), n, replace = TRUE, prob = c(p_co, p_at, p_nt))
  X[which(co_status == 0),] = mvrnorm(length(which(co_status == 0)),rep(0,d),diag(d))
  X[which(co_status == 1),] = mvrnorm(length(which(co_status == 1)),rep(0.25,d),diag(d))
  X[which(co_status == 2),] = mvrnorm(length(which(co_status == 2)),rep(0.5,d),diag(d))
  Treatment = treatment_received(Z, co_status) # Generate treatment received
  e1 = rnorm(n, mean = 0, sd = 1)
  U = 1 + Treatment + delta * Z + e1
  U_c = 1 + (1 - Treatment) + delta * Z + e1
  if(model == 'homo'){
    e2 = rnorm(n, mean = 0, sd = 1)
    Y = Treatment * theta + apply(X[,1:d], 1, sum) + U + e2
    Y_c = (1 - Treatment) * theta + apply(X[,1:d], 1, sum) + U_c + e2
  }else{
    e2 = rnorm(n, mean = 0, sd = 1)
    Y = Treatment * (2 * sin(X[,1]/2) + 1.5) + apply(X[,1:d], 1, sum) + U + e2
    Y_c = (1 - Treatment) * (2 * sin(X[,1]/2) + 1.5) + apply(X[,1:d], 1, sum) + U_c + e2
  }
  
  df = cbind(Y,Treatment,Z,co_status,X, Y_c)
  return(df)
}



main = function(n_sim = 500){
  
  for (n in c(1000,2000,5000)) {
    for (d in c(5,10,20)) {
      for (m in c('homo','heter')) {
        for (delta in c(0,1)) {
          for (p in c(0.7,0.8,0.9)) {
            for (i in seq(1,n_sim,1)) {
              path = 'C:/Users/chenk/Dropbox/kan_bingkai_dylan/Differential_Effect/Simulation/'
              file_dic = paste0(path,'n_', n, '/d_', d, '/', m, '/delta_', delta,'/p_0', p*10, '/')
              dataset = DPG(p_at = (1 - p)/2, p_nt = (1 - p)/2 , p_co = p, d = d, 
                            n = n, theta = 2.0, delta = delta, model = m )
              
              write.csv(dataset, file = paste0(file_dic, i,'.csv'),row.names = FALSE)
              print(paste0('Current n is ', n, '. Current d is ', d, '. Current model is '
                           , m, '. Current delta is ', delta,
                           '. Current p is ', p, '. Current folder is ', i))
        }

            
          }

        }
        
      }
      
    }
  }
  
}


# Generate Data
#main()










