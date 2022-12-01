# This is the analysis of the simulation results
library(xtable)
library(ggplot2)
source('Estimation_DE.R')

generate_table = function(n = 1000, target = 'ATE'){
    s = 0
    for (d in c(5,10,20)) {
      for (m in c('homo','heter')) {
        
        subtable = matrix(rep(0,6*1), 1,6)
        k = 1
        
        for (delta in c(0,1)) {

          for (p in c(0.7,0.8,0.9)) {
              path = 'C:/Users/chenk/Dropbox/kan_bingkai_dylan/Differential_Effect/Simulation/'
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

italic <- function(x){
  paste0('{\\emph{ ', x, '}}')
}

large <- function(x){
  paste0('{\\Large{\\bfseries ', x, '}}')
}



generate_plot = function(n = 5000,d = 20, delta = 1, p = 0.9, folder = 1){
  m1 = 'homo'
  m2 = 'heter'
  path = 'C:/Users/chenk/Dropbox/kan_bingkai_dylan/Differential_Effect/Simulation/'
  file_dic1 = paste0(path,'n_', n, '/d_', d, '/', m1, 
                    '/delta_', delta,'/p_0', p*10, '/')
  file_dic2 = paste0(path,'n_', n, '/d_', d, '/', m2, 
                     '/delta_', delta,'/p_0', p*10, '/')
  dataset1 = read.csv(paste0(file_dic1, folder,'.csv'))
  dataset2 = read.csv(paste0(file_dic2, folder,'.csv'))
  
  xlist = seq(0,0.5,0.01)
  
  catelist1 = rep(0,length(xlist))
  lower.catelist1 = rep(0,length(xlist))
  upper.catelist1 = rep(0,length(xlist))
  
  catelist2 = rep(0,length(xlist))
  lower.catelist2 = rep(0,length(xlist))
  upper.catelist2 = rep(0,length(xlist))
  
  i = 1
  for (x in xlist) {
    catelist1[i] = get_CATE(x,dataset1)
    lower.catelist1[i] = (bound_CATE(x,dataset1))[1,1]
    upper.catelist1[i] = (bound_CATE(x,dataset1))[2,1]
    
    catelist2[i] = get_CATE(x,dataset2)
    lower.catelist2[i] = (bound_CATE(x,dataset2))[1,1]
    upper.catelist2[i] = (bound_CATE(x,dataset2))[2,1]
    
    i = i + 1
    print(i)
  }
  
  df = as.data.frame(cbind(xlist, catelist1, upper.catelist1, lower.catelist1,
             catelist2, upper.catelist2, lower.catelist2))
  
  par(mfrow=c(1,2),xpd=TRUE)
  
  
    
  plot(xlist, catelist1, ylim = c(-4,7),xlab = 'x', ylab = 'Treatment Effect',type = "l", lty=1, lwd=2,
       main = 'Homogeneous Effect',cex.lab = 1, cex.axis = 1)
  lines(xlist,upper.catelist1,type = "l", lty=2, lwd=2)
  lines(xlist,lower.catelist1,type = "l", lty=3, lwd=2)
  
  legend("topright", inset=c(0.12, 0.01), legend=c("True", "Upper","Lower"),
         col=c("black", "black","black"), lty=c(1,2,3), cex=0.7,box.col = "white")
  
  plot(xlist, catelist2, ylim = c(-4,7),xlab = 'x', ylab = 'Treatment Effect',type = "l", lty=1, lwd=2,
       main = 'Heterogeneous Effect',cex.lab = 1, cex.axis = 1)
  lines(xlist,upper.catelist2,type = "l", lty=2, lwd=2)
  lines(xlist,lower.catelist2,type = "l", lty=3, lwd=2)
  legend("topright", inset=c(0.12, 0.01), legend=c("True", "Upper","Lower"),
         col=c("black", "black","black"), lty=c(1,2,3), cex=0.7,box.col = "white")
  
}



#print(xtable(generate_table(n=5000,target = 'CATE')))  










