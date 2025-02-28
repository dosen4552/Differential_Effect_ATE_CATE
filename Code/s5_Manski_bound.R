NHANES <- read.csv("/Users/kanchen/Dropbox/kan_bingkai_jeff_dylan/Differential_Effect/Data/NHANES/NHANES_new.csv")
data_example_lead_age <- NHANES
colnames(data_example_lead_age)[1:3] <- c('Y',"Treatment",'Z')


p1 <- mean(data_example_lead_age[data_example_lead_age$Z == TRUE,]$Treatment )
p0 <- mean(data_example_lead_age[data_example_lead_age$Z == FALSE,]$Treatment )
m1 <- mean(data_example_lead_age[data_example_lead_age$Z == TRUE,]$Y )
m0 <- mean(data_example_lead_age[data_example_lead_age$Z == FALSE,]$Y )

Y_max <- max(data_example_lead_age$Y)
Y_min <- 0


lower <- (m1-m0)/(p1-p0) - max(c((1-p1)/(p1-p0) * (Y_max - Y_min) , p0/(p1-p0) * (Y_max - Y_min)  )  )

upper <- (m1-m0)/(p1-p0) + max(c((1-p1)/(p1-p0) * (Y_max - Y_min), p0/(p1-p0) * (Y_max - Y_min)  )  )
