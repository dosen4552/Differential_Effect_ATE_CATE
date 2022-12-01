library(foreign)
library(dplyr)


file_path <- 'C:/Users/chenk/Dropbox/kan_bingkai_dylan/Differential_Effect/Data/NHANES/'

DEMO_J <- read.xport(paste0(file_path, 'DEMO_J.XPT'))
DUQ_J <- read.xport(paste0(file_path, 'DUQ_J.XPT'))
PBCD_J <- read.xport(paste0(file_path, 'PBCD_J.XPT'))
SMQ_J <- read.xport(paste0(file_path, 'SMQ_J.XPT'))


mdat1 <- DEMO_J %>%
  mutate(age = RIDAGEYR, gender = RIAGENDR, race = RIDRETH3, 
         education = DMDEDUC2, income_ratio = INDFMPIR) %>%
  select(SEQN, age, gender, race, education, income_ratio)


mdat2 <- mdat1 %>%
  left_join(SMQ_J, by = 'SEQN') %>%
  mutate(smoker = ifelse(SMD650 >= 10, TRUE, ifelse(SMD650 == 1, 
                                                  FALSE, NA) ) ) %>%
  filter(!is.na(smoker)) %>%
  select(SEQN,smoker, age, gender, race, education, income_ratio)


mdat3 <- mdat2 %>%
  left_join(PBCD_J, by = 'SEQN') %>%
  mutate(lead = LBDBPBSI, cadmium = LBDBCDSI) %>%
  select(SEQN,lead,smoker,cadmium,age, gender, race, education, income_ratio)

mdat4 <- mdat3 %>%
  left_join(DUQ_J, by = 'SEQN') %>%
  filter(DUQ240 == 1 | DUQ240 == 2) %>%
  mutate(drug_use = ifelse(DUQ240 == 1, TRUE, FALSE)) %>%
  select(lead,smoker,drug_use,cadmium,age, gender, 
         race, education, income_ratio,SEQN)


NHANES <- mdat4[complete.cases(mdat4),]


write.csv(NHANES, paste0(file_path, 'NHANES.csv'),row.names = FALSE)



#### Implementation of real data ####
source('Estimation_DE.R')
################ Lead ##########################
data_example_lead_age <- NHANES
colnames(data_example_lead_age)[1:3] <- c('Y',"Treatment",'Z')


#### age ####
bound_ATE(data_example_lead_age)

bound_CATE(x = 50,data_example_lead_age)

#### income ratio ####
data_example_lead_income <- data_example_lead_age[c("Y","Treatment","Z","cadmium",
                                "income_ratio" , "gender","race",        
                                "education","age","SEQN")]

bound_CATE(x = 4,data_example_lead_income)

#### gender ####
data_example_lead_gender <- data_example_lead_age[c("Y","Treatment","Z","cadmium",
                                                     "gender","income_ratio","race",        
                                                    "education","age","SEQN")]


bound_CATE(x = 1,data_example_lead_gender)

#### race ####
data_example_lead_race <- data_example_lead_age[c("Y","Treatment","Z","cadmium",
                                                    "race","income_ratio","gender",        
                                                    "education","age","SEQN")]

bound_CATE(x = 3,data_example_lead_race)


#### education ####
data_example_lead_edu <- data_example_lead_age[c("Y","Treatment","Z","cadmium",
                                                  "education","income_ratio","gender",        
                                                  "race","age","SEQN")]
bound_CATE(x = 3,data_example_lead_edu)

#### Graph for age and income ratio ####
xlist_age = seq(20,69,1)
xlist_income = seq(0.1,5,0.1)

lower.catelist1 = rep(0,length(xlist_age))
upper.catelist1 = rep(0,length(xlist_age))

lower.catelist2 = rep(0,length(xlist_income))
upper.catelist2 = rep(0,length(xlist_income))



i = 1
for (x in xlist_age) {
  lower.catelist1[i] = (bound_CATE(x,data_example_lead_age))[1,1]
  upper.catelist1[i] = (bound_CATE(x,data_example_lead_age))[2,1]
  i = i + 1
  print(i)
}


i = 1
for (x in xlist_income) {
  lower.catelist2[i] = (bound_CATE(x,data_example_lead_income))[1,1]
  upper.catelist2[i] = (bound_CATE(x,data_example_lead_income))[2,1]
  i = i + 1
  print(i)
}

par(mfrow=c(1,2),xpd=FALSE)

plot(xlist_age, upper.catelist1,xlab = 'Age', ylab = 'Treatment Effect',type = "l", lty=2, lwd=2,
     main = 'Age v.s. Treatment Effect',cex.lab = 1, cex.axis = 1, ylim = c(-0.15,0.15), xlim = c(30,65))
lines(xlist_age,lower.catelist1,type = "l", lty=3, lwd=2, col = 'grey')

legend("topright", inset=c(0.12, 0.01), legend=c("Upper","Lower"),
       col=c("black",'grey'), lty=c(2,3), cex=0.7,box.col = "white")

plot(xlist_income, upper.catelist2, xlab = 'Ratio of Family Income to Poverty', ylab = 'Treatment Effect',type = "l", lty=2, lwd=2,
     main = 'Income Ratio v.s. Treatment Effect',cex.lab = 1, cex.axis = 1,ylim = c(-0.15,0.15), xlim = c(0.1,2.9))
lines(xlist_income,lower.catelist2,type = "l", lty=3, lwd=2, col = 'grey')
legend("topright", inset=c(0.12, 0.01), legend=c("Upper","Lower"),
       col=c("black",'grey'), lty=c(2,3), cex=0.7,box.col = "white")



library(gtsummary)

NHANES %>% tbl_summary(by = smoker) %>% 
  modify_caption("**Table 1.**")



