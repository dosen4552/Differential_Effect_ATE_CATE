library(dplyr)


file_path <- "C:/Users/chenk/Dropbox/kan_bingkai_dylan/Differential_Effect/Data/YRBSS/XXHq.xlsx"


YRBSS2019 <- readxl::read_xlsx(file_path)


####### related variables
# q75: soft drink consumption.

# q30: Have you ever tried cigarette smoking
# q31: How old were you when you first tried cigarette smoking, even one or two puffs?
# q32: During the past 30 days, on how many days did you smoke cigarettes?
# q33: During the past 30 days, on the days you smoked, how many cigarettes did you smoke per day?

# q17: During the past 12 months, how many times were you in a physical fight?
# q18: During the past 12 months, how many times were you in a physical fight on school property?
# q23: During the past 12 months, have you ever been bullied on school property?


########## demographics
# BMIPCT: The BMI percentile for age and sex
# raceeth: race and ethnity
# weight: A weight applied to each record to adjust for school and student nonresponse and oversampling of black and Hispanic students.


mdat1 <- YRBSS2019 %>%
  mutate(age = ifelse(q1==1,12,ifelse(q1==2,13,ifelse(q1==3,14,ifelse(q1==4,15,
          ifelse(q1==5,16,ifelse(q1==6,17,ifelse(q1==7,18,NA))))))), 
         gender = q2, race = raceeth, BMI = BMIPCT, grade = q3, soft_drink = ifelse(q75==1, FALSE, TRUE),
         smoking = ifelse(q30==1, TRUE, FALSE), violence1 = ifelse(q17==1, FALSE, TRUE), 
         violence2 = ifelse(q18==1, FALSE, TRUE)) %>%
  select(violence1,soft_drink, smoking, violence2, age, 
         gender, race, BMI, grade, record)

YRBSS <- mdat1[complete.cases(mdat1),]  
  
write.csv(YRBSS, paste0(file_path, 'YRBSS.csv'),row.names = FALSE)  

#### Implementation of real data ####
source('Estimation_DE.R')  
  
data_example_violence_age <- YRBSS
colnames(data_example_violence_age)[1:3] <- c('Y',"Treatment",'Z')

bound_ATE(data_example_violence_age)


#### Age ########
bound_CATE(x = 13,data_example_violence_age)



#### BMI percentile for age and sex ####
data_example_violence_BMI <- data_example_violence_age[c("Y","Treatment","Z","violence2",
                                                    "BMI" ,"age", "gender","race",        
                                                    "grade","age","record")]


bound_CATE(x = 49,data_example_violence_BMI)


#### gender ####


data_example_violence_gender <- data_example_violence_age[c("Y","Treatment","Z","violence2",
                                                            "gender","age", "BMI","race",        
                                                            "grade","age","record")]

bound_CATE(x = 1,data_example_violence_gender)

bound_CATE(x = 2,data_example_violence_gender)



#### race ####

data_example_violence_race <- data_example_violence_age[c("Y","Treatment","Z","violence2",
                                                          "race","age", "BMI","gender",        
                                                            "grade","age","record")]


bound_CATE(x = 3,data_example_violence_race)

bound_CATE(x = 5,data_example_violence_race)


#### grade ####

data_example_violence_grade <- data_example_violence_age[c("Y","Treatment","Z","violence2",
                                                          "grade","age", "BMI","gender",        
                                                          "race","age","record")]



bound_CATE(x = 3,data_example_violence_grade)


bound_CATE(x = 4,data_example_violence_grade)


#### Graph for age and BMI percentile for age and sex ####
xlist_age = seq(13,18,0.1)
xlist_BMI = seq(0.1,99,0.5)

lower.catelist1 = rep(0,length(xlist_age))
upper.catelist1 = rep(0,length(xlist_age))

lower.catelist2 = rep(0,length(xlist_BMI))
upper.catelist2 = rep(0,length(xlist_BMI))



i = 1
for (x in xlist_age) {
  lower.catelist1[i] = (bound_CATE(x,data_example_violence_age))[1,1]
  upper.catelist1[i] = (bound_CATE(x,data_example_violence_age))[2,1]
  i = i + 1
  print(i)
}


i = 1
for (x in xlist_BMI) {
  lower.catelist2[i] = (bound_CATE(x,data_example_violence_BMI))[1,1]
  upper.catelist2[i] = (bound_CATE(x,data_example_violence_BMI))[2,1]
  i = i + 1
  print(i)
}

par(mfrow=c(1,2),xpd=FALSE)

plot(xlist_age, upper.catelist1,xlab = 'Age', ylab = 'Treatment Effect',type = "l", lty=2, lwd=2,
     main = 'Age v.s. Treatment Effect',cex.lab = 1, cex.axis = 1, ylim = c(-0.3,0.21), xlim = c(13,18))
lines(xlist_age,lower.catelist1,type = "l", lty=3, lwd=2, col = 'grey')

legend("topright", inset=c(0.12, 0.01), legend=c("Upper","Lower"),
       col=c("black",'grey'), lty=c(2,3), cex=0.7,box.col = "white")

plot(xlist_BMI, upper.catelist2, xlab = 'BMI Percentile for Age and Sex', ylab = 'Treatment Effect',type = "l", lty=2, lwd=2,
     main = 'BMI Percentile v.s. Treatment Effect',cex.lab = 1, cex.axis = 1,ylim = c(-0.3,0.21), xlim = c(39,53))
lines(xlist_BMI,lower.catelist2,type = "l", lty=3, lwd=2, col = 'grey')
legend("topright", inset=c(0.12, 0.01), legend=c("Upper","Lower"),
       col=c("black",'grey'), lty=c(2,3), cex=0.7,box.col = "white")



library(gtsummary)

YRBSS %>% tbl_summary(by = soft_drink) %>% 
  modify_caption("**Table 1.**")


  
  
  

