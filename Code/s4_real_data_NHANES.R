library(foreign)
library(dplyr)
library(MASS)

file_path <- "kan_bingkai_dylan/Differential_Effect/Data/NHANES/"

DEMO_J <- read.xport(paste0(file_path, 'DEMO_J.XPT'))
DUQ_J <- read.xport(paste0(file_path, 'DUQ_J.XPT'))
PBCD_J <- read.xport(paste0(file_path, 'PBCD_J.XPT'))
SMQ_J <- read.xport(paste0(file_path, 'SMQ_J.XPT'))
PAQ_J <- read.xport(paste0(file_path, 'PAQ_J.XPT'))
ALQ_J <- read.xport(paste0(file_path, 'ALQ_J.XPT'))
DR1TOT_J <- read.xport(paste0(file_path, 'DR1TOT_J.XPT'))


mdat1 <- DEMO_J %>%
  mutate(age = RIDAGEYR, gender = RIAGENDR, race = RIDRETH3, 
         education = DMDEDUC2, income_ratio = INDFMPIR) %>%
  dplyr::select(SEQN, age, gender, race, education, income_ratio)

# change definition of smoker to match Rosenbaum Chance
mdat2 <- mdat1 %>%
  left_join(SMQ_J, by = 'SEQN') %>%
  mutate(smoker = case_when(
    SMD650 >= 10 ~ TRUE,      # Define smoker
    SMQ020 == 2 ~ FALSE,      # Define nonsmoker
    TRUE ~ NA                  # Assign NA to others
  )) %>%
  filter(!is.na(smoker)) %>%    # Remove rows that don't fit either category
  dplyr::select(SEQN, smoker, age, gender, race, education, income_ratio)


mdat3 <- mdat2 %>%
  left_join(PBCD_J, by = 'SEQN') %>%
  mutate(lead = LBDBPBSI, cadmium = LBDBCDSI) %>%
  dplyr::select(SEQN,lead,smoker,cadmium,age, gender, race, education, income_ratio)

mdat4 <- mdat3 %>%
  left_join(DUQ_J, by = 'SEQN') %>%
  filter(DUQ240 == 1 | DUQ240 == 2) %>%
  mutate(drug_use = ifelse(DUQ240 == 1, TRUE, FALSE)) %>%
  dplyr::select(cadmium,smoker,drug_use,lead,age, gender, 
                race, education, income_ratio,SEQN)

#mdat5 <- mdat4 %>%
#  left_join(PAQ_J, by = 'SEQN') %>%
#  filter(PAQ650 == 1 | PAQ650 == 2) %>%
#  mutate(physical_activity = ifelse(PAQ650 == 1, FALSE, TRUE)) %>%
#  dplyr::select(cadmium,smoker,drug_use,lead,age, gender, 
#                race, education,physical_activity, income_ratio,SEQN)

# used to be ALQ130 >= 130
mdat6 <- mdat4 %>%
  left_join(ALQ_J, by = 'SEQN') %>%
  filter(ALQ130  != 777 & ALQ130  != 999) %>%
  mutate(alcohol_abuse = ifelse(ALQ130 >= 5, TRUE, FALSE)) %>%
  dplyr::select(cadmium,smoker,drug_use,lead,age, gender, 
                race, education, alcohol_abuse,
                income_ratio,SEQN)
# used to be DR1TTFAT >= 150
mdat7 <- mdat6 %>%
  left_join(DR1TOT_J, by = 'SEQN') %>%
  mutate(fat_intake = ifelse(DR1TTFAT >= 150, TRUE, FALSE)) %>%
  dplyr::select(cadmium,smoker,drug_use,lead,age, gender, 
                race, education, alcohol_abuse,fat_intake,
                income_ratio,SEQN)



NHANES <- mdat4[complete.cases(mdat4),]


write.csv(NHANES, paste0(file_path, 'NHANES_new.csv'),row.names = FALSE)


NHANES_2PL <- mdat7[complete.cases(mdat7),]


write.csv(NHANES_2PL, paste0(file_path, 'NHANES_2PL_new.csv'),row.names = FALSE)



#### Implementation of real data ####
source('kan_bingkai_dylan/Differential_Effect/Code/Estimation_two_stage_DE.R')
source('kan_bingkai_dylan/Differential_Effect/Code/bound_ATE_new.R')
################ Lead ##########################
data_example_lead_age <- NHANES
colnames(data_example_lead_age)[1:3] <- c('Y',"Treatment",'Z')

# Convert categorical variable  to one-hot encoding
data_example_lead_age$race <- as.factor(data_example_lead_age$race)
race_onehot <- model.matrix(~ race - 1, data = data_example_lead_age)
race_onehot_numeric <- as.data.frame(data.matrix(race_onehot))

data_example <- bind_cols(data_example_lead_age %>% dplyr::select(Y, Treatment, Z, lead,
                age, gender, education, income_ratio), race_onehot_numeric)
data_example$SEQN <- data_example_lead_age$SEQN
#### age ####
bound_ATE_new_bootstrap(data_example)

#### gender ####
data_example_lead_gender <- data_example[c("Y","Treatment","Z","lead",
                                           "gender","income_ratio","race1",
                                           "race2","race3","race4","race6",
                                           "race7",
                                           "education","age","SEQN")]

# Female
bound_CATE(x = 2,data_example_lead_gender)
# Male
bound_CATE(x = 1,data_example_lead_gender)

#### education ####
data_example_lead_edu <- data_example[c("Y","Treatment","Z","lead",
                                                 "education","income_ratio",
                                                 "race1","race2","race3",
                                                 "race4","race6","race7",
                                                 "gender","age","SEQN")]
# High school graduate/GED or equivalent
bound_CATE(x = 3,data_example_lead_edu)
# College graduate or above
bound_CATE(x = 4,data_example_lead_edu)


library(gtsummary)

table1 <- NHANES %>% tbl_summary(by = smoker) %>% 
  modify_caption("**Table 1.**")

# Convert the gtsummary object to kable in LaTeX format
kable_latex <- as_kable(table1, format = "latex")

# Print the raw LaTeX code to the console
cat(kable_latex)



