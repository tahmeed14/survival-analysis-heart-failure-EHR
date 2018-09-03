setwd("/home/tureen/heart_ICD")
source("/exports/EHR/data/read_diag_file.R")
library(dplyr)
library(devtools)
library(readr)
library(PheWAS)
library(tidyr)
library(ggplot2)

# Step 1
# Reading raw tables from EHR team data source
## Data Manipulation(s) Done ## 
# Age, DaysSinceBirth were changed to numeric values

outcome.icd9 <- "428"
outcome.icd10 <- "I50"

lab.data <- read.csv("/exports/EHR/BDSI_labs.csv", header = T)
lab.data$DaysSinceBirth <- as.numeric(lab.data$DaysSinceBirth)
lab.data$StudyID <- as.character(lab.data$StudyID)

#age.data <- read.csv("/exports/EHR/data/Freeze2/BDSI_MGI_Freeze2_Age.csv",header=F, skip=1, quote="",comment="",sep=",", fill=TRUE)
#age.data <- read.csv("/exports/EHR/data/Freeze2/BDSI_MGI_Freeze2_Age_Sex.csv", header = T)
age.data <- read.csv("/exports/EHR/BDSI_MGI_Freeze2_Age_Sex.csv", header = T)
colnames(age.data) <- c("StudyID", "Age", "Sex")
age.data$Age <- as.numeric(as.character(age.data$Age))
age.data$StudyID <- as.character(age.data$StudyID)

#diag.data <- read.csv("/exports/EHR/data/Freeze2/BDSI_MGI_Freeze2_Diagnosis.csv", header = T)
diag.data <- read.csv("/exports/EHR/BDSI_MGI_Freeze2_Diagnosis.csv", header = T)

diag.data$DaysSinceBirth <- as.numeric(as.character(diag.data$DaysSinceBirth))
diag.data$StudyID <- as.character(diag.data$StudyID)
diag.data$Description <- as.character(diag.data$Description)
diag.data$Description <- toupper(diag.data$Description)

# Drop Source Column and filter out duplicates, we are not interested in this column at all
diag.data$Source <-  NULL
diag.data <-  distinct(diag.data)

#lab.data <- read.csv("/exports/EHR/data/Freeze2/BDSI_MGI_Freeze2_Labs.tmp.csv", header = T, quote = "")

# Step 2: Combine age.data & diag.data to create full data frame
# Matching methods: Match by StudyID for age and diag, Match by Study ID & DaysSB for lab data
temp0.data <- diag.data %>% left_join(age.data, by = "StudyID")
temp0.data <- distinct(temp0.data) #There were no repetiting rows

#Step 3.5: Checking the percentage of missing value for each patient
miss.data <- filter(temp0.data, is.na(DaysSinceBirth))

miss.data <- as.data.frame(table(miss.data$StudyID))
colnames(miss.data) <- c("StudyID","Freq")
miss.data$StudyID <- as.character(miss.data$StudyID)
encounter.data <- as.data.frame(table(temp0.data$StudyID))
colnames(encounter.data) <- c("StudyID", "Freq")
encounter.data$StudyID <- as.character(encounter.data$StudyID)
encount.data <- miss.data %>% left_join(encounter.data, by = "StudyID")
colnames(encount.data) <- c("StudyID","MissingEncounters","TotalEncounters")
encount.data = encount.data %>% mutate(PercentMissing = round((MissingEncounters/TotalEncounters)*100,2))

#Omit rows with NA's
temp1.data <- na.omit(temp0.data) #We omit this missing value because the percentage of missing value is very small
# in the overall data set. The removal of the missing data should not deviate/alter our model by too much

# Step 3: Missing Data Imputation for DaysSinceBirth using MICE, missForest and Hmisc

# Step 4: Joining Lab Data to our Data Frame
# temp0.5.data <- temp0.data %>% left_join(lab.data, by = c("StudyID", "DaysSinceBirth"))
# temp0.5.data$UNIT <- as.character(temp0.5.data$UNIT)
# temp0.5.data$RANGE <- as.character(temp0.5.data$RANGE)
# temp0.5.data$RESULT_NAME <- as.character(temp0.5.data$RESULT_CODE)


# Step 5: Create columns/covariates we are interested in

# Create minimum column that holds the minimum DaysSinceBirth value for each patient
minimum.data <-  temp1.data %>% group_by(StudyID) %>%  summarise(Minimum = min(DaysSinceBirth, na.rm = TRUE))
temp1.data <-  temp1.data %>% left_join(minimum.data, by = "StudyID")
prior.data <-  temp1.data %>% mutate(DaysSinceFirstSeen = DaysSinceBirth - Minimum)

#create follow up which is maximum of days since first seen
prior.data <- prior.data %>% group_by(StudyID) %>% mutate(FollowUp = max(DaysSinceFirstSeen, na.rm = TRUE) )

#remove people with 0 follow up, meaning they were only in the UMICH system for just 1 day
prior.data <- prior.data %>% filter(FollowUp != 0)
# Remove patients who was diagnosed with heartfailure on their first ever diagnosis since they already had HF
y <- prior.data %>% filter(DaysSinceFirstSeen == 0 & (grepl("428", Code, fixed = FALSE) | grepl("I50", Code, fixed = FALSE)))
excludeID.vec <- as.vector(y$StudyID) #217 people with inital heart failure in their first diagnosis
temp2.data <- prior.data %>% filter(!StudyID %in% excludeID.vec)

#Create temp heart failure data set including only patients with heart failure
heart.fail.temp <- temp2.data %>% filter((grepl("428.", Code, fixed = FALSE) | grepl("I50", Code, fixed = FALSE)))
heart.fail.temp <- heart.fail.temp %>% group_by(StudyID)
# We have two Codes that were unfortunately pulled from temp2.data with which we are not interested in
# Codes = T84.428A & T84.428D
heart.fail.temp <- heart.fail.temp %>% filter(!(Code == "T84.428A" | Code == "T84.428D"))
heart.fail.temp <- heart.fail.temp %>% mutate(FirstEncount_HF = min(DaysSinceBirth, na.rm = TRUE))
heart.fail.temp <- heart.fail.temp %>% mutate(HF_FollowUp = max(DaysSinceBirth,na.rm = TRUE) - FirstEncount_HF)
#Safety Check, there should be 13910 observations in heart.fail.temp. This holds.

# Visualize of the Age and DaysSinceBirth and FirstHeartEncounter for Heart Failure Patients
# ggplot(r, aes(x = factor(0), y = FirstEncounter)) + geom_violin(color = "gold2",fill = "gold1") + geom_boxplot(fill = "blue") + xlab("") + ylab("Days Til First Heart Failure") +
#   scale_x_discrete(breaks = NULL) + coord_flip())
ggplot(heart.fail.temp, aes(Age*365)) + geom_histogram(bins=100, fill = "steelblue2")
ggplot(r, aes(FirstEncounter/365.25)) +  geom_histogram(bins = 100, fill = "steelblue2") +
  theme(axis.text.x = element_text(size = 25, color = "black"), axis.text.y = element_text(size = 25, color = "black"), axis.title = element_text(size=20))

# Baseline is Age 50 which is approximately 50*365 = 18250 DaysSinceBirth
# Exclude people who do not have any encounters past age 50 or ~18250
temp3.data <- temp2.data %>% filter(max(DaysSinceBirth) < 14235)
young.IDs <- as.vector(temp3.data$StudyID)
temp4.data <- temp2.data %>% filter(!StudyID %in% young.IDs)

temp4.5.data <- temp4.data %>% filter(DaysSinceBirth >= 14235)
temp4.6.data <- temp4.5.data %>% mutate(StudyFollowUp = max(DaysSinceBirth) - min(DaysSinceBirth))

# Exclude people who have been diagnosed with heart failure by age 50 (18250)
# temp5.data <- temp4.6.data %>% filter((grepl("428", Code, fixed = FALSE) | grepl("I50", Code, fixed = FALSE)) & DaysSinceBirth <= 14235)
# already.HF <- as.vector(temp5.data$StudyID)
# temp6.data <- temp4.data %>% filter(!StudyID %in% already.HF)

#get unique people with follow up
followup.data <- as.data.frame(temp4.6.data$StudyID) %>% bind_cols(as.data.frame(temp6.data$StudyFollowUp))
colnames(followup.data) <- c("StudyID", "FollowUp")
followup.data <- followup.data %>% distinct()

# Exclude People who have had Heart Failure at their first encounter, this was done in line 82
ggplot(followup.data, aes(FollowUp)) + geom_histogram(bins = 50, fill = "steelblue3") + xlab("") + ylab("Follow Up Time") +
  scale_x_discrete(breaks = NULL)

ggplot(followup.data, aes(FollowUp)) + geom_density(fill = "steelblue3", color = "black") + xlab("Number of Patients") + ylab("Follow Up Time") +
  scale_x_discrete(breaks = NULL)

# Exclude people that have follow up time that is less than 2 years
temp5.data <- temp4.6.data %>% filter(StudyFollowUp < 730)
low.followup <- as.vector(temp5.data$StudyID)
temp6.data <- temp4.6.data %>% filter(!StudyID %in% low.followup)

#Get everyone's first year of follow up and put it in a separate data frame
#first create value of everyone's minimum days since first seen
firstday.data <-  temp6.data %>% group_by(StudyID) %>%  summarise(FirstDay = min(DaysSinceFirstSeen, na.rm = TRUE))
temp6.data <-  left_join(temp6.data, firstday.data, by = "StudyID")
yearone <- temp6.data %>% filter((DaysSinceFirstSeen - FirstDay) <= 365)
#So year one contains everyone's data from their first year of when our study started (or their first year if greater than age 49). 
#Good thing is that everyone in the year one data set
#matches to everyone in the big data set-- so everyone has data in their first year

#Delete year one rows from temp6.data to create analysis data
analysis.data <- temp6.data %>% anti_join(yearone, by = c("StudyID", "DaysSinceBirth"))

interview.data <- yearone %>% group_by(StudyID, DaysSinceBirth) %>% 
  mutate(
    HeartFail_Gen = (grepl("428.", Code, fixed = FALSE) | grepl("I50.", Code, fixed = FALSE)),
    Female = "F" %in% Sex,
    #Arrythmia = grepl("427.9", Code, fixed = FALSE),
    Arrythmia = ("427.9" %in% Code | "I49.9" %in% Code),
    # Hyper.Malig = "401.0" %in% Code,
    # Hyper.Benign = "401.1" %in% Code,
    # Hyper.Unspec = "401.9" %in% Code,
    Hypertension = (grepl("401.", Code, fixed = FALSE) | "I10" %in% Code),
    CardioMyopathy = (grepl("425", Code, fixed = FALSE)| grepl("I42.", Code, fixed = FALSE)),
    Emphysema = (grepl("492.", Code, fixed = FALSE)|grepl("J43.", Code, fixed = FALSE)),
    Thyrotoxicosis = (grepl("242.", Code, fixed = FALSE)|grepl("E05.", Code, fixed = FALSE)),
    #Myocarditis = grepl("422", Code, fixed = FALSE),
    #Obesity = (grepl("278.",Code, fixed = FALSE) | grepl("E66.", Code, fixed = FALSE)),
    Morbid.Obesity = ("278.01" %in% Code|"E66.01" %in% Code),
    Overweight = ("278.02" %in% Code|"E66.3" %in% Code),
    #Obesity.Hypoventilation = "278.03" %in% Code,
    Hypercholesterol = ("272.0" %in% Code | "E78.0" %in% Code),
    #Depression = (grepl("311", Code, fixed = FALSE)|grepl("F32.9", Code, fixed = FALSE)),
    Type2DiabMel = ("250.00" %in% Code | "E11.9" %in% Code),
    Hepatitis_C = ("B19.20" %in% Code | "070.70" %in% Code),
    CKD_Mild = ("585.2" %in% Code | "N18.2" %in% Code ),
    CKD_Moderate = ("585.3" %in% Code | "N18.3" %in% Code),
    CKD_Severe = ("585.4" %in% Code | "N18.4" %in% Code),
    CKD_Failure = ("585.5" %in% Code | "N18.5" %in% Code),
    Sleep_Disorder = (grepl("780.5", Code, fixed = FALSE) | grepl("G47.", Code, fixed = FALSE)),
    Phobic_Anxiety = (grepl("300.2", Code, fixed = FALSE | grepl("F40.", Code, fixed = FALSE)))
    #Alc_Depen_Syn = (grepl("303.", Code, fixed = FALSE) | grepl("F10.2", Code, fixed = FALSE))
  )

#Identify Patient ID's that have heart failure encounter in the interview data, then eliminate from interview and analysis data
IDs.HF.y1 <- interview.data %>% filter(grepl("428.", Code, fixed = FALSE) | grepl("I50.", Code, fixed = FALSE))
IDs.HF.vec <- as.vector(IDs.HF.y1$StudyID)
interview1.data <- interview.data %>% filter(!StudyID %in% IDs.HF.vec)

analysis1.data <- analysis.data %>% filter(!StudyID %in% IDs.HF.vec)

# Create People's study time for follow up in our survival analysis
#analysis1.data <- analysis1.data %>% group_by(StudyID)  %>% mutate(SAFollowUP = max(DaysSinceBirth) - min(DaysSinceBirth))

#Calclulate first event of HF
#subset of analysis1.data to start with
# Grabs everyone that has heart failure and creates column that has their first encounter day for heart failure
t <-  filter(analysis1.data, grepl("428.", Code, fixed = FALSE) | grepl("I50.", Code, fixed = FALSE))
t <- t %>% group_by(StudyID) %>% mutate(FirstHF = min(DaysSinceBirth))
#create Start Date = FirstDay + 365
analysis1.data <-  analysis1.data %>% group_by(StudyID) %>% mutate(Start = min(DaysSinceBirth))
analysis1.data <-  analysis1.data %>% group_by(StudyID) %>% mutate(End = max(DaysSinceBirth))

End.data <-  analysis1.data %>% select(StudyID, Start, End) %>% distinct() 
HF.EndData <-  t %>% select(StudyID, FirstHF) %>% distinct()
End.data <-  End.data %>% left_join(HF.EndData, by = "StudyID")
Time.data <-  End.data %>% mutate(EHR.End = ifelse(is.na(FirstHF), End, FirstHF)) # Create Data Frame that consists of End of Study Time and First HF Encounter if exists, otherwise NA value

interview2.data <- interview1.data %>% ungroup() %>% select(StudyID, Code, Female:Phobic_Anxiety)
interview2.data <- interview2.data %>% group_by(StudyID) %>% summarize_at(vars(Female:Phobic_Anxiety), any)

Time.data <- Time.data %>% select(StudyID, Start, EHR.End)

SA.Data <- interview2.data %>% left_join(Time.data, by = "StudyID")


#Create Heart_Fail_Gen Column 
heart.fail <- analysis1.data %>% select(StudyID, Code) %>% group_by(StudyID) %>% 
  mutate(
    HeartFail_Gen = (grepl("428.", Code, fixed = FALSE) | grepl("I50.", Code, fixed = FALSE))
  )
heart.fail <- heart.fail %>% group_by(StudyID) %>% summarize_at(vars(HeartFail_Gen), any)

# Create Survival Analysis Data Set with columns StudyID ~Covariates Start End Y-Outcome
SA.Data <- SA.Data %>% left_join(heart.fail, by = "StudyID")

#We called it peanuts cause we were losing our sanity, peanuts has the difference time for start and end times in the original SA data set
peanuts <- SA.Data %>% group_by(StudyID) %>% 
  mutate(
    Difference = EHR.End - Start
  )

peanuts1 <- peanuts %>% filter(Difference == 0) # Scrape all the people that have a difference of 0
peanuts.ID <- peanuts1$StudyID
SA.Data <- SA.Data %>% filter(!StudyID %in% peanuts.ID) #Filter our SA Data Frame of these people, we do not want them in our analysis


# Initialize survival package to begin survival analysis
#*********************************************** SURVIVAL ANALYSIS **********************************#
library(survival)

# Run Surv function
SA.Data <- SA.Data %>% mutate(Followup = EHR.End - Start)
SA.Data <- SA.Data %>% mutate(Age = as.integer(Start/365.25))
SA.Data <- SA.Data %>% mutate(diab_age = Type2DiabMel*Age)

S <- Surv(time = SA.Data$Followup, event = SA.Data$HeartFail_Gen)
#SA.Formula.Data = SA.Data %>% ungroup() %>% select(Female:Start)
#cox.model <- coxph(as.formula(paste0("S~", paste0(colnames(SA.Formula.Data),collapse = "+"))), data = SA.Data)

cox.model.female <- coxph(S~+Arrythmia+Hypertension+CardioMyopathy+Emphysema+Morbid.Obesity+Hypercholesterol+Type2DiabMel+Hepatitis_C+CKD_Mild+CKD_Moderate+CKD_Severe+CKD_Failure+Sleep_Disorder+Age+diab_age, subset = (Female == TRUE) ,data = SA.Data)
cox.model.male <- coxph(S~+Arrythmia+Hypertension+CardioMyopathy+Emphysema+Morbid.Obesity+Hypercholesterol+Type2DiabMel+Hepatitis_C+CKD_Mild+CKD_Moderate+CKD_Severe+CKD_Failure+Sleep_Disorder+Age+diab_age, subset = (Female == FALSE) ,data = SA.Data)

library(broom)
library(lazyeval)
library(forestmodel)
library(survminer)
forest_model(cox.model.male, format_options = list(colour = "steelblue2", shape = "|", text_size = 5, banded = TRUE))
forest_model(cox.model.female, format_options = list(colour = "firebrick2", shape = "|", text_size = 5, banded = TRUE))

model.female <- cox.model.female
model.male <- cox.model.male

plot(cox.zph(cox.model, transform="km", global=TRUE), col = c("firebrick2", "steelblue3"), cex = 2, cex.axis = 2)
plot(cox.zph(model.female, transform="km", global=TRUE))
plot(cox.zph(model.male, transform="km", global=TRUE), col = c("firebrick2", "steelblue3"))


#ggforest(model.female, data = SA.Data, cpositions = c(0.02,
#0.22, 0.4))
#survival distribution for baseline people (person with no risk factors)
#assumes all covariates are 0

plot(survfit(model.female, conf.type = "none"), 
     xscale = 365.25,
     xlab = "Years of FollowUp",
     ylab = "Proportion Free of Heart Failure",
     main = "Baseline Survival Distribution",
     col = "firebrick2")
lines(survfit(model.male, conf.type = "none"), col = "steelblue3")
legend(0, .50, c("Female", "Male"), col= c("firebrick2", "steelblue3"),lwd=2)

#####
#kaplan meier curve
#"average" for all covariates
plot(survfit(S~1), 
     xscale = 365.25,
     xlab = "Years of Follow Up",
     ylab = "Proportion Free of Heart Failure",
     main = "Kaplan Meier Curve ", las = 1,
     col = c("steelblue3","firebrick2","firebrick2"))
# objNpsurv <- npsurv(formula = S ~ Female, data = SA.Data)
# survplot(objNpsurv, label.curves=list(keys="lines"), xlab = "Days of Follow Up", 
#          ylab = "Proportion Free of Heart Failure", col = c(1,12), 
#          type = "kaplan-meier", n.risk = TRUE, 
#          y.n.risk='auto', adj.n.risk = 0.5) 
#get survival curve for specific model
newdata <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata2 <- data.frame(
  #Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = TRUE,
  Emphysema = FALSE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60) 
newdata3 <- data.frame(
  Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = TRUE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60) 
newdata4 <- data.frame(
  Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = TRUE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)
newdata5 <- data.frame(
  Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = TRUE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata6 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = TRUE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata7 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = TRUE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = TRUE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata8 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = FALSE,
  CKD_Moderate = TRUE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata9 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = FALSE,
  CKD_Moderate = FALSE,
  CKD_Severe = TRUE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata10 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = FALSE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = TRUE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)
newdata11 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata12 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = TRUE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata13 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = TRUE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 60)

newdata14 <- data.frame(
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  #Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  #Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Hepatitis_C = FALSE,
  #Phobic_Anxiety = FALSE,
  diab_age = FALSE,
  Age = 61)

original.person <- survfit(model.male, newdata, se.fit=TRUE, conf.type = "none", 
                           individual=FALSE, censor=TRUE)
cardiomyoapthy.p1 <- survfit(model.male, newdata2, se.fit=TRUE, conf.type = "none", 
                             individual=FALSE, censor=TRUE)
diabetes <- survfit(model.male, newdata3, se.fit=TRUE, conf.type = "none", 
                    individual=FALSE, censor=TRUE)
hypertension <- survfit(model.male, newdata4, se.fit=TRUE, conf.type = "none", 
                        individual=FALSE, censor=TRUE)
emphysema <- survfit(model.male, newdata5, se.fit=TRUE, conf.type = "none", 
                     individual=FALSE, censor=TRUE)
morbid.obesity <- survfit(model.male, newdata6, se.fit=TRUE, conf.type = "none", 
                          individual=FALSE, censor=TRUE)
ckd.moderate <- survfit(model.male, newdata8, se.fit=TRUE, conf.type = "none", 
                        individual=FALSE, censor=TRUE)
ckd.severe <- survfit(model.male, newdata9, se.fit=TRUE, conf.type = "none", 
                      individual=FALSE, censor=TRUE)
ckd.failure <- survfit(model.male, newdata10, se.fit=TRUE, conf.type = "none", 
                       individual=FALSE, censor=TRUE)
sleep.disorder <- survfit(model.male, newdata12, se.fit=TRUE, conf.type = "none", 
                          individual=FALSE, censor=TRUE)
hep.c <- survfit(model.male, newdata13, se.fit=TRUE, conf.type = "none", 
                 individual=FALSE, censor=TRUE)
age.change <- survfit(model.male, newdata14, se.fit=TRUE, conf.type = "none", 
                      individual=FALSE, censor=TRUE)

plot(original.person, 
     xscale = 365.25,
     xlab = "Years of FollowUp",
     ylab = "Proportion Free of Heart Failure",
     main = "Survival Distribution with Varying Covariates (Male)", las = 1, col = c(1),
     cex = 2, cex.axis = 1.5, lwd = 2) 
lines(cardiomyoapthy.p1, col = "firebrick2")
lines(diabetes, col = "navy")
lines(hypertension, col = "steelblue3")
#lines(emphysema, col = "YELLOW")
#lines(morbid.obesity, col = "PURPLE")
#lines(ckd.moderate, col = "BROWN")
lines(ckd.severe, col = "seagreen4")
#lines(ckd.failure, col = "PINK")
#lines(sleep.disorder, col = "indianred4")
#lines(hep.c, col = "seagreen1")
#lines(age.change, col = "orangered")

legend("bottomleft", text.width = .05, bty = "n",  lwd=c(2,2), cex = 1.3, y.intersp=.3, 
       c("Reference","+ Hypertension", "+ Diabetes", "+ CKD Severe", "+ Cardiomyopathy"), 
       col = c("black", "steelblue3", "navy", "seagreen4", "firebrick2")) 


original.person <- survfit(model.female, newdata, se.fit=TRUE, conf.type = "none", 
                           individual=FALSE, censor=TRUE)
cardiomyoapthy.p1 <- survfit(model.female, newdata2, se.fit=TRUE, conf.type = "none", 
                             individual=FALSE, censor=TRUE)
diabetes <- survfit(model.female, newdata3, se.fit=TRUE, conf.type = "none", 
                    individual=FALSE, censor=TRUE)
hypertension <- survfit(model.female, newdata4, se.fit=TRUE, conf.type = "none", 
                        individual=FALSE, censor=TRUE)
emphysema <- survfit(model.female, newdata5, se.fit=TRUE, conf.type = "none", 
                     individual=FALSE, censor=TRUE)
morbid.obesity <- survfit(model.female, newdata6, se.fit=TRUE, conf.type = "none", 
                          individual=FALSE, censor=TRUE)
ckd.moderate <- survfit(model.female, newdata8, se.fit=TRUE, conf.type = "none", 
                        individual=FALSE, censor=TRUE)
ckd.severe <- survfit(model.female, newdata9, se.fit=TRUE, conf.type = "none", 
                      individual=FALSE, censor=TRUE)
ckd.failure <- survfit(model.female, newdata10, se.fit=TRUE, conf.type = "none", 
                       individual=FALSE, censor=TRUE)
sleep.disorder <- survfit(model.female, newdata12, se.fit=TRUE, conf.type = "none", 
                          individual=FALSE, censor=TRUE)
hep.c <- survfit(model.female, newdata13, se.fit=TRUE, conf.type = "none", 
                 individual=FALSE, censor=TRUE)
age.change <- survfit(model.female, newdata14, se.fit=TRUE, conf.type = "none", 
                      individual=FALSE, censor=TRUE)

plot(original.person, 
     xscale = 365.25,
     xlab = "Years of Follow Up",
     ylab = "Proportion Free of Heart Failure",
     main = "Survival Distribution with Varying Covariates (Female)", las = 1, col = c(1),
     cex = 2, cex.axis = 1.5, lwd = 2) 
lines(cardiomyoapthy.p1, col = "firebrick2")
legend("bottomleft", text.width = .1, bty = "n",  lwd=c(2,2), cex = 1.5, y.intersp=.5, 
       c("Reference", "+ Cardiomyopathy"), 
       col = c("black","firebrick2")) 
lines(diabetes, col = "navy")
lines(hypertension, col = "steelblue3")
#lines(emphysema, col = "YELLOW")
#lines(morbid.obesity, col = "PURPLE")
#lines(ckd.moderate, col = "BROWN")
lines(ckd.severe, col = "seagreen4")
#lines(ckd.failure, col = "PINK")
#lines(sleep.disorder, col = "indianred4")
#lines(hep.c, col = "seagreen1")
#lines(age.change, col = "orangered")


legend("bottomleft", text.width = .05, bty = "n",  lwd=c(2,2), cex = 1.2, y.intersp=.3, 
       c("+ Cardiomyopathy", "Reference", "+ CKD Severe", "+ Hypertension", "+ Diabetes"), 
       col = c("firebrick2", "black", "seagreen4","steelblue3", "navy")) 

#***************************************************************************************#



## cox.model table results
## The following code is courtesy of https://github.com/NikNakk, Repository: https://github.com/NikNakk/forestmodel
## Stack Overflow Thread: https://stackoverflow.com/questions/31289111/optimal-efficient-plotting-of-survival-regression-analysis-results


#survival distribution for baseline people (person with no risk factors)
#assumes all covariates are 0

#survival distribution for baseline people (person with no risk factors)
#assumes all covariates are 0
plot(survfit(cox.model), 
     xscale = 365.25,
     xlab = "Follow Up After Age 40",
     ylab = "Proportion Free of Heart Failure",
     main = "",
     col = c("red", "darkblue", "darkblue"),
     las = 1)

#kaplan meier curve
#"average" for all covariates
plot(survfit(S~1), 
     xscale = 365.25,
     xlab = "Years of Follow Up",
     ylab = "Proportion Free of Heart Failure",
     main = "", col = c("blue", 1, 1), las = 1)

fit <- survfit(S ~ Female, data = SA.Data)
plot(fit, lty = 2:3, 
     xscale = 365.25,
     xlab = "Years of Follow Up",
     ylab = "Proportion Free of Heart Failure",
     main = "", las = 1)
legend(50, .4, c("Male", "Female"), lty = 2:3)

fit2 <- survfit(S ~ Hypertension, data = SA.Data)
plot(fit2, lty = 2:3, 
     xscale = 365.25,
     xlab = "Years of Follow Up",
     ylab = "Proportion Free of Heart Failure",
     main = "", las = 1, col = c("steelblue2", "gold3")) 
legend(50, .4, c("No Hypertension", "Hypertension"), lty = 2:3, col = c("steelblue2","gold3"), box.col = "steelblue")

fit3 <- survfit(S ~ CardioMyopathy, data = SA.Data)
plot(fit3, lty = 2:3, 
     xscale = 365.25,
     xlab = "Years of Follow Up",
     ylab = "Proportion Free of Heart Failure",
     main = "", las = 1, col = c("steelblue2", "gold3")) 
legend(40, .3, c("No CardioMyopathy", "CardioMyopathy"), lty = 2:3, col = c("steelblue2","gold3"), box.col = "steelblue")

fit4 <- survfit(S ~ Type2DiabMel, data = SA.Data)
plot(fit4, lty = 2:3, 
     xscale = 365.25,
     xlab = "Years of Follow Up",
     ylab = "Proportion Free of Heart Failure",
     main = "", las = 1, col = c("steelblue2", "gold3")) 
legend(50, .4, c("No Type II Diabetes", "Type II Diabetes"), lty = 2:3, col = c("steelblue2","gold3"), box.col = "steelblue")

fit5 <- survfit(S ~ CKD_Severe, data = SA.Data)
plot(fit5, lty = 2:3, 
     xscale = 365.25,
     xlab = "Years of Follow Up",
     ylab = "Proportion Free of Heart Failure",
     main = "", las = 1, col = c("steelblue2", "gold3")) 
legend(50, .4, c("No Severe CKD", "Severe CKD"), lty = 2:3, col = c("steelblue2","gold3"), box.col = "steelblue")

# objNpsurv <- npsurv(formula = S ~ Female, data = SA.Data)
# survplot(objNpsurv, label.curves=list(keys="lines"), xlab = "Days of Follow Up", 
#          ylab = "Proportion Free of Heart Failure", col = c(1,12), 
#          type = "kaplan-meier", n.risk = TRUE, 
#          y.n.risk='auto', adj.n.risk = 0.5)

#Create Hypothetical Person with specific covariates for model visualization and cov. effect
newdata <- data.frame(
  Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  Hepatitis_C = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Phobic_Anxiety = FALSE,
  Age = 45)  

original.person <- plot(survfit(formula = cox.model, newdata = newdata, se.fit=TRUE, conf.int= FALSE, 
                                conf.type = "none", individual= FALSE, censor=TRUE), 
                        xscale = 365.25,
                        xlab = "Years of FollowUp",
                        ylab = "Proportion Free of Heart Failure",
                        main = "", las = 1, lty = 2) 

#This person has cardiomyopathy
newdata1 <- data.frame(
  Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = TRUE,
  Emphysema = FALSE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE ,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Phobic_Anxiety = FALSE,
  Age = 45)  

cardiomyopathy <- survfit(formula = cox.model, newdata = newdata1, se.fit=TRUE, conf.int= FALSE, 
                          conf.type = "none", individual= FALSE, censor=TRUE)

lines(cardiomyopathy)

newdata2 <- data.frame(
  Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE,
  Hypercholesterol = TRUE,
  Type2DiabMel = TRUE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Phobic_Anxiety = FALSE,
  Age = 45)

type2.diabetes <- survfit(formula = cox.model, newdata = newdata2, se.fit=TRUE, conf.int= FALSE, 
                          conf.type = "none", individual= FALSE, censor=TRUE)
lines(type2.diabetes, col = 3)

newdata3 <- data.frame(
  Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = TRUE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Phobic_Anxiety = FALSE,
  Age = 45)

hyper.tension <- survfit(formula = cox.model, newdata = newdata3, se.fit=TRUE, conf.int= FALSE, 
                         conf.type = "none", individual= FALSE, censor=TRUE)
lines(hyper.tension, col = 4)

newdata3 <- data.frame(
  Female = TRUE,
  Arrythmia = TRUE,
  Hypertension = FALSE,
  CardioMyopathy = FALSE,
  Emphysema = FALSE,
  Thyrotoxicosis = FALSE,
  Morbid.Obesity = FALSE,
  Overweight = FALSE,
  Hypercholesterol = TRUE,
  Type2DiabMel = FALSE,
  CKD_Mild = TRUE,
  CKD_Moderate = FALSE,
  CKD_Severe = FALSE,
  CKD_Failure = FALSE,
  Sleep_Disorder = FALSE,
  Phobic_Anxiety = FALSE,
  Age = 45)

hyper.tension <- survfit(formula = cox.model, newdata = newdata3, se.fit=TRUE, conf.int= FALSE, 
                         conf.type = "none", individual= FALSE, censor=TRUE)
lines(hyper.tension, col = 5)

