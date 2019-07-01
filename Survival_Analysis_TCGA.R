##############################################################
 "Survival Analysis of Breast Cancer Data from the TCGA Dataset"
author: "Ramaa Nathan"
##############################################################


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, cache=TRUE)

#The Following libraries need to be installed once.
# Load the bioconductor installer. 
# Try http:// if https:// doesn't work.
#source("https://bioconductor.org/biocLite.R")

# Install the main RTCGA package
#biocLite("RTCGA")

# Install the clinical and mRNA gene expression data packages
#biocLite("RTCGA.clinical")
#biocLite("RTCGA.mRNA")

# Test CRAN package installation:
library(tidyverse)
library(stringr) #for string manipulations
library(forcats) # for factors

library(knitr) # for kable
library(kableExtra) # for formatting tables

library(gridExtra) #for plotting in a grid

library(xtable) 
library(JM)
library(cmprsk)

#load survival libraries
library(survminer)
library(survival)

# Test RTCGA: 
library(RTCGA)
library(RTCGA.clinical)
library(RTCGA.mRNA)

#modeling
library(AICcmodavg)  #to extract AICc
library(parmsurvfit)
library(flexsurv)
library(SurvRegCensCov)


# Create the clinical data and extract data for all three cancers
clin <- survivalTCGA(BRCA.clinical, OV.clinical, GBM.clinical, 
                     extract.cols=c("admin.disease_code","patient.drugs.drug.therapy_types.therapy_type"))


# Create the clinical data
brca_clin_orig <- survivalTCGA(BRCA.clinical, 
                     extract.cols=c("patient.gender", "patient.race",
                                    "patient.ethnicity","patient.days_to_birth","patient.vital_status",
                                    "patient.drugs.drug.therapy_types.therapy_type",
                                    "patient.stage_event.pathologic_stage",
                                    "patient.stage_event.tnm_categories.pathologic_categories.pathologic_t",
                                    "patient.stage_event.tnm_categories.pathologic_categories.pathologic_n",
                                    "patient.stage_event.tnm_categories.pathologic_categories.pathologic_m"))
#sapply(brca_clin_orig,class)


## Checks for Missing and Invalid Data

# Check for number of missing data in each column
brca_clin_orig %>% sapply(function(x) sum(is.na(x)))

#Check number of rows with missing data
na_rows <- brca_clin_orig %>% apply(MARGIN=1, FUN=function(x) sum(is.na(x))) 
sum(na_rows>0)

brca_clin_orig %>% sapply(summary)
brca_clin_orig %>% filter(times < 0) %>% nrow()


## Data Cleaning and Transformation
#Next, the following cleaning and transformations are applied to the TCCGA BRCA.clinical data to get a clean and compact dataset:

#1. Rename the long variable names to short names.
#2. Filter out the 12 observations corresponding to males diagnosed with breast cancer.
#3. Filter out the 2 observations with negative "times" value.
#4. Create a "age" variable that contains the number of days at first diagnosis to age in years at first diagnosis.
#5. Create a "years_to_event" variable which is the "times"" variable   converted from days to years.
#6. Data in the pathology columns contain information on both stage and sub-stage. Transform the data to only contain the high level stage information.
#7. Modify the "therapy_type"" to contain three types - chemotherapy, hormone therapy and Other (lump all the other infrequent types into Other)
#8. Modify the "race"" to contain three types - black or african american and white (lump the other two types into Other)


#helper Functions
clean_pathologic_stage <- function(x) {
  x %>% str_replace_all(c(
    "stage iv[a-d]*"="stage4",
    "stage [i]{3}[a-d]*"="stage3",
    "stage i{2}[a-d]*"="stage2",
    "stage i{1}[a-d]*"="stage1",
    "stage x"="stageX"))
}

clean_pathologyTstage <- function(x) {
  x %>% str_replace_all(c(
    "\\s*tx"="tx",
    "\\s*t1[a-z]*"="t1", 
    "\\s*t2[a-z]*"="t2",
    "\\s*t3[a-z]*"="t3",
    "\\s*t4[a-z]*"="t4"))
}

clean_pathologyMstage <- function(x) {
  x %>% str_replace_all(c("cm0\\s+\\(i[\\+,-]\\)"="cm0"))
}

clean_pathologyNstage <- function(x) {
  x %>% str_trim(.) %>% 
    str_replace_all(c(
      "nx"="nx", 
      "n1[a-z]*"="n1", 
      "n2[a-z]*"="n2",
      "n3[a-z]*"="n3",
      "n4[a-z]*"="n4",
      "n0"="n0",
      "n0\\s+\\([a-z]*[\\+|-]\\)"="n0"))
}


#transform the data
brca_clin <- brca_clin_orig %>% 
  rename(gender=patient.gender,
         race=patient.race,
         ethnicity=patient.ethnicity,
         vital_status=patient.vital_status,
         therapy_type=patient.drugs.drug.therapy_types.therapy_type,
         pathologic_stage=patient.stage_event.pathologic_stage,
         pathologyTstage=patient.stage_event.tnm_categories.pathologic_categories.pathologic_t,
         pathologyNstage=patient.stage_event.tnm_categories.pathologic_categories.pathologic_n,
         pathologyMstage=patient.stage_event.tnm_categories.pathologic_categories.pathologic_m) %>%
  filter(gender == "female") %>% 
  filter(times > 0) %>%
  mutate(age=abs(as.numeric(patient.days_to_birth))/365,
         therapy_type = ifelse(is.na(therapy_type),"No Info",therapy_type),
         therapy_type = fct_lump(therapy_type,3),
         race = fct_lump(race,2),
         pathologic_stage=str_trim(pathologic_stage),
         pathologyTstage=str_trim(pathologyTstage),
         pathologyNstage=str_trim(pathologyNstage),
         pathologyMstage=str_trim(pathologyMstage),
         pathologic_stage = clean_pathologic_stage(pathologic_stage),
         pathologyTstage = clean_pathologyTstage(pathologyTstage),
         pathologyNstage = clean_pathologyNstage(pathologyNstage),
         pathologyMstage = clean_pathologyMstage(pathologyMstage),
         years_to_event=times/365,
         agecat=cut(age, breaks=c(0, 40, 60, Inf), labels=c("young", "middle", "old"))
         )

#Convert specified columns from character to factor type.
convert_to_factor <- c("ethnicity", "pathologic_stage",
                       "pathologyTstage", "pathologyMstage","pathologyNstage")
brca_clin <- brca_clin %>% mutate_at(convert_to_factor,factor)

#remove unnecessary columns
brca_clin <- brca_clin %>% 
  dplyr::select(-starts_with("patient"))
#names(brca_clin)

#Verify the class types
#sapply(brca_clin,class)

#What are the dimensions of the brca_clin dataset?
#dim(brca_clin)


## Data Visualizations

### Histograms
**Age at First Diagnosis**

brca_clin %>% filter(!is.na(age)) %>% 
  ggplot(aes(x=age)) +
    geom_histogram(aes(y=..density..),color="black",fill="white")+
    geom_density(alpha=0.2,fill="#FF6666") +
    geom_vline(aes(xintercept=mean(age)), color="blue", linetype="dashed", size=1) +
    labs(title="Age at First Diagnosis")

brca_clin %>% filter(!is.na(age)) %>% 
  ggplot(aes(x=age)) +
    geom_histogram(color="black",fill="white")+
    #geom_density(alpha=0.2,fill="#FF6666") +
    geom_vline(aes(xintercept=mean(age)), color="blue", linetype="dashed", size=1) +
    labs(title="Age at First Diagnosis")
    

brca_clin %>% filter(!is.na(times)) %>%
  ggplot(aes(x=years_to_event)) +
    geom_histogram(color="black",fill="white") +
    geom_vline(aes(xintercept=mean(years_to_event)), color="blue", linetype="dashed", size=1) +
    labs(title="Time to Event")

brca_clin %>% dplyr::select(age) %>% summary()


#Contingency table
xtabs(~admin.disease_code+patient.vital_status,data=clin) %>% addmargins()


# Right Censoring Plot by Disease
clin <- clin %>%
  mutate(years_to_event=times/365) %>%
  filter(!is.na(admin.disease_code)) %>%
  arrange(admin.disease_code) 

clin %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, 
           yend = index, 
           colour = admin.disease_code,
           shape = factor(patient.vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - By Disease") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Disease", labels = c("BRCA (Breast Cancer)", "GBM (Glioblastoma Multiforme)", "OV (Ovarian Cancer)"))




###Censoring and Event Plots  by Age Category

#Contingency table
xtabs(~agecat+vital_status,data=brca_clin) %>% addmargins()

brca_clin_age <- brca_clin %>%
  filter(!is.na(agecat)) %>%
  arrange(agecat) 

brca_clin_age %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, #times, 
           yend = index, 
           colour = agecat,
           shape = factor(vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - BRCA by Age Categories") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Age Categories", labels = c("Young", "Middle", "Old"))

 

###Censoring and Event Plots  by Race

#Contingency table
xtabs(~race+vital_status,data=brca_clin) %>% addmargins()

brca_clin_race <- brca_clin %>%
  filter(!is.na(race)) %>%
  arrange(race) 

brca_clin_race %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, #times, 
           yend = index, 
           colour = race,
           shape = factor(vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - BRCA by Race") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Race", labels = c("Black or African American", "White", "Other"))


### Censoring and Event Plots by  Ethnicity

#Contingency table
xtabs(~ethnicity+vital_status,data=brca_clin) %>% addmargins()

brca_clin_ethnicity <- brca_clin %>%
  filter(!is.na(ethnicity)) %>%
  arrange(ethnicity)

brca_clin_ethnicity %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, #times, 
           yend = index, 
           colour = ethnicity,
           shape = factor(vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - BRCA by Ethnicity") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Ethnicity", labels = c("Hispanic or Latino", "Not Hispanic or Latino"))



### Censoring and Event Plots by Therapy Type

#Contingency table
xtabs(~therapy_type+vital_status,data=brca_clin) %>% addmargins()

brca_clin_therapy <- brca_clin %>%
  # mutate(therapy_type = ifelse(is.na(therapy_type),"No Info",therapy_type),
  #        therapy_type = fct_lump(therapy_type,3)) %>%
  arrange(therapy_type)

brca_clin_therapy %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, #times, 
           yend = index, 
           colour = therapy_type,
           shape = factor(vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - BRCA by Therapy Type") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Therapy Type", labels = c("Chemotherapy","Hormone Therapy", "No Info", "Other"))


### Censoring and Event Plots  by Pathological Stage
#Contingency table
xtabs(~pathologic_stage+vital_status,data=brca_clin) %>% addmargins()

brca_clin_stage <- brca_clin %>%
  filter (!is.na(pathologic_stage)) %>%
  arrange(pathologic_stage)

brca_clin_stage %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, #times, 
           yend = index, 
           colour = pathologic_stage,
           shape = factor(vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - BRCA by Pathologic (Cancer) Stage") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Cancer Stage", labels = c("Stage1","Stage2","Stage3","Stage4","StageX"))


### Censoring and Event Plots  by T stage
#Contingency table
xtabs(~pathologyTstage+vital_status,data=brca_clin) %>% addmargins()

brca_clin_Tstage <- brca_clin %>%
  filter (!is.na(pathologyTstage)) %>%
  arrange(pathologyTstage)

brca_clin_Tstage %>%
  arrange(pathologyTstage) %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, #times, 
           yend = index, 
           colour = pathologyTstage,
           shape = factor(vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - BRCA by Tumor (T) stage") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Tumor Stage", labels = c("t1","t2","t3","t4","tx"))


### Censoring and Event Plots  by M Stage

#Contingency table
xtabs(~pathologyMstage+vital_status,data=brca_clin) %>% addmargins()

brca_clin_Mstage <- brca_clin %>%
  filter (!is.na(pathologyMstage)) %>%
  arrange(pathologyMstage)

brca_clin_Mstage %>%
  arrange(pathologyMstage) %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, #times, 
           yend = index, 
           colour = pathologyMstage,
           shape = factor(vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - BRCA by Metastasis (M) stage") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Metastatis Stage", labels = c("cm0","m0","m1","mx"))


### Censoring and Event Plots  by N  Stage

#Contingency table
xtabs(~pathologyNstage+vital_status,data=brca_clin) %>% 
  addmargins() %>% 
  kable() %>% 
  kable_styling(full_width=F,bootstrap_options=c("bordered")) 

brca_clin_Nstage <- brca_clin %>%
  filter (!is.na(pathologyNstage)) %>%
  arrange(pathologyNstage)

brca_clin_Nstage %>%
  arrange(pathologyNstage) %>%
  mutate(index=1:n()) %>% 
  ggplot(
       aes(xend = 0, 
           y = index, 
           x = years_to_event, #times, 
           yend = index, 
           colour = pathologyNstage,
           shape = factor(vital_status))) + 
  geom_segment() + 
  geom_point() +
  ggtitle("Right Censoring in TCGA - BRCA by Lymph Nodes (N) stage") +
  labs(x="Years to Event", y="Subjects") +
  scale_shape_discrete(name = "Status", labels = c("Censored","Event")) +
  scale_color_discrete(name = "Lymph Nodes Stage", labels = c("n0","n1","n2","n3","nx"))
 


## Cancer Types
# Model the Kaplan Meier Survival Curve
surv_can_obj <- with(clin, Surv(years_to_event, patient.vital_status))
surv_can_fit<- survfit(surv_can_obj ~ admin.disease_code,data=clin)
surv_can_fit

#Plot the survival curves for the different groups
ggplot_age_surv <- ggsurvplot(surv_can_fit,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves")

#Plot the cumulative hazard rates for the different groups
ggplot_age_cumhaz <- ggsurvplot(surv_can_fit,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves",
           risk.table=TRUE,
           cumevents=TRUE,
           fontsize=5,
           tables.height=0.3,
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_age_surv, ggplot_age_cumhaz), ncol=1)



##Time to Event only


### Kaplan Meier Survival Curves 

# Model the Kaplan Meier Survival Curve
surv_obj <- with(brca_clin, Surv(years_to_event, vital_status))
surv_fit<- survfit(surv_obj ~ 1,data=brca_clin)
surv_fit

#Plot the survival curves for the different groups
ggplot_age_surv <- ggsurvplot(surv_fit,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves")

#Plot the cumulative hazard rates for the different groups
ggplot_age_cumhaz <- ggsurvplot(surv_fit,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves",
           risk.table=TRUE,
           cumevents=TRUE,
           fontsize=5,
           tables.height=0.3,
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_age_surv, ggplot_age_cumhaz), ncol=1)

### Parametric Modeling
#Fit Weibull 
surv_wei <- survreg(surv_obj ~ 1, data=brca_clin, dist="weibull")
surv_wei_conv <- ConvertWeibull(surv_wei)
surv_wei_conv

#Fit Exponential
surv_exp <- update(surv_wei,dist="exponential")
summary(surv_exp)
lambda <- exp(-1.0 * unname(surv_exp$icoef))
paste("Constant Hazard Rate for Exponential Distribution without any predictor:",round(lambda,2), " events per year")

# Compare the weibull and exponential distributions
cmp_wei_exp <- anova(surv_wei, surv_exp)
p_val <- cmp_wei_exp$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp$"-2*LL"[1]
deviance_exp <- cmp_wei_exp$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)

ifelse(better_dist == "weibull", surv_wei_conv, summary(surv_exp))


## Age Category

### Kaplan Meier Survival Curves for Age Categories
# Model the Kaplan Meier Survival Curve
surv_obj_age <- with(brca_clin_age, Surv(years_to_event, vital_status))
surv_fit_age<- survfit(surv_obj_age ~ agecat,data=brca_clin_age)#times,
surv_fit_age

#Plot the survival curves for the different groups
ggplot_age_surv <- ggsurvplot(surv_fit_age,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves by Age Categories")

#Plot the cumulative hazard rates for the different groups
ggplot_age_cumhaz <- ggsurvplot(surv_fit_age,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves by Age Categories",
           risk.table=TRUE,
           cumevents=TRUE,
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_age_surv, ggplot_age_cumhaz), ncol=1)

**Compare Survival Curves**
#Check for Difference between the survival curves
survdiff(surv_obj_age ~ agecat,data=brca_clin_age)


**Check for  pairs of curves that are different**
surv_pairdiff_age <- pairwise_survdiff(Surv(years_to_event,vital_status)~agecat,data=brca_clin)
surv_pairdiff_age
symnum(surv_pairdiff_age$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
   symbols = c("****", "***", "**", "*", "+", " "),
   abbr.colnames = FALSE, na = "")


### Cox Proportional Hazard (PH) Model for Age

# Fit the model
cox_age <- coxph(surv_obj_age ~ age, data=brca_clin_age)
#plot the Confidence Intervals of the Hazard Ratio - no need here
#ggforest(cox_age, data=brca_clin_age)

#determining the survival rate for a new age value
#summary(survfit(cox_age,newdata=data.frame(age=60)))

# Compute the base hazard rate
#h0 <- basehaz(cox_age, centered=FALSE)

**Check for Proportional Hazard Assumptions**
test_ph_age <- cox.zph(cox_age)
test_ph_age

ggcoxzph(test_ph_age,ncol=2)


#**Interpret Model Output**
summary(cox_age)


### Parametric Modeling - Age as Categorical
#Fit Weibull 
surv_wei_age <- survreg(surv_obj_age ~ agecat, data=brca_clin_age, dist="weibull")
surv_wei_age_conv <- ConvertWeibull(surv_wei_age)
summary(surv_wei_age)


#Fit Exponential
surv_exp_age <- update(surv_wei_age,dist="exponential")
summary(surv_exp_age)
#lambda <- exp(-1.0 * unname(surv_exp_age$coefficients[1]))
#paste("Constant Hazard Rate for Exponential Distribution with Age as predictor:",round(lambda,2), " events/year")


# Compare the weibull and exponential distributions
cmp_wei_exp_age <- anova(surv_wei_age, surv_exp_age)
kable(as.data.frame(cmp_wei_exp_age))
p_val <- cmp_wei_exp_age$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp_age$"-2*LL"[1]
deviance_exp <- cmp_wei_exp_age$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)


ifelse(better_dist == "weibull", print(surv_wei_age_conv), summary(surv_exp_age))



### Parameter Modeling Age as continuous
#Fit Weibull 
surv_wei_age_c <- survreg(surv_obj_age ~ age, data=brca_clin_age, dist="weibull")
surv_wei_age_c_conv <- ConvertWeibull(surv_wei_age_c)
summary(surv_wei_age_c)

#Fit Exponential
surv_exp_age <- update(surv_wei_age,dist="exponential")
summary(surv_exp_age)
#lambda <- exp(-1.0 * unname(surv_exp_age$coefficients[1]))
#paste("Constant Hazard Rate for Exponential Distribution with Age as predictor:",round(lambda,2), " events/year")


**Analyse the weibull distribution using Age as a continuous variable)**
age_wei <- survreg(surv_obj_age ~ age, data=brca_clin_age, dist="weibull")
ConvertWeibull(age_wei)

We get the exact same results as was obtained with Cox PH model. This further validates the fits of the models. 


## Race

### Kaplan Meier Survival Curves
# Model the Kaplan Meier Survival Curve
surv_obj_race <- with(brca_clin_race,Surv(years_to_event, vital_status))
surv_fit_race <- survfit(surv_obj_race ~ race,data=brca_clin_race)#times,
surv_fit_race

#Plot the survival curves for the different groups
ggplot_race_surv <- ggsurvplot(surv_fit_race,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves by Race",
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

#Plot the cumulative hazard rates for the different groups
ggplot_race_cumhaz <- ggsurvplot(surv_fit_race,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves by Race",
           risk.table=TRUE,
           cumevents=TRUE,
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)


arrange_ggsurvplots(list(ggplot_race_surv, ggplot_race_cumhaz), ncol=1)
#Check for Difference between the survival curves
survdiff(surv_obj_race~race,data=brca_clin_race)



### Cox Proportional Hazard Model for Race
# Fit the model
cox_race <- coxph(surv_obj_race ~ race, data=brca_clin_race)
**Check for Proportional Hazard Assumptions**
test_ph_race <- cox.zph(cox_race)
test_ph_race

ggcoxzph(test_ph_race, data=brca_clin_race)
#plot(test_ph_race)

#ggcoxdiagnostics(cox_race,type="schoenfeld", data=brca_clin_race)



**Interpret the model**
summary(cox_race)

#Plot the CIs of the hazard ratio
ggforest(cox_race,data=brca_clin_race)

#Plot Survival curves for each group
#ggadjustedcurves(cox_race, var="race", data = brca_clin_race)



### Parametric modeling - Race
#Fit Weibull 
surv_wei_race <- survreg(surv_obj_race ~ race, data=brca_clin_race, dist="weibull")
surv_wei_race_conv <- ConvertWeibull(surv_wei_race)
summary(surv_wei_race)

#Fit Exponential
surv_exp_race <- update(surv_wei_race,dist="exponential")
summary(surv_exp_race)


## Ethnicity

### Kaplan Meier Survival Curves for Ethnicity
# Model the Kaplan Meier Survival Curve
surv_obj_ethnicity <- with(brca_clin_ethnicity, Surv(years_to_event, vital_status))
surv_fit_ethnicity<- survfit(surv_obj_ethnicity ~ ethnicity,data=brca_clin_ethnicity)#times,
surv_fit_ethnicity

#Plot the survival curves for the different groups
ggplot_ethnicity_surv <- ggsurvplot(surv_fit_ethnicity,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves by Ethnicity",
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

#Plot the cumulative hazard rates for the different groups
ggplot_ethnicity_cumhaz <- ggsurvplot(surv_fit_ethnicity,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves by Ethnicity",
           risk.table=TRUE,
           cumevents=TRUE,
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_ethnicity_surv, ggplot_ethnicity_cumhaz), ncol=1)

#Check for Difference between the survival curves
surv_fit_ethnicity <- survdiff(surv_obj_ethnicity ~ ethnicity,data=brca_clin_ethnicity)



### Cox Proportional Hazard (PH) Model for Ethnicity
# Fit the model
cox_ethnicity <- coxph(surv_obj_ethnicity ~ ethnicity, data=brca_clin_ethnicity)

### Parametric Modeling - Ethnicity
#Fit Weibull 
surv_wei_ethnicity <- survreg(surv_obj_ethnicity ~ ethnicity, data=brca_clin_ethnicity, dist="weibull")
surv_wei_ethnicity_conv <- ConvertWeibull(surv_wei_ethnicity)
summary(surv_wei_ethnicity)

#Fit Exponential
surv_exp_ethnicity <- update(surv_wei_ethnicity,dist="exponential")
summary(surv_exp_ethnicity)
#lambda <- exp(-1.0 * unname(surv_exp_age$coefficients[1]))
#paste("Constant Hazard Rate for Exponential Distribution with Age as predictor:",round(lambda,2), " events/year")


# Compare the weibull and exponential distributions
cmp_wei_exp_ethnicity <- anova(surv_wei_ethnicity, surv_exp_ethnicity)
p_val <- cmp_wei_exp_ethnicity$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp_ethnicity$"-2*LL"[1]
deviance_exp <- cmp_wei_exp_ethnicity$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)



## Therapy Type

### Kaplan Meier Survival Curves for Therapy Type
# Model the Kaplan Meier Survival Curve
surv_obj_therapy <- with(brca_clin_therapy, Surv(years_to_event, vital_status))
surv_fit_therapy<- survfit(surv_obj_therapy~therapy_type,data=brca_clin_therapy)#times,
surv_fit_therapy

#Plot the survival curves for the different groups
ggplot_therapy_surv <- ggsurvplot(surv_fit_therapy,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves by Therapy Type",
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

#Plot the cumulative hazard rates for the different groups
ggplot_therapy_cumhaz <- ggsurvplot(surv_fit_therapy,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves by Therapy Type",
           risk.table=TRUE,
           cumevents = TRUE,
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_therapy_surv, ggplot_therapy_cumhaz), ncol=1)

#Check for Difference between the survival curves
survdiff(surv_obj_therapy ~ factor(therapy_type,exclude=NULL),data=brca_clin_therapy)
surv_pairdiff_therapy <- pairwise_survdiff(Surv(years_to_event,vital_status)~therapy_type,data=brca_clin)
surv_pairdiff_therapy

symnum(surv_pairdiff_therapy$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
   symbols = c("****", "***", "**", "*", "+", " "),
   abbr.colnames = FALSE, na = "")


### Cox Proportional Hazard (PH) Model for Therapy
## Fit the model
cox_therapy <- coxph(surv_obj_therapy ~ therapy_type, data=brca_clin_therapy)

**Check for Proportional Hazards assumption**
#check for proportional hazard assumptions
test_ph_therapy <- cox.zph(cox_therapy)
test_ph_therapy

ggcoxzph(test_ph_therapy)
#plot(test_ph_therapy)

#ggcoxdiagnostics(cox_therapy,type="schoenfeld", data=brca_clin_therapy)


**Model Interpretation**
summary(cox_therapy)
ggforest(cox_therapy,data=brca_clin_therapy)



#**Plots of Survival Curves by Therapy Type** 
ggadjustedcurves(cox_therapy, var="therapy_type", data = brca_clin_therapy)

### Parametric modeling - Therapy
#Fit Weibull 
surv_wei_therapy <- survreg(surv_obj_therapy ~ therapy_type, data=brca_clin_therapy, dist="weibull")
surv_wei_therapy_conv <- ConvertWeibull(surv_wei_therapy)
summary(surv_wei_therapy)


#Fit Exponential
# We  cannot consider nested models here as weibull is not a good fit.
surv_exp_therapy <- survreg(surv_obj_therapy ~ therapy_type, data=brca_clin_therapy,dist="exponential")
summary(surv_exp_therapy)
lambda <- exp(-1.0 * unname(surv_exp_therapy$coefficients[1]))
paste("Constant Hazard Rate for Exponential Distribution with therapy as predictor:",round(lambda, digits=2) , "events/year")




# Compare the weibull and exponential distributions
cmp_wei_exp_therapy <- anova(surv_wei_therapy, surv_exp_therapy)
cmp_wei_exp_therapy
p_val <- cmp_wei_exp_therapy$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp_therapy$"-2*LL"[1]
deviance_exp <- cmp_wei_exp_therapy$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)



## Cancer Stage (Cancer Stage)

### Kaplan Meier Survival Curves for Cancer Stage
# Model the Kaplan Meier Survival Curve
surv_obj_stage <- with(brca_clin_stage, Surv(years_to_event, vital_status))
surv_fit_stage<- survfit(surv_obj_stage ~ pathologic_stage,data=brca_clin_stage)#times,
surv_fit_stage

#Plot the survival curves for the different groups
ggplot_stage_surv <- ggsurvplot(surv_fit_stage,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves by Cancer Stage",
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

#Plot the cumulative hazard rates for the different groups
ggplot_stage_cumhaz <- ggsurvplot(surv_fit_stage,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves by Cancer Stage",
           risk.table=TRUE,
           cumevents=TRUE,
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_stage_surv, ggplot_stage_cumhaz), ncol=1)

#Check for Difference between the survival curves
survdiff(surv_obj_stage ~ factor(pathologic_stage),data=brca_clin_stage)

brca_clin_stage_2 <- brca_clin_stage %>%
  filter(!is.na(pathologic_stage)) %>%
  mutate(pathologic_stage = factor(pathologic_stage))

surv_pairdiff_stage <- pairwise_survdiff(Surv(years_to_event,vital_status)~pathologic_stage,data=brca_clin)
surv_pairdiff_stage
symnum(surv_pairdiff_stage$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
   symbols = c("****", "***", "**", "*", "+", " "),
   abbr.colnames = FALSE, na = "")




### Cox Proportional Hazard (PH) Model for Pathologic Stage (Cancer Stage)
cox_stage <- coxph(surv_obj_stage ~ pathologic_stage, data=brca_clin_stage)

**Check for Proportional Hazards Assumption**
#check for proportional hazard assumptions
test_ph_stage <- cox.zph(cox_stage)
test_ph_stage

ggcoxzph(test_ph_stage)
#plot(test_ph_stage)
#ggcoxdiagnostics(cox_stage,type="schoenfeld", data=brca_clin_stage)


When tested for proportional hazards, it seems that the hazard ratios of stage 3 and stage 4 with stage 1 are not proportional. We will attempt by splitting the survival data into multiple groups by time cutoff.

**Split Survival Data into Time Groups**
brca_clin_stage_2 <- survSplit(Surv(years_to_event,vital_status) ~ pathologic_stage,
                               data=brca_clin_stage,
                               cut=c(8),
                               episode="time_group")
cox_stage_2 <- coxph(Surv(tstart,years_to_event,vital_status) ~pathologic_stage:strata(time_group),
                     data=brca_clin_stage_2)

test_ph_stage_2 <- cox.zph(cox_stage_2)
#ggcoxzph(test_ph_stage_2,ncol=1)


### Parametric modeling - Cancer Stage
#Fit Weibull 
surv_wei_stage <- survreg(surv_obj_stage ~ pathologic_stage, data=brca_clin_stage, dist="weibull")
summary(surv_wei_stage)
surv_wei_stage_conv <- ConvertWeibull(surv_wei_stage)


#Fit Exponential
surv_exp_stage <- update(surv_wei_stage,dist="exponential")
summary(surv_exp_stage)
lambda <- exp(-1.0 * unname(surv_exp_stage$coefficients[1]))
paste("Constant Hazard Rate for Exponential Distribution with stage as predictor:", round(lambda,2), " events/year")



# Compare the weibull and exponential distributions
cmp_wei_exp_stage <- anova(surv_wei_stage, surv_exp_stage)
cmp_wei_exp_stage
p_val <- cmp_wei_exp_stage$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp_stage$"-2*LL"[1]
deviance_exp <- cmp_wei_exp_stage$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)

ifelse(better_dist == "weibull", print(surv_wei_stage_conv), summary(surv_exp_stage))


## Tumor Stage (Pathology T)

### Kaplan Meier Survival Curves for Tumor Stage
# Model the Kaplan Meier Survival Curve
surv_obj_Tstage <- with(brca_clin_Tstage, Surv(years_to_event, vital_status))
surv_fit_Tstage<- survfit(surv_obj_Tstage ~ pathologyTstage,data=brca_clin_Tstage)#times,
surv_fit_Tstage

#Plot the survival curves for the different groups
ggplot_Tstage_surv <- ggsurvplot(surv_fit_Tstage,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves by Tumor Stage",
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

#Plot the cumulative hazard rates for the different groups
ggplot_Tstage_cumhaz <- ggsurvplot(surv_fit_Tstage,
           pval = FALSE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves by Tumor Stage",
           risk.table=TRUE,
           cumevents=TRUE,
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_Tstage_surv, ggplot_Tstage_cumhaz), ncol=1)

#Check for Difference between the survival curves
survdiff(surv_obj_Tstage ~ pathologyTstage,data=brca_clin_Tstage)
surv_pair_Tstage <- pairwise_survdiff(Surv(years_to_event,vital_status)~pathologyTstage,data=brca_clin)
symnum(surv_pair_Tstage$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
   symbols = c("****", "***", "**", "*", "+", " "),
   abbr.colnames = FALSE, na = "")



### Cox Proportional Hazard (PH) Model for Pathologic Stage (T Stage)
cox_Tstage <- coxph(surv_obj_Tstage ~ pathologyTstage, data=brca_clin_Tstage)

**Check for Proportional Hazards Assumption**
#check for proportional hazard assumptions
test_ph_Tstage <- cox.zph(cox_Tstage)
test_ph_Tstage

ggcoxzph(test_ph_Tstage)
#plot(test_ph_Tstage)
##ggcoxdiagnostics(cox_Tstage,type="schoenfeld", data=brca_clin_Tstage)



#**Split the Survival Data into time groups**
brca_clin_Tstage_2 <- survSplit(Surv(years_to_event,vital_status) ~ pathologyTstage,
                               data=brca_clin_Tstage,
                               cut=c(2),
                               episode="time_group")
cox_Tstage_2 <- coxph(Surv(tstart,years_to_event,vital_status) ~pathologyTstage:strata(time_group),
                     data=brca_clin_Tstage_2)

test_ph_Tstage_2 <- cox.zph(cox_Tstage_2)
#ggcoxzph(test_ph_Tstage_2)



### Parametric modeling - Tstage
#Fit Weibull 
surv_wei_Tstage <- survreg(surv_obj_Tstage ~ pathologyTstage, data=brca_clin_Tstage, dist="weibull")
surv_wei_Tstage_conv <- ConvertWeibull(surv_wei_Tstage)
summary(surv_wei_Tstage)

#Fit Exponential
surv_exp_Tstage <- update(surv_wei_Tstage,dist="exponential")
summary(surv_exp_Tstage)

lambda <- exp(-1.0 * unname(surv_exp_Tstage$coefficients[1]))
paste("Constant Hazard Rate for Exponential Distribution with Tstage as predictor:",round(lambda,2), " events/year")



# Compare the weibull and exponential distributions
cmp_wei_exp_Tstage <- anova(surv_wei_Tstage, surv_exp_Tstage)
p_val <- cmp_wei_exp_Tstage$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp_Tstage$"-2*LL"[1]
deviance_exp <- cmp_wei_exp_Tstage$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)

ifelse(better_dist == "weibull", print(surv_wei_Tstage_conv), summary(surv_exp_Tstage))




## Metastasis Status (Pathology M)

### Kaplan Meier Survival Curves for Metastasis Stage
# Model the Kaplan Meier Survival Curve
surv_obj_Mstage <- with(brca_clin_Mstage, Surv(years_to_event, vital_status))
surv_fit_Mstage<- survfit(surv_obj_Mstage ~ pathologyMstage,data=brca_clin_Mstage)#times,
surv_fit_Mstage

#Plot the survival curves for the different groups
ggplot_Mstage_surv <- ggsurvplot(surv_fit_Mstage,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves by Metastasis Stage",
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

#Plot the cumulative hazard rates for the different groups
ggplot_Mstage_cumhaz <- ggsurvplot(surv_fit_Mstage,
           pval = FALSE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves by Metastasis Stage",
           risk.table=TRUE,
           cumevents=TRUE,
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_Mstage_surv, ggplot_Mstage_cumhaz), ncol=1)

#Check for Difference between the survival curves
survdiff(surv_obj_Mstage ~ factor(pathologyMstage),data=brca_clin_Mstage)
surv_pair_Mstage <- pairwise_survdiff(Surv(years_to_event,vital_status)~pathologyMstage,data=brca_clin)
surv_pair_Mstage
symnum(surv_pair_Mstage$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
   symbols = c("****", "***", "**", "*", "+", " "),
   abbr.colnames = FALSE, na = "")


### Cox Proportional Hazard (PH) Model for Pathologic Stage (M Stage)
cox_Mstage <- coxph(surv_obj_Mstage ~ pathologyMstage, data=brca_clin_Mstage)


### Parametric modeling - Mstage
#Fit Weibull 
surv_wei_Mstage <- survreg(surv_obj_Mstage ~ pathologyMstage, data=brca_clin_Mstage, dist="weibull")
surv_wei_Mstage_conv <- ConvertWeibull(surv_wei_Mstage)
summary(surv_wei_Mstage)

#Fit Exponential
surv_exp_Mstage <- update(surv_wei_Mstage,dist="exponential")
summary(surv_exp_Mstage)
lambda <- exp(-1.0 * unname(surv_exp_Mstage$coefficients[1]))
paste("Constant Hazard Rate for Exponential Distribution with Mstage as predictor:",round(lambda,2), " events/year")


# Compare the weibull and exponential distributions
cmp_wei_exp_Mstage <- anova(surv_wei_Mstage, surv_exp_Mstage)
p_val <- cmp_wei_exp_Mstage$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp_Mstage$"-2*LL"[1]
deviance_exp <- cmp_wei_exp_Mstage$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)

ifelse(better_dist == "weibull", surv_wei_Mstage_conv, summary(surv_exp_Mstage))




## Lymph Nodes Status (Pathology N)

### Kaplan Meier Survival Curves for Lymph Nodes Stage
# Model the Kaplan Meier Survival Curve
surv_obj_Nstage <- with(brca_clin_Nstage, Surv(years_to_event, vital_status))
surv_fit_Nstage<- survfit(surv_obj_Nstage ~ pathologyNstage,data=brca_clin_Nstage)#times,
surv_fit_Nstage

#Plot the survival curves for the different groups
ggplot_Nstage_surv <- ggsurvplot(surv_fit_Nstage,
           pval = TRUE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           surv.median.line = "hv",
           legend="right",
           title="Survival Curves by Lymph Nodes Stage",
           fontsize=2,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

#Plot the cumulative hazard rates for the different groups
ggplot_Nstage_cumhaz <- ggsurvplot(surv_fit_Nstage,
           pval = FALSE,
           conf.int = TRUE,
           conf.int.fill="strata",
           conf.int.alpha=0.2,
           fun="cumhaz",
           legend="right",
           title="Cumulative Hazard Curves by Lymph Nodes Stage",
           fontsize=2,
           risk.table=TRUE,
           cumevents = TRUE,
           tables.height=0.3,
           tables.col="strata",
           surv.plot.height=1)

arrange_ggsurvplots(list(ggplot_Nstage_surv, ggplot_Nstage_cumhaz), ncol=1)

#Check for Difference between the survival curves
survdiff(surv_obj_Nstage ~ factor(pathologyNstage),data=brca_clin_Nstage)
surv_pair_Nstage <- pairwise_survdiff(Surv(years_to_event,vital_status)~pathologyNstage,data=brca_clin)
surv_pair_Nstage
symnum(surv_pair_Nstage$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
   symbols = c("****", "***", "**", "*", "+", " "),
   abbr.colnames = FALSE, na = "")


### Cox Proportional Hazard (PH) Model for Pathologic Stage (N Stage)
cox_Nstage <- coxph(surv_obj_Nstage ~ pathologyNstage, data=brca_clin_Nstage)

**Check for Proportional Hazards Assumption**
#check for proportional hazard assumptions
test_ph_Nstage <- cox.zph(cox_Nstage)
test_ph_Nstage

ggcoxzph(test_ph_Nstage)
#plot(test_ph_Nstage)


#ggcoxdiagnostics(cox_Nstage,type="schoenfeld", data=brca_clin_Nstage)

#**Interpretation of the Model**
summary(cox_Nstage)


ggforest(cox_Nstage,data=brca_clin_Nstage)
ggadjustedcurves(cox_Nstage, var="pathologyNstage", data = brca_clin_Nstage)

### Parametric modeling - Nstage
#Fit Weibull 
surv_wei_Nstage <- survreg(surv_obj_Nstage ~ pathologyNstage, data=brca_clin_Nstage, dist="weibull")
surv_wei_Nstage_conv <- ConvertWeibull(surv_wei_Nstage)
summary(surv_wei_Nstage)

#Fit Exponential
surv_exp_Nstage <- update(surv_wei_Nstage,dist="exponential")
summary(surv_exp_Nstage)
lambda <- exp(-1.0 * unname(surv_exp_Nstage$coefficients[1]))
paste("Constant Hazard Rate for Exponential Distribution with Nstage as predictor:",round(lambda,2), " events/year")


# Compare the weibull and exponential distributions
cmp_wei_exp_Nstage <- anova(surv_wei_Nstage, surv_exp_Nstage)
cmp_wei_exp_Nstage
p_val <- cmp_wei_exp_Nstage$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp_Nstage$"-2*LL"[1]
deviance_exp <- cmp_wei_exp_Nstage$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)

ifelse(better_dist == "weibull", print(surv_wei_Nstage_conv), summary(surv_exp_Nstage))




## Cox Proportional Hazard Model for All - no interactions

### Cox model with Age, Therapy and Cancer Stages

surv_obj_all <- with(brca_clin, Surv(years_to_event, vital_status))
cox_all <- coxph(surv_obj_all ~ age  + therapy_type + pathologic_stage , data=brca_clin)
cox_all_cat <- coxph(surv_obj_all ~ agecat  + therapy_type + pathologic_stage , data=brca_clin)


#**Check for Proportional Hazards Assumption**
# Test for PH significance
test_ph_all <- cox.zph(cox_all)
test_ph_all
plot(test_ph_all)



summary(cox_all)
plot(cox_all$residuals)



### Cox model with Age, Therapy, Tumor, and Lymph Node Stages
cox_all_2 <- coxph(surv_obj_all ~ age + therapy_type +  pathologyTstage + pathologyNstage, data=brca_clin)

**Check for Proportional Hazards Assumption**
test_ph_all_2 <- cox.zph(cox_all_2)
test_ph_all_2

#ggcoxzph(test_ph_all_2)
plot(test_ph_all_2)


**Interpretation of the Model**
```{r}
summary(cox_all_2)
plot(cox_all_2$residuals)


### Parametric modeling - All
#Fit Weibull 
surv_wei_all <- survreg(surv_obj_all ~ age + therapy_type + pathologic_stage, data=brca_clin, dist="weibull")
surv_wei_all_conv <- ConvertWeibull(surv_wei_all)
summary(surv_wei_all)

#Fit Exponential
surv_exp_all <- update(surv_wei_all,dist="exponential")
summary(surv_exp_all)
lambda <- exp(-1.0 * unname(surv_exp_all$coefficients[1]))
paste("Constant Hazard Rate for Exponential Distribution with all as predictor:",round(lambda,2), " events/year")

# Compare the weibull and exponential distributions
cmp_wei_exp_all <- anova(surv_wei_all, surv_exp_all)
p_val <- cmp_wei_exp_all$"Pr(>Chi)"[2]
deviance_wei <- cmp_wei_exp_all$"-2*LL"[1]
deviance_exp <- cmp_wei_exp_all$"-2*LL"[2]
better_dist <- ifelse(p_val < 0.05, ifelse(deviance_wei < deviance_exp, "weibull", "exp"), "both same")
paste("The better fit parametric model:",better_dist)

ifelse(better_dist == "weibull", surv_wei_all_conv, summary(surv_exp_all))


