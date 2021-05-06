
library(dplyr)
library(survival)
library(survminer)

setwd("~/github/FLOMICS/")

clusters <- read.csv("Cluster Labels/InfiniumClust_SNF_tSeq_Labels_10Feb2021.csv")
survival <- read.csv("metadata/clinical_data_upd7Apr2021.csv")[,c(1:27)] %>% mutate(TTP = as.numeric(TTP))

survival.df <- clusters %>%
  filter(!is.na(SNFClust10Feb2021)) %>%
  mutate(ID = substr(ID, 1, 9)) %>%
  left_join(survival, by = c("ID" = "LY_FL_ID")) %>%
  mutate(SNFClust10Feb2021 = factor(SNFClust10Feb2021, levels = c("2", "1"))) %>%
  mutate(FLIPI_BINARY = factor(FLIPI_BINARY, levels = c("LOW_INTERMEDIATE", "HIGH"))) %>%
  mutate(TYPE = factor(TYPE, levels = c("LIMITED", "ADVANCED")))

# exclude <- c("LY_FL_046_T1", "LY_FL_316_T1", "LY_FL_018_T1", "LY_FL_248_T1", "LY_FL_181_T1", "LY_FL_318_T1", "LY_FL_119_T1",
#              "LY_FL_276_T1", "LY_FL_273_T1", "LY_FL_136_T1", "LY_FL_445_T1", "LY_FL_196_T1", "LY_FL_179_T1", "LY_FL_167_T1")
# exclude <- substr(exclude, 1, 9)
# 
# survival.df <- survival.df %>%
#   filter(!ID %in% exclude)

# median follow-up of living patients
survival.df %>%
  filter(CODE_OS == "0") %>%
  .$OS %>%
  median()

# Univariate analysis  
fit <- survfit(Surv(TTP, CODE_TTP) ~ SNFClust10Feb2021, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)
summary(fit, times = c(5, 7))

fit <- survfit(Surv(TTP, CODE_PFS) ~ SNFClust10Feb2021, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

fit <- survfit(Surv(OS, CODE_OS) ~ SNFClust10Feb2021, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

fit <- survfit(Surv(OS, CODE_DSS) ~ SNFClust10Feb2021, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)
summary(fit, times = c(5, 7))

fit <- survfit(Surv(TTT, CODE_TRANSF) ~ SNFClust10Feb2021, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE,
           fun = "cumhaz", ylim = c(0,1))
summary(fit, times = c(5, 7))

# Cox multivariate analysis for DSS
surv_object <- Surv(time = survival.df$OS, event = survival.df$CODE_DSS)
fit.coxph <- coxph(surv_object ~ SNFClust10Feb2021 + FLIPI_BINARY, data = survival.df)
cox.zph(fit.coxph)
ggforest(fit.coxph, data = survival.df)

surv_object <- Surv(time = survival.df$OS, event = survival.df$CODE_DSS)
fit.coxph <- coxph(surv_object ~ SNFClust10Feb2021 + TYPE, data = survival.df)
cox.zph(fit.coxph)
ggforest(fit.coxph, data = survival.df)

surv_object <- Surv(time = survival.df$OS, event = survival.df$CODE_DSS)
fit.coxph <- coxph(surv_object ~ SNFClust10Feb2021 + FLIPI_BINARY + TYPE, data = survival.df)
cox.zph(fit.coxph)
ggforest(fit.coxph, data = survival.df)

# Cox multivariate analysis for TTT
surv_object <- Surv(time = survival.df$TTT, event = survival.df$CODE_TRANSF)
fit.coxph <- coxph(surv_object ~ SNFClust10Feb2021 + FLIPI_BINARY, data = survival.df)
cox.zph(fit.coxph)
ggforest(fit.coxph, data = survival.df)

surv_object <- Surv(time = survival.df$TTT, event = survival.df$CODE_TRANSF)
fit.coxph <- coxph(surv_object ~ SNFClust10Feb2021 + TYPE, data = survival.df)
cox.zph(fit.coxph)
ggforest(fit.coxph, data = survival.df)

surv_object <- Surv(time = survival.df$TTT, event = survival.df$CODE_TRANSF)
fit.coxph <- coxph(surv_object ~ SNFClust10Feb2021 + FLIPI_BINARY + TYPE, data = survival.df)
cox.zph(fit.coxph)
ggforest(fit.coxph, data = survival.df)

# Plotting KM with 4 groups (SNF cluster / TYPE)

fit <- survfit(Surv(TTP, CODE_TTP) ~ SNFClust10Feb2021 + TYPE, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

fit <- survfit(Surv(TTP, CODE_PFS) ~ SNFClust10Feb2021 + TYPE, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

fit <- survfit(Surv(OS, CODE_OS) ~ SNFClust10Feb2021 + TYPE, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

fit <- survfit(Surv(OS, CODE_DSS) ~ SNFClust10Feb2021 + TYPE, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

# Plotting KM with 4 groups (SNF cluster / FLIPI)

fit <- survfit(Surv(TTP, CODE_TTP) ~ SNFClust10Feb2021 + FLIPI_BINARY, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

fit <- survfit(Surv(TTP, CODE_PFS) ~ SNFClust10Feb2021 + FLIPI_BINARY, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

fit <- survfit(Surv(OS, CODE_OS) ~ SNFClust10Feb2021 + FLIPI_BINARY, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)

fit <- survfit(Surv(OS, CODE_DSS) ~ SNFClust10Feb2021 + FLIPI_BINARY, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE)
pairwise_survdiff(Surv(OS, CODE_DSS) ~ SNFClust10Feb2021 + FLIPI_BINARY, data = survival.df)

fit <- survfit(Surv(TTT, CODE_TRANSF) ~ SNFClust10Feb2021 + FLIPI_BINARY, data = survival.df)
ggsurvplot(fit, data = survival.df, risk.table = TRUE, pval = TRUE,
           fun = "cumhaz", ylim = c(0,1))


