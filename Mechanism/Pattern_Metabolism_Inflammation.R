##Associations between ML_score and metabolites
metabolites <- fread("Metabolites.csv")
ML_score <- fread("ML_score.csv")
covariate <- fread("covariate_mice.csv")
outcome <- fread("Dementia.csv")
data <- merge(ML_score,metabolites,by="eid")
data <- merge(data,covariate,by = "eid")
data <- merge(outcome,data,by = "eid")
##Linear regression
library(purrr)
reg_result <- data.frame(exposure=character(),outcome=character(),estimate = numeric(),p=numeric())
for (i in c(1:249)){
  exposure <- as.data.frame(scale(data[,c(13)]))
  outcome_name <- colnames(data)[i+13]
  outcome <- scale(data[[i+13]])
  reg <- map(exposure,~ lm(outcome~.x + energy + age + sex + ethnic + townsend + 
                             education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,data))
  pvalue <- map_df(reg,~summary(.x)[["coefficients"]][2,4])
  estimate <- map_df(reg,~summary(.x)[["coefficients"]][2,1])
  reg_result0 <- rbind(estimate,pvalue)
  reg_result0 <- as.data.frame(t(reg_result0))
  colnames(reg_result0) <- c("estimate","p")
  reg_result0 <- mutate(reg_result0,exposure = rownames(reg_result0),outcome = outcome_name,
                        .before = estimate)
  reg_result <- rbind(reg_result,reg_result0)
  print(i)
}

reg_result$fdr <- p.adjust(reg_result$p,method = "fdr")
write.csv(reg_result,"Pattern_metabolite.csv",row.names = F)

##Associations between ML_score and inflammation indexes
inflammation <- fread("Inflammation.csv")
inflammation <- filter(inflammation,!if_all(.cols = 2:10,~ is.na(.) == T))
ML_score <- fread("ML_score.csv")
covariate <- fread("covariate_mice.csv")
outcome <- fread("Dementia.csv")
data <- merge(ML_score,inflammation,by="eid")
data <- merge(data,covariate,by = "eid")
data <- merge(outcome,data,by = "eid")
#Linear regression
reg_result <- data.frame(exposure=character(),outcome=character(),estimate = numeric(),p=numeric())
for (i in c(1:9)){
  exposure <- as.data.frame(scale(data[,c(13)]))
  outcome_name <- colnames(data)[i+13]
  outcome <- scale(data[[i+13]])
  reg <- map(exposure,~ lm(outcome~.x + energy + age + sex + ethnic + townsend + 
                             education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,data))
  pvalue <- map_df(reg,~summary(.x)[["coefficients"]][2,4])
  estimate <- map_df(reg,~summary(.x)[["coefficients"]][2,1])
  reg_result0 <- rbind(estimate,pvalue)
  reg_result0 <- as.data.frame(t(reg_result0))
  colnames(reg_result0) <- c("estimate","p")
  reg_result0 <- mutate(reg_result0,exposure = rownames(reg_result0),outcome = outcome_name,
                        .before = estimate)
  reg_result <- rbind(reg_result,reg_result0)
  print(i)
}

reg_result$fdr <- p.adjust(reg_result$p,method = "fdr")
write.csv(reg_result,"Pattern_inflammation.csv",row.names = F)

##SEM
library(lavaan)
library(semPlot)
data <- merge(ML_score,inflammation,by="eid")
data <- merge(data,metabolites,by = "eid")
data <- merge(outcome,data,by = "eid")
#Variable selection
sig_infla <- read.xlsx("Infla_sem.xlsx",sheet = "Inflammation")
sig_metab <- read.xlsx("Metab_sem.xlsx",sheet = "Metabolite")
infla_New <- sig_infla %>% select(1:5) %>% filter(fdr < 0.05) %>% arrange(fdr)
metab_New <- sig_metab %>% select(1:5) %>% filter(rank(fdr) < 11) %>% arrange(fdr)
data_sem <- select(data,all_of(infla_New$outcome),all_of(metab_New$outcome),ML_score,ACD_status)
data_sem <- mutate(data_sem,across(.cols = 1:18,scale))
colnames(data_sem)[8:17] <- paste0("Metabolite",1:10)
#analysis
model <- '# measurement model
          Inflammation =~ Neutrophils + SII + CRP + Platelets + NLR + Monocytes + Lymphocytes
          Metabolites =~ Metabolite1 + Metabolite2 + Metabolite3 + Metabolite4 + Metabolite5 + Metabolite6 + Metabolite7 + Metabolite8 + Metabolite9 + Metabolite10
          New =~ New_score
          ACD =~ ACD_status
          
          # structural model
          Inflammation ~ New + Metabolites
          Metabolites ~ New
          ACD ~ New + Inflammation + Metabolites'
fit <- sem(model, data = data_sem, bounds = T)
lavInspect(fit, 'optim.gradient')
result <- summary(fit, standardized = TRUE) 
write.table(result,"Pattern_infla_metab_sem.txt")