##Associations between ML_score and serum proteins
prot <- fread("ProteomicsData_raw_new.csv")
#quality control
func_miss <- function(x){
  a <- length(which(is.na(x) == T))/length(x)
  return(a)
}
b <- apply(prot, MARGIN = 2, func_miss)
prot1 <- dplyr::select(prot,all_of(which(as.numeric(b < 0.25)==T)))
ML_score <- fread("ML_score.csv")
covariate <- fread("covariate_mice.csv")
outcome <- fread("Dementia.csv")
data <- merge(ML_score,prot1,by="eid")
data <- merge(data,covariate,by = "eid")
metabolic_out <- merge(outcome,data,by = "eid")
#Linear regression
library(purrr)
reg_result <- data.frame(exposure=character(),outcome=character(),estimate = numeric(),p=numeric())
for (i in c(1:2919)){
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
}

reg_result$fdr <- p.adjust(reg_result$p,method = "fdr")
write.csv(reg_result,"Pattern_prot.csv",row.names = F)

##Associations between serum proteins and incident dementia
result <- data.frame(exposure=character(),outcome=character(),estimate = numeric(),p=numeric())
for (i in colnames(prot1)[-1]) {
  print(i)
  temp <- dplyr::select(data,ACD_status,all_of(i),energy:CVD)
  colnames(temp)[2] <- "pheno"
  fit <- glm(ACD_status~pheno + energy + age + sex + ethnic + townsend + 
               education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,family = binomial(),data = temp)
  fit0 <- summary(fit)
  pvalue <- fit0[["coefficients"]][2,4]
  estimate <- fit0[["coefficients"]][2,1]
  result0 <- data.frame(exposure = i,outcome = "dementia",estimate = estimate,p=pvalue)
  result <- rbind(result,result0)
}
result$FDR <- FDR(result$p,1)
library(mediation)
mediation <- function(pheno,mediator,outcome){
  data <- dplyr::select(data0,all_of(c(pheno,mediator,outcome)),energy:CVD)
  colnames(data)[1:3] <- c("pheno","mediator","outcome")
  fit_a <- lm(mediator~pheno + energy + age + sex + ethnic + townsend + 
                education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,data)##线性
  fit_a_summary <- summary(fit_a)
  coef_a <- as.data.frame(fit_a_summary$coefficients) %>% slice(2:2) %>%
    mutate(pheno=pheno,mediator=mediator,.before="Estimate")
  colnames(coef_a)[3:6] <- c("a.estimate","a.std.error","a.t.value","a.p.value")
  fit_b <- glm(outcome~ mediator + energy + age + sex + ethnic + townsend + 
                 education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,family = binomial(),data)##线性
  fit_b_summary <- summary(fit_b)
  coef_b <- as.data.frame(fit_b_summary$coefficients) %>% slice(2:2) %>%
    mutate(mediator=mediator,outcome=outcome,.before="Estimate")
  colnames(coef_b)[3:6] <- c("b.value","b.std.error","b.z.value","b.p.value")
  fit_c <- glm(outcome~pheno + energy + age + sex + ethnic + townsend + 
                 education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,family = binomial(),data)##线性
  fit_c_summary <- summary(fit_c)
  coef_c <- as.data.frame(fit_c_summary$coefficients) %>% slice(2:2) %>%
    mutate(mediator=mediator,outcome=outcome,.before="Estimate")
  colnames(coef_c)[3:6] <- c("c.value","c.std.error","c.z.value","c.p.value")
  fit_c_adj <- glm(outcome~pheno + mediator + energy + age + sex + ethnic + townsend + 
                     education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,family = binomial(),data)##线性
  contcont <- mediate(fit_a, fit_c_adj,sims=100,treat = 'pheno',mediator ='mediator')
  contsuma <- summary(contcont)
  effect <- c("TE","ACME","ADE","prop")
  estimate <- c(contsuma$tau.coef,contsuma$d.avg,contsuma$z.avg,contsuma$n.avg)
  LCI <- c(contsuma$tau.ci[1],contsuma$d.avg.ci[1],contsuma$z.avg.ci[1],contsuma$n.avg.ci[1])
  HCI <- c(contsuma$tau.ci[2],contsuma$d.avg.ci[2],contsuma$z.avg.ci[2],contsuma$n.avg.ci[2])
  p.value <- c(contsuma$tau.p,contsuma$d.avg.p,contsuma$z.avg.p,contsuma$n.avg.p)
  output <- data.frame(group="low",pheno=pheno,mediator=mediator,disease=outcome,effect=effect,
                       estimate=estimate,LCI=LCI,HCI=HCI,p.value=p.value,a.estimate=coef_a$a.estimate,
                       a.p.value=coef_a$a.p.value,b.estimate=coef_b$b.value,b.p.value=coef_b$b.p.value,c.estimate=coef_c$c.value,c.p.value=coef_c$c.p.value)
  return(output)
}
#analysis
prot_ML <- filter(Pattern_prot,FDR<0.05)
prot_dementia <- filter(Prot_dementia,FDR<0.05)
prot_selec <- filter(prot_ML,outcome %in% prot_dementia$exposure)
write.table(prot_selec$outcome,"PPI_protein.txt",row.names = F)
output0 <- data.frame()
for (i in prot_selec$outcome) {
  output <- mediation("ML_score",i,"ACD_status")
  output0 <- rbind(output0,output)
}
mediator <- filter(output0,effect == "prop",p.value < 0.05,estimate>0)
write.csv(output0,"Result/Prot_mediation.csv",row.names = F)
write.table(mediator$mediator,"PPI_protein.txt",row.names = F)