##Associations between dietary patterns and incident dementia
library(survival)
ML_score <- fread("ML_score.csv")
MIND_score <- fread("MIND_score.csv")
covariate <- fread("covariate_mice.csv")
outcome <- fread("Dementia.csv")
data <- merge(ML_score,MIND_score,by = "eid",sort = F)
data <- merge(data,covariate,by = "eid",sort = F)
data <- merge(outcome,data,by = "eid",sort = F)
#Multivariate Cox proportional hazards model one
func_quan <- function(exposure,out,model){
  if(model == 1){data <- select(data,days,all_of(out),all_of(exposure),energy:APOE4)
  }else if(model == 2){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD)}
  setnames(data,c(out,exposure),c("outcome","pheno"))
  data$pheno <- cut(data$pheno,breaks = quantile(data$pheno),include.lowest = T)
  total <- group_by(data,pheno) %>% summarise(n=n())
  case <- group_by(data,pheno) %>% summarise(n=sum(outcome))
  fit <- coxph(Surv(days,outcome) ~ .,data = data)
  fit0 <- summary(fit)
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  n=2;a=n+2;b=n+1
  output <- data.frame(exposure=rep(exposure,a),
                       outcome=rep(out,a),
                       model=rep(paste0("model_",model),a),
                       class=c(1:a),
                       class_value=total$pheno[1:a],
                       total=total$n[1:a],
                       case=case$n[1:a],
                       HR=c(1,coef$`exp(coef)`[1:b]),
                       LCI=c(1,confint$`lower .95`[1:b]),
                       HCI=c(1,confint$`upper .95`[1:b]),
                       P=c(1,coef$`Pr(>|z|)`[1:b]))
  return(output)
}
func_linear <- function(exposure,out,model){
  if(model == 1){data <- select(data,days,all_of(out),all_of(exposure),energy:APOE4)
  }else if(model == 2){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD)}
  setnames(data,c(out,exposure),c("outcome","pheno"))
  # data$pheno <- scale(data$pheno,scale = 0.2*7)
  data$pheno <- (data$pheno-mean(data$pheno))/1.4
  n <- nrow(data);case <- sum(data$outcome,na.rm = T)
  fit <- coxph(Surv(days,outcome) ~ .,data = data)
  fit0 <- summary(fit)
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  output <- data.frame(exposure=exposure,
                       outcome=out,
                       model=paste0("model_",model),
                       class=5,
                       class_value="Per 20% increment",
                       total=n,
                       case=case,
                       HR=c(coef$`exp(coef)`[1]),
                       LCI=c(confint$`lower .95`[1]),
                       HCI=c(confint$`upper .95`[1]),
                       P=c(coef$`Pr(>|z|)`[1]))
  return(output)
}

exposure <- c("ML_score","MIND_score")
out <- c("ACD_status","AD_status")
model <- 1:2
output <- data.frame(exposure=character(),outcome=character(),model=as.character(),class=numeric(),class_value=factor(),
                     total=integer(),case=integer(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric())
for (i in exposure) {
  for (j in out) {
    for (k in model) {
      output0 <- func_quan(i,j,k)
      output <- rbind(output,output0)
      output0 <- func_linear(i,j,k)
      output <- rbind(output,output0)
    }
  }
}
write.csv(output,"Result/Pattern_dementia.csv",row.names = F)