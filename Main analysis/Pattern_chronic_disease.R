##Associations between ML_score and other health-related outcomes
library(survival)
ML_score <- fread("ML_score.csv")
covariate <- fread("covariate_mice.csv")
data <- merge(ML_score,covariate,by = "eid",sort = F)
output <- data.frame(exposure=character(),outcome=character(),model=as.character(),class=numeric(),class_value=factor(),
                     total=integer(),case=integer(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric())
output2 <- data.frame(exposure=character(),outcome=character(),model=as.character(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric())
start_date <- read.csv("UKB_dietary_date_max.csv",header = T)
start_date$start_date <- as.Date(start_date$start_date)
death <- read.csv("UKB_death.csv",header = T)
death$death <- pmax(death$X40000.0.0,death$X40000.1.0,na.rm = T)
end_date <- mutate(death,end_date=as.Date(ifelse(is.na(death),"2022-12-31",death))) %>% dplyr::select("eid","end_date")
outco <- merge(start_date,end_date,by="eid")
name <- fread("Target_code.csv")
for (out_name in name$Disease_code) {
  outcome <- read.csv(paste0("chronic/Targets_RAW/",out_name,".csv"),header = T)[,c("eid","target_date","target_y")]
  outcome <- merge(outco,outcome,by="eid")
  outcome <- dplyr::select(outcome,eid,start_date,end_date,target_date) %>% mutate(target_date=as.Date(target_date)) %>%
    mutate(days=ifelse(is.na(.data$target_date),difftime(end_date,start_date,units = "days"),difftime(.data$target_date,start_date,units = "days")),
           status=ifelse(is.na(.data$target_date),0,1),.keep="unused",.after=eid)
  data1 <- merge(outcome,data,by = "eid",sort = F)
  data1 <- filter(data1,days > 0)
  exposure <- "ML_score"
  out <- "status"
  func_quan <- function(exposure,out,model,input){
    data <- select(input,days,all_of(out),all_of(exposure),energy:BMI)
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
                         outcome=rep(name$Disease[name$Disease_code == out_name],a),
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
  for (i in exposure) {
    for (j in out) {
      for (k in model) {
        output0 <- func_quan(i,j,k,dat)
        output <- rbind(output,output0)
      }
    }
  }
  
  func_linear <- function(exposure,out,model,input){
    if(model == 1){data <- select(input,days,all_of(out),all_of(exposure),energy:BMI)
    }else if(model == 2){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD)}
    setnames(data,c(out,exposure),c("outcome","pheno"))
    data$pheno <- scale(data$pheno,scale = 0.2*7)
    n <- nrow(data);case <- sum(data$outcome,na.rm = T)
    fit <- coxph(Surv(days,outcome) ~ .,data = data)
    fit0 <- summary(fit)
    coef <- as.data.frame(fit0$coefficients)
    confint <- as.data.frame(fit0$conf.int)
    output <- data.frame(exposure=exposure,
                         outcome=name$Disease[name$Disease_code == out_name],
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
  for (i in exposure) {
    for (j in out) {
      for (k in model) {
        output0 <- func_linear(i,j,k,dat)
        output <- rbind(output,output0)
      }
    }
  }
}
output$FDR <- FDR(output$P,5)
write.csv(output,paste0("pattern_result_linear_New","",".csv"),row.names = F)
write.csv(output2,paste0("pattern_result_New","",".csv"),row.names = F)