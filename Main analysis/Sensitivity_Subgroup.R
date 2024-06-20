##Sensitivity analysis
library(survival)
ML_score <- fread("ML_score.csv")
covariate <- fread("covariate_mice.csv")
outcome <- fread("Dementia.csv")
data0 <- merge(ML_score,covariate,by = "eid",sort = F)
data0 <- merge(outcome,data0,by = "eid",sort = F)

#Exclude any component from ML_score
ML_score <- select(ML_score,-`potatoes`)
# ML_score <- select(ML_score,-`green leafy vegetables`)
# ML_score <- select(ML_score,-`olive oil`)
# ML_score <- select(ML_score,-`berries and citrus fruits`)
# ML_score <- select(ML_score,-`eggs`)
# ML_score <- select(ML_score,-`poultry`)
# ML_score <- select(ML_score,-`sweetened beverages`)
ML_score$ML_score <- rowSums(ML_score[,2:7])
note <- "minus potatoes"
note <- "minus green leafy vegetables"
note <- "minus olive oil"
note <- "minus berries and citrus fruits"
note <- "minus eggs"
note <- "minus poultry"
note <- "minus sweetened beverages"
#Excluding events in the first 2 or 5 years
data <- filter(data0,days > 730)
note <- "2-year"
data <- filter(data0,days > 1826)
note <- "5-year"
#Model 2 + Fluid intelligence score at baseline
cognition <- fread("cognition_baseline.csv")
data <- merge(data0,cognition,by = "eid")
summary(select(data,FI:Pairs))
note <- "2006-2010_cognition"
#Model 2 + History of depression at baseline
depression <- fread("Depression.csv")
data <- merge(data0,depression[,c("eid","depression")],by = "eid")
note <- "depression"
#Model 2 + 24 other food groups
group <- fread("Food_group.csv")
in_group <- c("Green leafy vegetables","Olive oil","Berries","Citrus fruits","Potatoes","Eggs","Poultry","Sweetened beverages")
other_group <- select(group,-all_of(in_group))
data <- merge(data0,other_group,by = "eid")
note <- "other group"

##Subgroup analysis
#stratified by sex
data <- filter(data0,sex == 1)
note <- "Man"
data <- filter(data0,sex == 2)
note <- "Woman"
#stratified by age
data <- filter(data0,age < 65)
note <- "<65"
data <- filter(data0,age >= 65)
note <- ">=65"
#stratified by APOE4 allele
data <- filter(data0,APOE4 == 0)
note <- "non-carrier"
data <- filter(data0,APOE4 == 1)
note <- "Carrier"
#stratified by BMI
data <- filter(data0,BMI <= 25)#
note <- "normal weight"
data <- filter(data0,BMI > 25)
note <- "overweight"
#stratified by smoking status
data <- filter(data0,smoking == 0)
note <- "never"
data <- filter(data0,smoking == 1)
note <- "former"
data <- filter(data0,smoking == 2)
note <- "current"
#stratified by prs
prs <- fread("PRS_base_IGAP_stage_1_all_20211126.txt")
prs <- select(prs,IID,pT_0.5) %>% 
  mutate(PRS = cut(pT_0.5,breaks=quantile(pT_0.5,c(0,0.2,0.8,1)),include.lowest = T)) %>% 
  mutate(PRS = as.integer(PRS))
data0 <- merge(data0,prs,by.x = "eid",by.y = "IID")
data <- filter(data0,PRS == 1)
note <- "Low genetic risk"
data <- filter(data0,PRS == 2)
note <- "Intermediate genetic risk"
data <- filter(data0,PRS == 3)
note <- "High genetic risk"

#analysis
output <- data.frame(exposure=character(),outcome=character(),model=as.character(),class=numeric(),class_value=factor(),
                     total=integer(),case=integer(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric(),note=character())

func_quan <- function(exposure,out,model){
  if(model == 1){data <- select(data,days,all_of(out),all_of(exposure),energy:APOE4)
  }else if(model == 2){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD)
  }else if(model == 3){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD,FI)
  }else if(model == 4){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD,depression)
  }else if(model == 5){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD,`Other fruits`:`Fried and fast foods`)}
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
                       P=c(1,coef$`Pr(>|z|)`[1:b]),
                       note=rep(note,a))
  return(output)
}


func_linear <- function(exposure,out,model){
  if(model == 1){data <- select(data,days,all_of(out),all_of(exposure),energy:APOE4)
  }else if(model == 2){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD)
  }else if(model == 3){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD,FI)
  }else if(model == 4){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD,depression)
  }else if(model == 5){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD,`Other fruits`:`Fried and fast foods`)}
  setnames(data,c(out,exposure),c("outcome","pheno"))
  data$pheno <- scale(data$pheno,scale = 0.2*7)
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
                       P=c(coef$`Pr(>|z|)`[1]),
                       note=note)
  return(output)
}

exposure <- "ML_score"
out <- c("ACD_status","AD_status")
model <- 2:5#assignment
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
write.csv(output,"Result/ML_dementia_sense_subgroup.csv",row.names = F)