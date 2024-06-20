##Associations between food items and incident dementia
library(survival)
Food <- fread("Food_items.csv")
covariate <- fread("covariate_mice.csv")
outcome <- fread("Dementia.csv")[,-1]
data <- merge(data,covariate,by = "eid",sort = F)
data <- merge(data,outcome,by = "eid",sort = F)
#Divided into intake and non-intake groups
func_dich <- function(x){
  a <- ifelse(x == 0,0,1)
  return(a)
}
data1 <- mutate(data,across(.cols = 2:ncol(Food),func1))
#Multivariate Cox proportional hazards model one
output <- data.frame(exposure=character(),class=numeric(),outcome=character(),model=character(),total=integer(),
                     case=integer(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric())
for (p in colnames(group)[-1]) {
  data4 <- select(data3,eid,days,ACD_status,as.character(p),energy:APOE4)
  colnames(data4)[4] <- "pheno"
  data4$pheno <- scale(data4$pheno)
  data4 <- na.omit(data4)
  total <- group_by(data4,pheno) %>% summarise(n=n())
  case <- group_by(data4,pheno) %>% summarise(n=sum(ACD_status))
  fit <- coxph(Surv(days,ACD_status) ~ .,data = data4[,-1])
  fit0 <- summary(fit)
  res <- cox.zph(fit)
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  output0 <- data.frame(exposure=rep(p,2),
                        class=1:2,
                        outcome=rep("dementia",2),
                        total=total$n[1:2],
                        case=case$n[1:2],
                        HR=c(1,coef$`exp(coef)`[1]),
                        LCI=c(1,confint$`lower .95`[1]),
                        HCI=c(1,confint$`upper .95`[1]),
                        P=c(1,coef$`Pr(>|z|)`[1]),
                        res.zph = res$table[1,3])
  output <- rbind(output,output0)
}
FDR <- function(p,n){
  row <- length(p)
  for (i in 1:n) {
    p[seq(i,row,n)] <- p.adjust(p[seq(i,row,n)],method = "fdr")
  }
  return(p)
}
output$FDR <- FDR(output$P,2)
write.csv(output,"Result/Food_dementia_model1.csv",row.names = F)
#Multivariate Cox proportional hazards model two
output <- data.frame(exposure=character(),class=numeric(),outcome=character(),model=character(),total=integer(),
                     case=integer(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric())
for (p in colnames(group)[-1]) {
  data4 <- select(data3,eid,days,ACD_status,as.character(p),energy:CVD)
  colnames(data4)[4] <- "pheno"
  data4$pheno <- scale(data4$pheno)
  data4 <- na.omit(data4)
  total <- group_by(data4,pheno) %>% summarise(n=n())
  case <- group_by(data4,pheno) %>% summarise(n=sum(ACD_status))
  fit <- coxph(Surv(days,ACD_status) ~ .,data = data4[,-1])
  fit0 <- summary(fit)
  res <- cox.zph(fit)
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  output0 <- data.frame(exposure=rep(p,2),
                        class=1:2,
                        outcome=rep("dementia",2),
                        total=total$n[1:2],
                        case=case$n[1:2],
                        HR=c(1,coef$`exp(coef)`[1]),
                        LCI=c(1,confint$`lower .95`[1]),
                        HCI=c(1,confint$`upper .95`[1]),
                        P=c(1,coef$`Pr(>|z|)`[1]),
                        res.zph = res$table[1,3])
  output <- rbind(output,output0)
}
FDR <- function(p,n){
  row <- length(p)
  for (i in 1:n) {
    p[seq(i,row,n)] <- p.adjust(p[seq(i,row,n)],method = "fdr")
  }
  return(p)
}
output$FDR <- FDR(output$P,2)
write.csv(output,"Result/Food_dementia_model2.csv",row.names = F)