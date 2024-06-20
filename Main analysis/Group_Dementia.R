##Associations between Food groups and incident dementia
library(survival)
group <- fread("Food_group.csv")
covariate <- fread("covariate_imp.csv")
outcome <- fread("Dementia.csv")
data <- merge(group,covariate,by = "eid",sort = F)
data <- merge(data,outcome,by = "eid",sort = F)
#Divided by quintiles
func_quin <- function(x){
  breaks = quantile(x,probs = seq(0,1,1/5))
  if(length(unique(breaks)) == 2){a <- factor(ifelse(x>0,1,0))
  }else if(length(unique(breaks)) < 6){a <- cut(x,breaks = c(min(x)-1,unique(breaks)),include.lowest = T)
  }else{a <- cut(x,breaks = unique(breaks),include.lowest = T)}
  return(a)
}
data1 <- mutate(data,across(.cols = 2:ncol(group),func_quin))
#Statistics
stat <- function(x,data_1,data_2){
  temp <- data_1[,1:37]
  setnames(temp,x,"quantile")
  temp <- merge(temp,select(data_2,eid,all_of(x)),by = "eid")
  setnames(temp,x,"Food")
  b <- temp %>% group_by(quantile) %>% 
    summarise(Median = round(median(Food),2),Range = paste(round(range(Food),2),collapse = "~"),
              n = n(),case = sum(ACD_status)) %>%
    mutate(pheno = x,.before = "quantile") %>%
    mutate(`Consumption level` = paste0(Median," (",Range,")"),.after = "quantile")
  return(b)
}
stat_data <- data.frame()
for (i in colnames(group)[-1]) {
  stat0 <- stat(i,data1,data0)
  stat_data <- rbind(stat_data,stat0)
}
write.csv(stat_data,"Statistic/Group_stat.csv",row.names = F)
#Multivariate Cox proportional hazards model
data2 <- mutate(data1,across(.cols = 2:ncol(group),as.integer))
func_cox <- function(component,ref){
  data3 <- select(data2,eid,days,ACD_status,all_of(component),energy:APOE4)
  setnames(data3,component,"pheno")
  n = length(unique(data3$pheno))
  if(ref==1){data3$pheno <- factor(data3$pheno)
  }else if(ref==2){
    if(n==2){data3$pheno <- factor(data3$pheno,levels = c(2,1))
    }else if(n==3){data3$pheno <- factor(data3$pheno,levels = c(2,1,3))
    }else if(n==4){data3$pheno <- factor(data3$pheno,levels = c(2,1,3,4))
    }else if(n==5){data3$pheno <- factor(data3$pheno,levels = c(2,1,3,4,5))}
  }else if(ref==3){
    if(n==2){data3$pheno <- factor(data3$pheno,levels = c(2,1))
    }else if(n==3){data3$pheno <- factor(data3$pheno,levels = c(3,1,2))
    }else if(n==4){data3$pheno <- factor(data3$pheno,levels = c(3,1,2,4))
    }else if(n==5){data3$pheno <- factor(data3$pheno,levels = c(3,1,2,4,5))}
  }else if(ref==4){
    if(n==2){data3$pheno <- factor(data3$pheno,levels = c(2,1))
    }else if(n==3){data3$pheno <- factor(data3$pheno,levels = c(3,1,2))
    }else if(n==4){data3$pheno <- factor(data3$pheno,levels = c(4,1,2,3))
    }else if(n==5){data3$pheno <- factor(data3$pheno,levels = c(4,1,2,3,5))}
  }else if(ref==5){
    if(n==2){data3$pheno <- factor(data3$pheno,levels = c(2,1))
    }else if(n==3){data3$pheno <- factor(data3$pheno,levels = c(3,1,2))
    }else if(n==4){data3$pheno <- factor(data3$pheno,levels = c(4,1,2,3))
    }else if(n==5){data3$pheno <- factor(data3$pheno,levels = c(5,1,2,3,4))}
  }
  class(data3$pheno);levels(data3$pheno)
  fit <- coxph(Surv(days,ACD_status) ~ .,data = data3[,-1])
  res <- cox.zph(fit)
  fit0 <- summary(fit)
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  a = n;b=n-1
  output <- data.frame(Ref=rep(ref,a),
                       exposure=rep(component,a),
                       outcome=rep("Depression",a),
                       class=levels(data3$pheno),
                       HR=c(1,coef$`exp(coef)`[1:b]),
                       LCI=c(1,confint$`lower .95`[1:b]),
                       HCI=c(1,confint$`upper .95`[1:b]),
                       P=c(1,coef$`Pr(>|z|)`[1:b]),
                       res.zph = rep(res$table[1,3],a))
  return(output)
}
output <- data.frame()
for (ref in 1:5) {
  for (i in colnames(group)[2:ncol(group)]) {
    output0 <- func_cox(i,ref)
    output <- rbind(output,output0)
  }
}
write.csv(output,"Result/Group_dementia_model1.csv",row.names = F)
#Multivariate Cox proportional hazards model
data2 <- mutate(data1,across(.cols = 2:ncol(group),as.integer))
func_cox <- function(component,ref){
  data3 <- select(data2,eid,days,ACD_status,all_of(component),energy:CVD)
  setnames(data3,component,"pheno")
  n = length(unique(data3$pheno))
  if(ref==1){data3$pheno <- factor(data3$pheno)
  }else if(ref==2){
    if(n==2){data3$pheno <- factor(data3$pheno,levels = c(2,1))
    }else if(n==3){data3$pheno <- factor(data3$pheno,levels = c(2,1,3))
    }else if(n==4){data3$pheno <- factor(data3$pheno,levels = c(2,1,3,4))
    }else if(n==5){data3$pheno <- factor(data3$pheno,levels = c(2,1,3,4,5))}
  }else if(ref==3){
    if(n==2){data3$pheno <- factor(data3$pheno,levels = c(2,1))
    }else if(n==3){data3$pheno <- factor(data3$pheno,levels = c(3,1,2))
    }else if(n==4){data3$pheno <- factor(data3$pheno,levels = c(3,1,2,4))
    }else if(n==5){data3$pheno <- factor(data3$pheno,levels = c(3,1,2,4,5))}
  }else if(ref==4){
    if(n==2){data3$pheno <- factor(data3$pheno,levels = c(2,1))
    }else if(n==3){data3$pheno <- factor(data3$pheno,levels = c(3,1,2))
    }else if(n==4){data3$pheno <- factor(data3$pheno,levels = c(4,1,2,3))
    }else if(n==5){data3$pheno <- factor(data3$pheno,levels = c(4,1,2,3,5))}
  }else if(ref==5){
    if(n==2){data3$pheno <- factor(data3$pheno,levels = c(2,1))
    }else if(n==3){data3$pheno <- factor(data3$pheno,levels = c(3,1,2))
    }else if(n==4){data3$pheno <- factor(data3$pheno,levels = c(4,1,2,3))
    }else if(n==5){data3$pheno <- factor(data3$pheno,levels = c(5,1,2,3,4))}
  }
  class(data3$pheno);levels(data3$pheno)
  fit <- coxph(Surv(days,ACD_status) ~ .,data = data3[,-1])
  res <- cox.zph(fit)
  fit0 <- summary(fit)
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  a = n;b=n-1
  output <- data.frame(Ref=rep(ref,a),
                       exposure=rep(component,a),
                       outcome=rep("Depression",a),
                       class=levels(data3$pheno),
                       HR=c(1,coef$`exp(coef)`[1:b]),
                       LCI=c(1,confint$`lower .95`[1:b]),
                       HCI=c(1,confint$`upper .95`[1:b]),
                       P=c(1,coef$`Pr(>|z|)`[1:b]),
                       res.zph = rep(res$table[1,3],a))
  return(output)
}
output <- data.frame()
for (ref in 1:5) {
  for (i in colnames(group)[2:ncol(group)]) {
    output0 <- func_cox(i,ref)
    output <- rbind(output,output0)
  }
}
write.csv(output,"Result/Group_dementia_model2.csv",row.names = F)