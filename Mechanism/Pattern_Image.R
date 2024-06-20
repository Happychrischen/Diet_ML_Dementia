##Linear associations between ML_score and brain structures
image1 <- fread("UKB_2_0_aparc_volume+area+thick/UKB_2_0_aparc_ThickAvg.csv")
# image2 <- fread("UKB_2_0_aparc_volume+area+thick/UKB_2_0_aparc_SurfArea.csv")
# image3 <- fread("UKB_2_0_aparc_volume+area+thick/UKB_2_0_aparc_GrayVol.csv")
# image4 <- fread("UKB_Subcortical.csv")
# image5 <- fread("UKB_DTI_FA.csv")
ML_score <- fread("New_score2.csv")
TIC <- fread("UKB_Brain_TIC.csv")
covariate <- fread("covariate_mice.csv")
outcome <- fread("Dementia.csv")
data <- merge(ML_score,image1,by="eid")
# data <- merge(ML_score,image2,by="eid")
# data <- merge(ML_score,image3,by="eid")
# data <- merge(ML_score,image4,by="eid")
# data <- merge(ML_score,image5,by="eid")
data <- merge(data,TIC,by = "eid")
data <- merge(data,covariate,by = "eid")
data <- merge(outcome,data,by = "eid")
#Linear regression
library(purrr)
reg_result <- data.frame(exposure=character(),outcome=character(),estimate = numeric(),p=numeric())
for (i in c(1:27)){
  exposure <- as.data.frame(scale(data[,c(13)]))
  outcome_name <- colnames(data)[i+13]
  outcome <- scale(data[[i+13]])
  reg <- map(exposure,~ lm(outcome~.x + TIC + energy + age + sex + ethnic + townsend + 
                             education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,data))
  # reg <- map(exposure,~ lm(outcome~.x + energy + age + sex + ethnic + townsend + 
  #                            education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,data))
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
reg_result <- mutate(reg_result,IDP = "Thickness")
# reg_result <- mutate(reg_result,IDP = "Area")
# reg_result <- mutate(reg_result,IDP = "Volume")
# reg_result <- mutate(reg_result,IDP = "Subcortex")
# reg_result <- mutate(reg_result,IDP = "DTI-FA")
write.csv(reg_result,"Pattern_image_thickness.csv",row.names = F)
# write.csv(reg_result,"Pattern_image_area.csv",row.names = F)
# write.csv(reg_result,"Pattern_image_volume.csv",row.names = F)
# write.csv(reg_result,"Pattern_image_subcortex.csv",row.names = F)
# write.csv(reg_result,"Pattern_image_FA.csv",row.names = F)

#Associations between ML components and brain structure
reg_result <- data.frame(exposure=character(),outcome=character(),estimate = numeric(),p=numeric())
for (i in c(1:27)){
  exposure <- as.data.frame(scale(data[,6:12]))
  outcome_name <- colnames(data)[i+13]
  outcome <- scale(data[[i+13]])
  reg <- map(exposure,~ lm(outcome~.x + TIC + energy + age + sex + ethnic + townsend + 
                             education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,data))
  # reg <- map(exposure,~ lm(outcome~.x + energy + age + sex + ethnic + townsend + 
  #                            education + smoking + activity + BMI + APOE4 + diabetes + hypertension + stroke + CVD,data))
  pvalue <- map_df(reg,~summary(.x)[["coefficients"]][2,4])
  estimate <- map_df(reg,~summary(.x)[["coefficients"]][2,1])
  reg_result0 <- rbind(estimate,pvalue)
  reg_result0 <- as.data.frame(t(reg_result0))
  colnames(reg_result0) <- c("estimate","p")
  reg_result0 <- mutate(reg_result0,exposure = rownames(reg_result0),outcome = outcome_name,
                        .before = estimate)
  reg_result <- rbind(reg_result,reg_result0)
  
}
FDR <- function(p,n){
  row <- length(p)
  for (i in 1:n) {
    p[seq(i,row,n)] <- p.adjust(p[seq(i,row,n)],method = "fdr")
  }
  return(p)
}
reg_result$fdr <- FDR(reg_result$p,7)
write.csv(reg_result,"Component_image_thickness.csv",row.names = F)
# write.csv(reg_result,"Component_image_area.csv",row.names = F)
# write.csv(reg_result,"Component_image_volume.csv",row.names = F)
# write.csv(reg_result,"Component_image_subcortex.csv",row.names = F)
# write.csv(reg_result,"Component_image_FA.csv",row.names = F)