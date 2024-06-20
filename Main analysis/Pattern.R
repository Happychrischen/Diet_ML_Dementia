##Scoring rules
#ML assisted diet score
group <- fread("Food_group.csv")
data <- select(group,eid,`Green leafy vegetables`,`Olive oil`,`Berries and citrus fruits`,`Potatoes`,`Eggs`,`Poultry`,`Sweetened beverages`)
func_ML <- function(data){
  data$`Green leafy vegetables` = ifelse(data$`Green leafy vegetables` > 0.25,1,0)
  data$`Olive oil` = ifelse(data$`Olive oil` > 0,1,0)
  data$`Berries and citrus fruits` = ifelse(data$`Berries and citrus fruits` > 0 & data$`Berries and citrus fruits` <= 2,1,0)
  data$Potatoes = ifelse(data$Potatoes > 0 & data$Potatoes <= 0.75,1,0)
  data$Eggs = ifelse(data$Eggs >0 & data$Eggs <= 1,1,0)
  data$Poultry = ifelse(data$Poultry >0 & data$Poultry <= 0.5,1,0)
  data$`Sweetened beverages` = ifelse(data$`Sweetened beverages` == 0,1,0)
  return(data)
}
score <- func_ML(data)
score$ML_score <- apply(data_score[,-1],MARGIN = 1,sum)
write.csv(score,"ML_score.csv",row.names = F)
#MIND score
data <- dplyr::select(group,eid,Berries,`Green leafy vegetables`:`Other vegetables`,Legumes,Nuts,
                           `Whole grains`,Poultry,`Fish and other seafood`,`Olive oil`,Wine,`Red meats`,`Processed meats`,
                           `Butter and margarine`,`Sweets and desserts`,`Fried and fast foods`)
meat <- rowSums(select(data_MIND,`Red meats`,`Processed meats`))
data <- select(data,-c(`Red meats`,`Processed meats`)) %>%
  mutate(`Red and processed meats` = meat,.after = Wine)
cheese <- rowSums(select(data,`102810`:`102910`))
data <- mutate(data,Cheese = cheese,.after = `Red and processed meats`)
func_ter <- function(x){
  breaks = quantile(x,probs = seq(0,1,1/3))
  if(length(unique(breaks)) == 2){a <- factor(ifelse(x>0,1,0))
  }else if(length(unique(breaks)) < 4){a <- cut(x,breaks = c(min(x)-1,unique(breaks)),include.lowest = T)
  }else{a <- cut(x,breaks = unique(breaks),include.lowest = T)}
  return(a)
}
data1 <- mutate(data,across(.cols = 2:16,func_ter))
data2 <- mutate(data1,across(.cols = 2:16,as.integer))
func_for1 <- function(x){
  a <- ifelse(x == 1,0,1)
  return(a)
}
func_for2 <- function(x){
  a <- ifelse(x == 1,0,ifelse(x == 2,0.5,1))
  return(a)
}
func_aginst1 <- function(x){
  a <- ifelse(x == 1,1,0)
  return(a)
}
func_aginst2 <- function(x){
  a <- ifelse(x == 1,1,ifelse(x == 2,0.5,0))
  return(a)
}
func_wine <- function(x){
  a <- ifelse(x == 1,1,ifelse(x == 0 | x >1,0,0.5))
  return(a)
}
data3 <- mutate(data2,across(.cols = c(2,6,10),func_for1)) %>%
  mutate(across(.cols = c(3,4,5,7,8,9,11),func_for2)) %>%
  mutate(across(.cols = c(12,13,14,15,16),func_aginst2))

data3 <- mutate(data3,Wine = func_wine(data$Wine))
data3$MIND_score <- apply(data3[,-1],MARGIN = 1,sum,na.rm = T)
write.csv(data_MIND3,"MIND_score.csv",row.names = F)