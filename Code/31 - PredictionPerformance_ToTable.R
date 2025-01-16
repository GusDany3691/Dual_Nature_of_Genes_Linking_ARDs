library(dplyr)
library(tidyr)

AlgFrm_Avg = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/Performances/AlgFrm_Mean.csv") %>%
  filter(Dts == 'UadPinNgb') %>%
  filter(Typ == 'All') %>%
  filter(Gen == 'Age') %>%
  select(Alg, Auc,Grp) %>% 
  mutate(Typ="Avg") %>%
  mutate(Col=paste(Grp,"-",Typ,sep="")) %>%
  select(Alg,Auc,Col) %>%
  mutate(Auc=round(Auc,2))
  
AlgFrm_ML = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/Performances/AlgFrm_ML.csv") %>%
  select(Alg, Auc, Grp) %>% 
  mutate(Typ="ML") %>%
  mutate(Col=paste(Grp,"-",Typ,sep="")) %>%
  select(Alg,Auc,Col) %>%
  mutate(Auc=Auc/100)  %>%
  mutate(Auc=round(Auc,2))

ColArr_Avg = AlgFrm_Mean$Col %>% unique()
ColArr_ML = AlgFrm_ML$Col %>% unique()

AlgFrm = rbind(AlgFrm_Avg,AlgFrm_ML)

WidAlgFrm <- AlgFrm %>%
  pivot_wider(names_from = Col, values_from = Auc) %>% 
  as.data.frame()
row.names(WidAlgFrm) = WidAlgFrm$Alg
WidAlgFrm$Alg = NULL
WidAlgFrm = WidAlgFrm[c('GenPhnCls_PrxFrm','GenComCls_PrxFrm','GenPhnAvg_PrxFrm','GenComAvg_PrxFrm','GenPhnNgb_PrxFrm','GenComNgb_PrxFrm'),
                      c('EndPin-Avg','EndPin-ML','EndC90-Avg','EndC90-ML','EndC95-Avg','EndC95-ML','EndKgo-Avg','EndKgo-ML')] 
WidAlgFrm[['Mean-Avg']] = WidAlgFrm %>% select(ColArr_Avg) %>% rowMeans() %>% round(2)
WidAlgFrm[['Mean-ML']] = WidAlgFrm %>% select(ColArr_ML) %>% rowMeans() %>% round(2)
WidAlgFrm[['Mean-All']] = WidAlgFrm %>% select(ColArr_Avg,ColArr_ML) %>% rowMeans() %>% round(2)
WidAlgFrm['Mean',] = WidAlgFrm %>% colMeans() %>% round(2)

write.csv(WidAlgFrm,'C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/Performances/AlgFrm_All.csv')

