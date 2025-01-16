library(dplyr)
library(igraph)
library(rje)

#DirGenPhnComFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")
DirGenPhnComFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")

GenAgeHum_Genes = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeHum_Genes.rds")

DirIntPhnFrm = DirGenPhnComFrm %>% 
  select(Gen,Phn) %>% 
  pull(Gen) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename(Gen='.',DirIntPhn='Freq')

DirIntComFrm = DirGenPhnComFrm %>% 
  select(Gen,Com) %>% 
  unique() %>% 
  pull(Gen) %>% 
  table() %>% 
  as.data.frame %>% 
  rename(Gen='.',DirIntCom='Freq')

TypIntArr = c('EndPin','EndC90','EndC95','EndKgo')
TypIntLen = length(TypIntArr)
i=1
for(i in 1:TypIntLen){

  print(i)
  TypInt = TypIntArr[i]
  
  TxtPth = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/",sep="")
  
  IndGenPhnComFrm = paste(TxtPth,"GenPhn_IndFrmLng.rds",sep="") %>% readRDS()
  
  IndIntPhnFrm = IndGenPhnComFrm %>% 
    select(Gen,Phn) %>% 
    pull(Gen) %>% 
    table() %>% 
    as.data.frame() %>% 
    rename(Gen='.',IndIntPhn='Freq')
  
  IndIntComFrm = IndGenPhnComFrm %>% 
    select(Gen,Com) %>% 
    unique() %>% 
    pull(Gen) %>% 
    table() %>% 
    as.data.frame() %>% 
    rename(Gen='.',IndIntCom='Freq')

  DirIndIntFrm = DirIntPhnFrm %>% 
    merge(DirIntComFrm,by='Gen',all=TRUE) %>%
    merge(IndIntPhnFrm,by='Gen',all=TRUE) %>%
    merge(IndIntComFrm,by='Gen',all=TRUE) %>% 
    mutate(GenAge = ifelse(Gen %in% GenAgeHum_Genes, "Age", "NotAge"))
  
  DirIndIntFrm[is.na(DirIndIntFrm)] = 0
  
  DirIndIntFrm %>% saveRDS(paste(TxtPth,"DirIndIntFrm.rds"))
  
}
