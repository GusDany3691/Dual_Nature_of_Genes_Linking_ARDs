library(igraph)
library(rje)
library(dplyr)
library(tidyr)
library(tidyverse)

GenPhnComFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")
PhnComFrm = GenPhnComFrm %>% select(Phn,Com) %>% unique()

####################################################################
### COMPUTE INDIRECT INTERACTIONS ##################################
####################################################################

TypIntArr = c('EndPin','EndC90','EndC95','EndKgo')
TypIntLen = length(TypIntArr)
i=4
for(i in 1:TypIntLen){
  
  TypInt = TypIntArr[i]
  
  TxtPth = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/",sep="")
  
  NodAll = paste(TxtPth,"EndNod.rds",sep="") %>% readRDS()
  GrpAll = paste(TxtPth,"EndGrp.rds",sep="") %>% readRDS()
  EdgAll = paste(TxtPth,"EndEdg.rds",sep="") %>% readRDS()
  
  WrkGrp = GrpAll$GlbMpa %>% delete_vertices(NodAll$Phn)
  GenAdjFrm = WrkGrp %>% as_adjacency_matrix %>% as.matrix() %>% as.data.frame()
  
  GenArr = NodAll$GlbGen %>% unique()
  GenLen = length(GenArr)
  
  ### GEN-PHN ########################################################
  
  PhnArr = GenPhnComFrm$Phn %>% unique()
  PhnLen = length(PhnArr)
  GenPhn_IndFrm = matrix(0,GenLen,PhnLen) %>% as.data.frame()
  colnames(GenPhn_IndFrm) = PhnArr
  row.names(GenPhn_IndFrm) = row.names(GenAdjFrm)
  
  # Get the number of association of genes to each phenotype
  for(j in 1:PhnLen){
    paste(i,1,j) %>% print()
    PhnQry = PhnArr[j] 
    GenPhn = GenPhnComFrm %>% filter(Phn %in% PhnQry) %>% pull(Gen)
    GenPhn_IndFrm[, PhnQry] = GenAdjFrm[, GenPhn, drop = FALSE] %>% rowSums()
  }
  
  # LONG FORM OF INTERACTOR FRAME
  GenPhn_IndFrmLng <- GenPhn_IndFrm %>%
    rownames_to_column(var = "rowname") %>%   # Convert row names to a column
    pivot_longer(
      cols = -rowname,                        # Exclude the 'rowname' column from pivoting
      names_to = "variable",                  # Name for the column with original column names
      values_to = "value"                     # Name for the column with values
    ) %>% 
    filter(value>0) %>%
    rename(Gen=rowname,Phn=variable,Val=value) %>%
    merge(PhnComFrm) %>% 
    select(Gen, Phn, Com) %>%
    unique()
    
  GenPhn_IndFrmLng %>% filter(Gen %in% "PRKACA")
  
  GenPhn_IndFrmLng %>% saveRDS(paste(TxtPth,"GenPhn_IndFrmLng.rds",sep=""))
  
}

