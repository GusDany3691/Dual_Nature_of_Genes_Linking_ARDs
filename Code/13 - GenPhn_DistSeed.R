library(dplyr)
library(igraph)
library(rje)

IntArr = c("EndPin","EndC90","EndC95","EndKgo")
IntLen = length(IntArr)

for(i in 1:IntLen){
  print(i)
  SrcInt = IntArr[i]
  
  ### LOAD NETWROK DATA ##########################################################
  
  GrpAll = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/EndGrp.rds", sep='') %>% readRDS()
  
  ### GEN_GEN FRAME ##############################################################
  
  FulGenGen_DstFrm = GrpAll$FulGlbGen %>% distances() %>% as.matrix() %>% as.data.frame()
  saveRDS(FulGenGen_DstFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/FulGenGen_DstFrm.rds", sep='') )
  
  FulGenGen_DstFrm = readRDS(paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/FulGenGen_DstFrm.rds", sep='') )
  FulGenGen_PrxFrm = FulGenGen_DstFrm %>% ProximityFunction()
  saveRDS(FulGenGen_PrxFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/FulGenGen_PrxFrm.rds", sep='') )
  
}

################################################################################
# FUNCTIONS
################################################################################

ProximityFunction = function(x){
  return(1/(1+x))
}

