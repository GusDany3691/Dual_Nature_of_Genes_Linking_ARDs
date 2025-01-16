library(dplyr)
library(igraph)
library(rje)

IntArr = c("EndPin","EndC90","EndC95","EndKgo")
IntLen = length(IntArr)

#i=4
for(i in 1:IntLen){
  print(i)
  SrcInt = IntArr[i]
  
  ### LOAD NETWROK DATA ##########################################################
  
  #NodAll = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/EndNod.rds", sep='') %>% readRDS()
  GrpAll = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/EndGrp.rds", sep='') %>% readRDS()
  #EdgAll = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/EndEdg.rds", sep='') %>% readRDS()
  
  ### GEN_GEN FRAME ##############################################################
  
  FulGenGen_DstFrm = GrpAll$FulGlbGen %>% distances() %>% as.matrix() %>% as.data.frame()
  saveRDS(FulGenGen_DstFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/FulGenGen_DstFrm.rds", sep='') )
  
  FulGenGen_DstFrm = readRDS(paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/FulGenGen_DstFrm.rds", sep='') )
  FulGenGen_PrxFrm = FulGenGen_DstFrm %>% ProximityFunction()
  saveRDS(FulGenGen_PrxFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/FulGenGen_PrxFrm.rds", sep='') )
  
  ### GEN_PHN FRAME ##############################################################
  
  #GenArr = NodAll$FulGen
  #GenLen = length(GenArr)
  
  #PhnArr = NodAll$Phn %>% unique()
  #PhnLen = length(PhnArr)
  #for(j in 1:PhnLen){
  #  paste(i,j, PhnLen) %>% print()
  #  PhnQry = PhnArr[j]
  #  PhnDel = PhnArr[PhnArr != PhnQry]
  #  QryGrp = GrpAll$FulGlbMpa %>% delete_vertices(PhnDel)
  #  sGenPhn_DstFrm = distances(QryGrp, PhnQry) %>% as.matrix() %>% data.frame()
  #  sGenPhn_DstFrm[[PhnQry]] = NULL
  #  row.names(sGenPhn_DstFrm) = PhnQry
  #  if(j == 1){
  #    GenPhn_DstFrm = sGenPhn_DstFrm
  #  } else{
  #    GenPhn_DstFrm = rbind(GenPhn_DstFrm, sGenPhn_DstFrm)
  #  }
  #}
  
  #FulGenPhn_DstFrm = GenPhn_DstFrm %>% as.data.frame() %>% t() %>% as.data.frame()
  #saveRDS(FulGenPhn_DstFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/FulGenPhn_DstFrm.rds", sep='') )
}

################################################################################
# FUNCTIONS
################################################################################

ProximityFunction = function(x){
  return(1/(1+x))
}

