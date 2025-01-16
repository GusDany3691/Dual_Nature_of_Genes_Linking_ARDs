# Frm = Frame
# Com = Community (ARC)
# Prx = Proximity
# Phn = Phenotype (ARD)
# Dst = Distance
# Cls = Closest
# Ngb = Neighbours
# Avg = Average
# Pin = Protein-protein interaction
# C90 = Coexpression 90
# C95 = Coexpression 95
# Kgo = Kegg

library(dplyr)
library(rje)

################################################################################

#GenAgeHumArr = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeHum_Genes.rds")

IntArr = c("EndPin","EndC90","EndC95","EndKgo")
IntLen = length(IntArr)

for(x in 1:IntLen){
  SrcInt = IntArr[x]

  GenGen_DstFrm = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","FulGenGen_DstFrm.rds",sep="") %>% readRDS()
  GenGen_PrxFrm = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","FulGenGen_PrxFrm.rds",sep="") %>% readRDS()
  
  GenPhnComFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")
  
  # ARDs #########################################################################
  PhnArr = GenPhnComFrm$Phn %>% unique()
  PhnLen = length(PhnArr)
  for(i in 1:PhnLen){
    paste(x,1,i) %>% print()
    PhnQry = PhnArr[i]
    PhnQryGen = GenPhnComFrm %>% filter(Phn %in% PhnQry) %>% pull(Gen)
    QryGenGen_DstFrm = GenGen_DstFrm[,PhnQryGen,drop=FALSE]
    QryGenGen_PrxFrm = GenGen_PrxFrm[,PhnQryGen,drop=FALSE]
    sGenPhnCls_DstFrm = QryGenGen_DstFrm %>% compute_column(type="Min",colname=PhnQry)
    sGenPhnCls_PrxFrm = QryGenGen_PrxFrm %>% compute_column(type="Max",colname=PhnQry)
    sGenPhnAvg_PrxFrm = QryGenGen_PrxFrm %>% compute_column(type="Men",colname=PhnQry) 
    sGenPhnNgb_PrxFrm = QryGenGen_PrxFrm %>% compute_column(type="Ngb",colname=PhnQry) 
    if(i==1){
      GenPhnCls_DstFrm = sGenPhnCls_DstFrm
      GenPhnCls_PrxFrm = sGenPhnCls_PrxFrm
      GenPhnAvg_PrxFrm = sGenPhnAvg_PrxFrm
      GenPhnNgb_PrxFrm = sGenPhnNgb_PrxFrm
    } else{
      GenPhnCls_DstFrm = cbind(GenPhnCls_DstFrm,sGenPhnCls_DstFrm)
      GenPhnCls_PrxFrm = cbind(GenPhnCls_PrxFrm,sGenPhnCls_PrxFrm)
      GenPhnAvg_PrxFrm = cbind(GenPhnAvg_PrxFrm,sGenPhnAvg_PrxFrm)
      GenPhnNgb_PrxFrm = cbind(GenPhnNgb_PrxFrm,sGenPhnNgb_PrxFrm)
    }
    
  }
  
  #GenPhnCls_DstFrm$GenAgeHum = ifelse(row.names(GenPhnCls_DstFrm) %in% GenAgeHumArr,TRUE,FALSE)
  #GenPhnCls_PrxFrm$GenAgeHum = ifelse(row.names(GenPhnCls_PrxFrm) %in% GenAgeHumArr,TRUE,FALSE)
  #GenPhnAvg_PrxFrm$GenAgeHum = ifelse(row.names(GenPhnAvg_PrxFrm) %in% GenAgeHumArr,TRUE,FALSE)
  #GenPhnNgb_PrxFrm$GenAgeHum = ifelse(row.names(GenPhnNgb_PrxFrm) %in% GenAgeHumArr,TRUE,FALSE)
  
  #write.csv(GenPhnCls_DstFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenPhnCls_DstFrm.csv", sep=""))
  #write.csv(GenPhnCls_PrxFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenPhnCls_PrxFrm.csv", sep=""))
  #write.csv(GenPhnAvg_PrxFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenPhnAvg_PrxFrm.csv", sep=""))
  #write.csv(GenPhnNgb_PrxFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenPhnNgb_PrxFrm.csv", sep=""))
  
  
  # ARCs #########################################################################
  
  ComArr = GenPhnComFrm$Com %>% unique()
  ComLen = length(ComArr)
  for(i in 1:ComLen){
    paste(x,2,i) %>% print()
    ComQry = ComArr[i]
    ComQryGen = GenPhnComFrm %>% filter(Com %in% ComQry) %>% pull(Gen) %>% unique()
    QryGenGen_DstFrm = GenGen_DstFrm[,ComQryGen,drop=FALSE]
    QryGenGen_PrxFrm = GenGen_PrxFrm[,ComQryGen,drop=FALSE]
    sGenComCls_DstFrm = QryGenGen_DstFrm %>% compute_column(type="Min",colname=ComQry)
    sGenComCls_PrxFrm = QryGenGen_PrxFrm %>% compute_column(type="Max",colname=ComQry)
    sGenComAvg_PrxFrm = QryGenGen_PrxFrm %>% compute_column(type="Men",colname=ComQry) 
    sGenComNgb_PrxFrm = QryGenGen_PrxFrm %>% compute_column(type="Ngb",colname=ComQry) 
    if(i==1){
      GenComCls_DstFrm = sGenComCls_DstFrm
      GenComCls_PrxFrm = sGenComCls_PrxFrm
      GenComAvg_PrxFrm = sGenComAvg_PrxFrm
      GenComNgb_PrxFrm = sGenComNgb_PrxFrm
    } else{
      GenComCls_DstFrm = cbind(GenComCls_DstFrm,sGenComCls_DstFrm)
      GenComCls_PrxFrm = cbind(GenComCls_PrxFrm,sGenComCls_PrxFrm)
      GenComAvg_PrxFrm = cbind(GenComAvg_PrxFrm,sGenComAvg_PrxFrm)
      GenComNgb_PrxFrm = cbind(GenComNgb_PrxFrm,sGenComNgb_PrxFrm)
    }
    
  }
  
  #GenComCls_DstFrm$GenAgeHum = ifelse(row.names(GenComCls_DstFrm) %in% GenAgeHumArr,TRUE,FALSE)
  #GenComCls_PrxFrm$GenAgeHum = ifelse(row.names(GenComCls_PrxFrm) %in% GenAgeHumArr,TRUE,FALSE)
  #GenComAvg_PrxFrm$GenAgeHum = ifelse(row.names(GenComAvg_PrxFrm) %in% GenAgeHumArr,TRUE,FALSE)
  #GenComNgb_PrxFrm$GenAgeHum = ifelse(row.names(GenComNgb_PrxFrm) %in% GenAgeHumArr,TRUE,FALSE)
  
  #write.csv(GenComCls_DstFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenComCls_DstFrm.csv", sep=""))
  #write.csv(GenComCls_PrxFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenComCls_PrxFrm.csv", sep=""))
  #write.csv(GenComAvg_PrxFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenComAvg_PrxFrm.csv", sep=""))
  #write.csv(GenComNgb_PrxFrm, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenComNgb_PrxFrm.csv", sep=""))
  
  # SAVE #######################################################################
  
  GenPhn_PrxDstLst = list(GenPhnCls_DstFrm=GenPhnCls_DstFrm, GenPhnCls_PrxFrm=GenPhnCls_PrxFrm, 
                          GenPhnAvg_PrxFrm=GenPhnAvg_PrxFrm, GenPhnNgb_PrxFrm=GenPhnNgb_PrxFrm,
                          GenComCls_DstFrm=GenComCls_DstFrm, GenComCls_PrxFrm=GenComCls_PrxFrm, 
                          GenComAvg_PrxFrm=GenComAvg_PrxFrm, GenComNgb_PrxFrm=GenComNgb_PrxFrm)
  
  saveRDS(GenPhn_PrxDstLst, paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",SrcInt,"/","GenPhn_PrxDstLst_Mng.rds", sep="") )

}



################################################################################
# FUNCTIONS
################################################################################

ProximityFunction = function(x){
  return(1/(1+x))
}

compute_column <- function(data, type="Max", colname = "Col") {
  # Compute row sums
  
  if(type=="Max"){
    row_values <- rowMaxs(data)
  }
  if(type=="Min"){
    row_values <- rowMins(data)
  }
  if(type=="Men"){
    row_values <- rowMeans(data)
  }
  if(type=="Ngb"){
    row_values <- count_row_elements_equal_to_0.5(data)
  }
  
  
  # Create a data frame with the row sums as a single row
  result <- data.frame(row_values)
  
  # Set the column names to the original row names of the input data
  row.names(result) <- rownames(data)
  
  # Set the row name for the result
  colnames(result) <- colname
  
  return(result)
}


count_row_elements_equal_to_0.5 <- function(data) {
  # Apply a function to each row that counts the number of elements equal to 0.5
  row_counts <- apply(data, 1, function(row) sum(row == 0.5))
  return(row_counts)
}
