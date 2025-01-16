library(dplyr)
library(igraph)
library(grid)
library(pROC)
library(data.table)

### INITIAL PARAMETERS #########################################################################################

GenPhnComFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")

AgeGenAll = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeHum_Genes.rds")
PhnGenAll = GenPhnComFrm$Gen %>% unique()

HrcPhn = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Showcase/PhnHrc_FrmEnd.csv")

DtsTypArr = c("Mpa", "Pin", "UadPinNgb")

TypIntArr = c("EndKgo", "EndC95", "EndC90", "EndPin")

AucFrm = data.frame()
for(TypInt in TypIntArr){
  
  ##############################################################################
  
  TxtPth = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/",sep="")
  
  FrmLst = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/","GenPhn_PrxDstLst_Mng.rds", sep="") %>% readRDS()
  
  
  Nod = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/EndNod.rds", sep='') %>% readRDS()
  Grp = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/EndGrp.rds", sep='') %>% readRDS()
  
  AgeGen = AgeGenAll %>% intersect(Nod$GenPin)
  PhnGen = PhnGenAll %>% intersect(Nod$GenPin)
  
  AgeNgb = Grp$GlbGen %>% ego(order = 1, nodes = AgeGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()
  PhnNgb = Grp$GlbGen %>% ego(order = 1, nodes = PhnGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()

  aAllNgb = c(PhnNgb, AgeNgb) %>% unique()
  aExcNgb = aAllNgb %>% setdiff(AgeGen) %>% setdiff(PhnGen)
  
  AucGenLst = list("AgeGen"= AgeGen, "PhnGen"= PhnGen,
                   "AgeNgb"= AgeNgb, "PhnNgb"= PhnNgb,
                   "aAllNgb"=aAllNgb, "aExcNgb"= aExcNgb)
  
  GenMpa = Nod$GenMpa
  GenPin = Nod$GlbPin
  GenUad = c(PhnGen, GenAge) %>% unique()
  GenUadAll = c(PhnGenAll, AgeGenAll) %>% unique()
  GenUadPin = intersect(GenUad, GenPin)
  GenUadPinNgb = ego(Grp$GlbGen, order = 1, nodes = GenUadPin, mode = "all", mindist = 0) %>% unlist() %>% names() %>% unique()
  
  DtsFrmLst = list()
  
  ### GenPhnAvg_PrxFrm #################################################
  
  DtsFrmLst$GenPhnAvg_PrxFrm = list()
  
  GenPhnAvg_PrxFrm = FrmLst$GenPhnAvg_PrxFrm[GenMpa,]
  RowGenArr = row.names(GenPhnAvg_PrxFrm)
  GenPhnAvg_PrxFrm = GenPhnAvg_PrxFrm %>% mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "CR", "NotCR"))
  DtsFrmLst$GenPhnAvg_PrxFrm$Mpa = GenPhnAvg_PrxFrm
  DtsFrmLst$GenPhnAvg_PrxFrm$Pin = GenPhnAvg_PrxFrm[GenPin,]
  DtsFrmLst$GenPhnAvg_PrxFrm$UadPinNgb = GenPhnAvg_PrxFrm[GenUadPinNgb,]
  
  for(DtsTyp in DtsTypArr){
    sAucFrm = AucFcn(DtsFrmLst$GenPhnAvg_PrxFrm[[DtsTyp]], AucGenLst, TypInt, Alg="GenPhnAvg_PrxFrm", Dts=DtsTyp)
    AucFrm = rbind(AucFrm, sAucFrm)
  }
  
  ### GenPhnCls_PrxFrm #################################################
  
  DtsFrmLst$GenPhnCls_PrxFrm = list()
  
  GenPhnCls_PrxFrm = FrmLst$GenPhnCls_PrxFrm[GenMpa,]
  RowGenArr = row.names(GenPhnCls_PrxFrm)
  GenPhnCls_PrxFrm = GenPhnCls_PrxFrm %>% mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "CR", "NotCR"))
  DtsFrmLst$GenPhnCls_PrxFrm$Mpa = GenPhnCls_PrxFrm
  DtsFrmLst$GenPhnCls_PrxFrm$Pin = GenPhnCls_PrxFrm[GenPin,]
  DtsFrmLst$GenPhnCls_PrxFrm$UadPinNgb = GenPhnCls_PrxFrm[GenUadPinNgb,]
  
  for(DtsTyp in DtsTypArr){
    sAucFrm = AucFcn(DtsFrmLst$GenPhnCls_PrxFrm[[DtsTyp]], AucGenLst, TypInt, Alg="GenPhnCls_PrxFrm", Dts=DtsTyp)
    AucFrm = rbind(AucFrm, sAucFrm)
  }
  
  ### GenComCls_PrxFrm #############################################
  
  DtsFrmLst$GenComCls_PrxFrm = list()
  
  GenComCls_PrxFrm = FrmLst$GenComCls_PrxFrm[GenMpa,]
  RowGenArr = row.names(GenComCls_PrxFrm)
  GenComCls_PrxFrm = GenComCls_PrxFrm %>% mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "CR", "NotCR"))
  DtsFrmLst$GenComCls_PrxFrm$Mpa = GenComCls_PrxFrm
  DtsFrmLst$GenComCls_PrxFrm$Pin = GenComCls_PrxFrm[GenPin,]
  DtsFrmLst$GenComCls_PrxFrm$UadPinNgb = GenComCls_PrxFrm[GenUadPinNgb,]
  
  for(DtsTyp in DtsTypArr){
    sAucFrm = AucFcn(DtsFrmLst$GenComCls_PrxFrm[[DtsTyp]], AucGenLst, TypInt, Alg="GenComCls_PrxFrm", Dts=DtsTyp)
    AucFrm = rbind(AucFrm, sAucFrm)
  }
  
  
  ### GenComAvg_PrxFrm #################################################
  
  DtsFrmLst$GenComAvg_PrxFrm = list()
  
  GenComAvg_PrxFrm = FrmLst$GenComAvg_PrxFrm[GenMpa,]
  RowGenArr = row.names(GenComAvg_PrxFrm)
  GenComAvg_PrxFrm = GenComAvg_PrxFrm %>% mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "CR", "NotCR"))
  DtsFrmLst$GenComAvg_PrxFrm$Mpa = GenComAvg_PrxFrm
  DtsFrmLst$GenComAvg_PrxFrm$Pin = GenComAvg_PrxFrm[GenPin,]
  DtsFrmLst$GenComAvg_PrxFrm$UadPinNgb = GenComAvg_PrxFrm[GenUadPinNgb,]

  for(DtsTyp in DtsTypArr){
    sAucFrm = AucFcn(DtsFrmLst$GenComAvg_PrxFrm[[DtsTyp]], AucGenLst, TypInt, Alg="GenComAvg_PrxFrm", Dts=DtsTyp)
    AucFrm = rbind(AucFrm, sAucFrm)
  }
  
    ### Sec.Ord.Path.Phn ###############################################
  
  DtsFrmLst$GenPhnNgb_PrxFrm = list()
  
  GenPhnNgb_PrxFrm = FrmLst$GenPhnNgb_PrxFrm[GenMpa,]
  RowGenArr = row.names(GenPhnNgb_PrxFrm)
  GenPhnNgb_PrxFrm = GenPhnNgb_PrxFrm %>% mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "CR", "NotCR"))
  DtsFrmLst$GenPhnNgb_PrxFrm$Mpa = GenPhnNgb_PrxFrm
  DtsFrmLst$GenPhnNgb_PrxFrm$Pin = GenPhnNgb_PrxFrm[GenPin,]
  DtsFrmLst$GenPhnNgb_PrxFrm$UadPinNgb = GenPhnNgb_PrxFrm[GenUadPinNgb,]
  
  for(DtsTyp in DtsTypArr){
    sAucFrm = AucFcn(DtsFrmLst$GenPhnNgb_PrxFrm[[DtsTyp]], AucGenLst, TypInt, Alg="GenPhnNgb_PrxFrm", Dts=DtsTyp)
    AucFrm = rbind(AucFrm, sAucFrm)
  }
  
  ### Sec.Ord.Path.Com ###############################################
  
  DtsFrmLst$GenComNgb_PrxFrm = list()
  
  GenComNgb_PrxFrm = FrmLst$GenComNgb_PrxFrm[GenMpa,]
  RowGenArr = row.names(GenComNgb_PrxFrm)
  GenComNgb_PrxFrm = GenComNgb_PrxFrm %>% mutate(Class = ifelse(RowGenArr %in% AgeGenAll, "CR", "NotCR"))
  DtsFrmLst$GenComNgb_PrxFrm$Mpa = GenComNgb_PrxFrm
  DtsFrmLst$GenComNgb_PrxFrm$Pin = GenComNgb_PrxFrm[GenPin,]
  DtsFrmLst$GenComNgb_PrxFrm$UadPinNgb = GenComNgb_PrxFrm[GenUadPinNgb,]
  
  for(DtsTyp in DtsTypArr){
    sAucFrm = AucFcn(DtsFrmLst$GenComNgb_PrxFrm[[DtsTyp]], AucGenLst, TypInt, Alg="GenComNgb_PrxFrm", Dts=DtsTyp)
    AucFrm = rbind(AucFrm, sAucFrm)
  }
  
  ### SAVE #####################################################################
  
  DtsNamArr = DtsFrmLst %>% names()
  
  GenSetArr = c("Mpa", "Pin", "UadPinNgb")
  DtsNamLen = length(DtsNamArr)
  GenSetLen = length(GenSetArr)
  
  for(i in 1:DtsNamLen){
    DtsNamQry = DtsNamArr[i]
    for(j in 1:GenSetLen){
      paste(TypInt,i,j) %>% print()
      GenSetQry = GenSetArr[j]
      
      QryFrm = DtsFrmLst[[DtsNamQry]][[GenSetQry]]
      RowNam = row.names(QryFrm)
      
      ClsArr = QryFrm$Class
      QryFrm$Class = NULL
      QryFrm$AgeGen = NULL
      QryFrm$PhnGen = NULL
      MenFrm = rowMeans(QryFrm) %>% as.data.frame()
      MenFrm$Gen = row.names(MenFrm)
      MenFrm$Class = ClsArr
      colnames(MenFrm) = c("Scr","Gen","Class")
      MenFrm = MenFrm %>% select(Gen,Scr,Class)
      
      QryFrm = QryFrm %>% mutate(AgeGen = ifelse(RowNam %in% AgeGen, "Age", "NotAge"))
      

      AucScr1 = AucFrm %>% filter(Typ=="All", Dts==GenSetQry,  Grp==TypInt, Gen=="Age", Alg==DtsNamQry) %>% unique() %>% pull(Auc)
      AucScr = (AucScr1*1000) %>% round(0)
      
      QryFrm %>% write.csv( paste(TxtPth,"Dat/",GenSetQry,"/",DtsNamQry,".csv",sep="") )
      MenFrm %>% write.csv( paste(TxtPth,"Dat/",GenSetQry,"/AUC",AucScr,DtsNamQry,".csv",sep="") )
      
    }
  }
} 


AucFrm %>% filter(Gen == "Age", Typ=="All", Dts=="UadPinNgb") %>% arrange() %>% unique() 

write.csv(AucFrm, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Performances/Mean_AlgFrm.csv" )

##############################################################################
# FUNCTIONS
##############################################################################


AucFcn = function(AllFrm, AucGenLst, TypInt, Alg, Dts){
  
  AgeGen = AucGenLst$AgeGen
  PhnGen = AucGenLst$PhnGen
  aExcNgb = AucGenLst$aExcNgb

  AllFrm$Class = NULL
  
  AvgFrm = AllFrm %>% rowMeans() %>% as.data.frame()
  
  colnames(AvgFrm) = "Scr"
  
  AvgFrm$AgeGen = 0
  AvgFrm$PhnGen = 0
  
  AvgFrm[AgeGen,"AgeGen"] = 1
  AvgFrm[PhnGen,"PhnGen"] = 1
  AvgNamQry = "Scr"
  
  # AGE
  AucArr.AgeAll = auc(AvgFrm[,'AgeGen'], AvgFrm[,AvgNamQry]) %>% as.numeric()
  AucArr.AgePhn = auc(AvgFrm[c(AgeGen, PhnGen),'AgeGen'], AvgFrm[c(AgeGen, PhnGen),AvgNamQry]) %>% as.numeric()
  AucArr.AgeOth = auc(AvgFrm[c(AgeGen, aExcNgb),'AgeGen'], AvgFrm[c(AgeGen, aExcNgb),AvgNamQry]) %>% as.numeric()
  
  # PHN
  AucArr.PhnAll = auc(AvgFrm[,'PhnGen'], AvgFrm[,AvgNamQry]) %>% as.numeric()
  AucArr.aPhnOth = auc(AvgFrm[c(PhnGen, aExcNgb),'PhnGen'], AvgFrm[c(PhnGen, aExcNgb),AvgNamQry]) %>% as.numeric()
  
  AucFrm.AgeAll = AucArr.AgeAll %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate(Gen = "Age", Typ = "All", Grp = TypInt)
  AucFrm.AgePhn = AucArr.AgePhn %>% as.data.frame() %>% t() %>% as.data.frame() %>%  mutate(Gen = "Age", Typ = "Phn", Grp = TypInt)
  AucFrm.AgeOth = AucArr.AgeOth %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate(Gen = "Age", Typ = "Oth", Grp = TypInt)
  
  AucFrm.PhnAll = AucArr.PhnAll %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate(Gen = "Phn", Typ = "All", Grp = TypInt)
  AucFrm.aPhnOth = AucArr.aPhnOth %>% as.data.frame() %>% t() %>% as.data.frame() %>% mutate(Gen = "Phn", Typ = "aOth", Grp = TypInt)
  
  sAucFrm = AucFrm.AgeAll %>% rbind(AucFrm.AgePhn) %>% rbind(AucFrm.AgeOth) %>%
    rbind(AucFrm.PhnAll) %>% rbind(AucFrm.aPhnOth) 
  
  colnames(sAucFrm) = c("Auc", "Gen", "Typ", "Grp")        
  
  sAucFrm$Dts = Dts
  sAucFrm$Alg = Alg

  sAucFrm = sAucFrm %>% select(Alg, Auc, Gen, Grp, Typ, Dts)
  
  return(sAucFrm)
  
}
