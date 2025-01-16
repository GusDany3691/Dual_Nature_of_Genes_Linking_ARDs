# Wid - Wide
# Com - Community (ARC)
# Arc - ARC
# Avl = Available 

library(dplyr)
library(tidyr)
library(rstatix)
library(reshape2)

GenPhnCom = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")

GenAgeHum = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeHum_Genes.rds")
GenAgeMod = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeMod_Genes.rds")

GenAgeFrm = data.frame(Phn="GenAge", Mng="GenAge", Com="GenAge", Gen=GenAgeHum) %>% 
  rbind( data.frame(Phn="ModAge", Mng="ModAge", Com="ModAge", Gen=GenAgeMod) )

TxtPth = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/",sep="")

# Retreieve gene specificity frame
OrgTauFrm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Specificity/Tau_gene_V8.csv")
EnsHgnFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/EnsHgnFrm.rds")
row.names(EnsHgnFrm) = EnsHgnFrm$ensembl_gene_id
OrgTauFrm$Gen = EnsHgnFrm[OrgTauFrm$gene_id,"external_gene_name"]
TauFrm = OrgTauFrm %>% select(Gen, tau) %>% filter(!is.na(tau))
rm(OrgTauFrm)

#GenPhnCom$Gen %>% intersect(NewTauFrm$Gen) %>% unique() %>% length()
#GenAgeMod %>% intersect(NewTauFrm$Gen) %>% unique() %>% length()

#colnames(TauFrm)

# LAST GENE-COMMUNITY DATA
NtwGen = TauFrm$Gen
PhnGen = GenPhnCom$Gen %>% unique()
OthGen = setdiff(NtwGen, PhnGen)
ComArr = GenPhnCom$Com %>% unique()
ComLen = length(ComArr)


################################################################################
### PLEIOTROPIC ASSOCIATIONS ###################################################
################################################################################

### GENERATE DATAFRAMES WITH SUM DATA ##########################################

# GENERATE THREE COLUMN DATAFRAME
ComGenFrm = GenPhnCom %>% select(Com, Gen) %>% unique()
ComGenFrm$Val = 1

# CONVERT THREE COLUMNS FRAME TO RECTANGULAR
WidComGenFrm <- ComGenFrm %>%
  spread(Com, Val)

# Set row names using the 'row' column and remove it from the data frame
rownames(WidComGenFrm) <- WidComGenFrm$Gen
WidComGenFrm <- WidComGenFrm %>% select(-Gen)
WidComGenFrm[is.na(WidComGenFrm)] = 0

ArcArr = WidComGenFrm %>% colnames() 


# CREATE COLUMN WITH MERGED COMMUNITIES SEPPARATED BY " - "
ComSetArr = c()
for(i in 1:nrow(WidComGenFrm)){
  print(i)
  ComSetArr[i] = colnames(WidComGenFrm[i,])[WidComGenFrm[i,] == 1] %>% paste(collapse=" - ")
}



# INCORPORATE SUM AND MERGED COMMUNITIES
WidComGenFrm$ComNum = rowSums(WidComGenFrm)
WidComGenFrm$OrgComSet = ComSetArr
WidComGenFrm$Gen = row.names(WidComGenFrm)


# MAKE NAME OF COMMUNITIES SMALLER
GenSumCom = WidComGenFrm %>% select(Gen, ComNum, OrgComSet)
OrgComSetArr = GenSumCom$OrgComSet

AbvComSetArr = full_shorten_sub_elements(GenSumCom$OrgComSet)

ComComFrm = data.frame(OrgComSet=OrgComSetArr, AbvComSet=AbvComSetArr)
ComComFrm = unique(ComComFrm)

GenSumCom$AbbComSet = AbvComSetArr

GenSumCom = GenSumCom %>% arrange(desc(ComNum))
GenSumCom = GenSumCom %>% filter(ComNum>0)

################################################################################
### COMPUTE DATA ###############################################################
################################################################################

GenPltFrm = GenSumCom %>% rename(Mng=AbbComSet, Phn=ComNum) %>% mutate(Com=paste("Plt",Phn,sep="")) %>% 
  select(Phn, Mng, Com, Gen)

NewGenPhnCom = rbind(GenPhnCom,GenPltFrm,GenAgeFrm)

NewGenPhnComTau = NewGenPhnCom %>% merge(TauFrm,by='Gen') %>% filter(!is.na(tau))

#dim(ComTauFrm )
ComTauFrm = NewGenPhnComTau %>% select(Gen, Com, tau) %>% unique() %>% rename(Scr=tau)
ComTauFrm$color_group = 'Other_ARCs'
ComTauFrm$Lbl = 'Tissue specificity in genes associated to ARCs and GenAge'
ComTauFrm$color_group[ComTauFrm$Com %in% "immunological/systemic disorders"] = "immunological/systemic disorders"
ComTauFrm$color_group[ComTauFrm$Com %in% unique(GenAgeFrm$Com)] = "immunological/systemic disorders"
ComTauFrm$color_group[ComTauFrm$Com %in% unique(GenPltFrm$Com)] = "Pleiotropies"

NewComArr = ComTauFrm$Com %>% unique()
NewComLen = length(NewComArr)

ComComFrm = data.frame(Scr  =  ComTauFrm$Scr,
                       ExtCom = ComTauFrm$Com,
                       InnCom = ComTauFrm$Com,
                       Inv = 1/ComTauFrm$Scr,
                       LogInv = log(1/ComTauFrm$Scr)
)

ComComScrFrm = ComComFrm


################################################################################
### STATISTICS #################################################################
################################################################################

ComComTstFrm = data.frame()
for(i in 1:NewComLen){
  
  ### PREPARE FRAME ############################################################

  QryExtCom = NewComArr[[i]]
  
  QryExtComTxt = QryExtCom
  
  if(QryExtCom == "immunological/systemic disorders"){
    QryExtComTxt = "immunological/systemic\n disorders"
  }
  
  ### COMPUTE COM-COM DYNAMICS #################################################
  
  sComComTstFrm = data.frame()
  for(j in 1:NewComLen){
    
    paste(i,j) %>% print()
    
    QryInnCom = NewComArr[[j]]
    
    
    # EXT.EXT
    ExtScrComCom = ComComScrFrm %>% filter(ExtCom %in% QryExtCom, InnCom %in% QryExtCom) %>% pull(Scr)
    ExtMenComCom = ExtScrComCom %>% mean()
    ExtLenComCom = ExtScrComCom %>% length()
    
    # INN.INN
    InnScrComCom = ComComScrFrm %>% filter(ExtCom %in% QryInnCom, InnCom %in% QryInnCom) %>% pull(Scr)
    InnMenComCom = InnScrComCom %>% mean()
    InnLenComCom = InnScrComCom %>% length()
    
    # MIX - EXT.INN
    MixScrComCom = ComComScrFrm %>% filter(ExtCom %in% QryExtCom, InnCom %in% QryInnCom) %>% pull(Scr)
    MixMenComCom = MixScrComCom %>% mean()
    MixLenComCom = MixScrComCom %>% length()
    
    
    # STATISTICAL TESTS
    ExtInn.Tst = tryCatch({t.test(ExtScrComCom, InnScrComCom)$p.value}, error = function(e){NA})

    # CREATE FRAME
    ssComComTstFrm =
      data.frame(ExtCom=QryExtCom,InnCom=QryInnCom,
                 ExtMen=ExtMenComCom,InnMen=InnMenComCom, MixMen=MixMenComCom,
                 ExtInn.Dif = ExtMenComCom-InnMenComCom, 
                 ExtLen = ExtLenComCom, InnLen=InnLenComCom, MixLen=MixLenComCom, 
                 ExtInn.Tst
      )
    
    sComComTstFrm = rbind(sComComTstFrm, ssComComTstFrm) 
    
  }
  
  ComComTstFrm = rbind(sComComTstFrm, ComComTstFrm)
  
}

### TEESTING ADJUSTMENTS #######################################################

# BONFERRININ CORRECTION
ComComTstFrm$ExtInn.TstBnf = p.adjust(ComComTstFrm$ExtInn.Tst, method="bonferroni")

# ORIGINAL TEST LABEL
ComComTstFrm$ExtInn.TstTxt = TstTxtFcn(ComComTstFrm$ExtInn.Tst)

# CORRECTED TEST LABEL
ComComTstFrm$ExtInn.TstBnfTxt = TstTxtFcn(ComComTstFrm$ExtInn.TstBnf)

### MEAN ARRAYS ################################################################

ComComTstFrm %>% filter(ExtCom=="ModAge",InnCom=="ModAge")

dim(ComComTstFrm)

ExtComMenFrm = ComComTstFrm %>% select(ExtCom,ExtMen) %>% unique() %>% rename(Com=ExtCom,Men=ExtMen)
row.names(ExtComMenFrm) = ExtComMenFrm$Com
ExtComMenFrm$Com=NULL
ExtComMenFrm$Men = round(ExtComMenFrm$Men,2)

InnComMenFrm = ComComTstFrm %>% select(InnCom,InnMen) %>% unique() %>% rename(Com=InnCom,Men=InnMen)
row.names(InnComMenFrm) = InnComMenFrm$Com
InnComMenFrm$Com=NULL
InnComMenFrm$Men = round(InnComMenFrm$Men,2)

################################################################################
### PLOTTING FRAMES ############################################################
################################################################################

LngComComLst = list(ExtInn=list(), ExtMix=list(), InnMix=list())
WidComComLst = list(ExtInn=list(), ExtMix=list(), InnMix=list())

### THRE COLUMNS CASE ##########################################################

# ORIGINAL TEST
LngComComLst$ExtInn$Tst = ComComTstFrm %>% select(ExtCom,InnCom,ExtInn.Tst) %>% rename(Scr=ExtInn.Tst)

# CORRECTD TEST
LngComComLst$ExtInn$TstBnf = ComComTstFrm %>% select(ExtCom,InnCom,ExtInn.TstBnf) %>% rename(Scr=ExtInn.TstBnf)

# ORIGINAL TEST LABEL
LngComComLst$ExtInn$TstTxt = ComComTstFrm %>% select(ExtCom,InnCom,ExtInn.TstTxt) %>% rename(Scr=ExtInn.TstTxt)

# CORRECTD TEST LABEL
LngComComLst$ExtInn$TstBnfTxt = ComComTstFrm %>% select(ExtCom,InnCom,ExtInn.TstBnfTxt) %>% rename(Scr=ExtInn.TstBnfTxt)

# MEAN DIFFERENCE
LngComComLst$ExtInn$Dif = ComComTstFrm %>% select(ExtCom,InnCom,ExtInn.Dif) %>% rename(Scr=ExtInn.Dif)

# MIX
LngComComLst$ExtInn$Mix = ComComTstFrm %>% select(ExtCom,InnCom,MixMen) %>% rename(Scr=MixMen)

### RECTANGULAR CASE ###########################################################

# ORIGINAL TEST 
WidComComLst$ExtInn$Tst = spread(LngComComLst$ExtInn$Tst, key=ExtCom, value=Scr)

# CORRECTED TEST 
WidComComLst$ExtInn$TstBnf = spread(LngComComLst$ExtInn$TstBnf, key=ExtCom, value=Scr)

# ORIGINAL TEST LABEL
WidComComLst$ExtInn$TstTxt = spread(LngComComLst$ExtInn$TstTxt, key=ExtCom, value=Scr)

# CORRECTED TEST LABEL
WidComComLst$ExtInn$TstBnfTxt = spread(LngComComLst$ExtInn$TstBnfTxt, key=ExtCom, value=Scr)

# MEAN DIFFERENCE
WidComComLst$ExtInn$Dif = spread(LngComComLst$ExtInn$Dif, key=ExtCom, value=Scr)

# MIX
WidComComLst$ExtInn$Mix = spread(LngComComLst$ExtInn$Mix, key=ExtCom, value=Scr)


################################################################################
### SAVING #####################################################################
################################################################################

saveRDS(ComComTstFrm, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/ComComTstFrm.rds") 
saveRDS(WidComComLst, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/WidComComLst.rds") 
saveRDS(InnComMenFrm, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/InnComMenFrm.rds") 
saveRDS(ExtComMenFrm, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/ExtComMenFrm.rds") 
saveRDS(LngComComLst, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/LngComComLst.rds") 
saveRDS(ComComScrFrm, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/ComComScrFrm.rds") 
#saveRDS(IntComComFrm, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/IntComComFrm.rds") 


################################################################################
# FUNCTIONS
################################################################################

full_shorten_sub_elements = function(input_vector){
  output_vector <- sapply(input_vector, shorten_sub_elements)
  return(as.character(output_vector))
}

# Function to shorten each sub-element
shorten_sub_elements <- function(x) {
  sub_elements <- unlist(strsplit(x, " - "))
  shortened <- substr(sub_elements, 1, 4)
  paste(shortened, collapse = "-",sep="")
}


SumComFcn = function(GenSumCom, GenMat){ 
  
  InnNtwGen = row.names(GenMat)
  
  SetArr = GenSumCom$OrgComSet %>% unique()
  SetLen = length(SetArr)
  
  ComArrFrm = data.frame()
  ComComFrm = data.frame()
  for(x in 1:SetLen){
    ExtQry = SetArr[x]
    
    ExtComSetGen = GenSumCom %>% filter(OrgComSet %in% ExtQry) %>% pull(Gen) %>% unique()
    sComComFrm = data.frame()
    
    for(y in 1:SetLen){
      
      paste(x, y) %>% print()
      InnQry = SetArr[y]
      
      InnComSetGen = GenSumCom %>% filter(OrgComSet %in% InnQry) %>% pull(Gen) %>% unique()
      
      InnComSetGenNtw = intersect(InnComSetGen, InnNtwGen)
      ExtComSetGenNtw = intersect(ExtComSetGen, InnNtwGen)
      
      UppComCom = GenMat[ExtComSetGenNtw,InnComSetGenNtw]
      LowComCom = GenMat[InnComSetGenNtw,ExtComSetGenNtw]
      
      UppComCom.Arr = UppComCom %>% unlist() %>% as.numeric()
      LowComCom.Arr = LowComCom %>% unlist() %>% as.numeric()
      
      ComCom.Arr = c(UppComCom.Arr, LowComCom.Arr)
      
      ComCom.Len = length(ComCom.Arr)
      if(ComCom.Len==0){
        ComCom.Arr = NA
      }
      
      ssComComFrm = data.frame(Scr = ComCom.Arr, ExtComSet = ExtQry, InnComSet = InnQry, Dbs=TypInt)
      
      sComComFrm = rbind(sComComFrm, ssComComFrm)
      
    }
    
    ComComFrm = rbind(ComComFrm, sComComFrm)
    
  }
  
  return(ComComFrm)
}

# Apply the function to each element of the input vector
full_count_dashes <- function(input_vector) {
  dash_counts <- sapply(input_vector, count_dashes)
  return(as.numeric(dash_counts)+1)
}

# Function to count the number of "-" symbols in a string
count_dashes <- function(x) {
  sum(strsplit(x, "")[[1]] == "-")
}

# Pvalues to astherisks
TstTxtFcn = function(Tst){
  TstTxt = rep("ns",length(Tst))
  TstTxt = ifelse(Tst  <= 5e-2, "*", TstTxt)
  TstTxt = ifelse(Tst  <= 1e-2, "**", TstTxt)
  TstTxt = ifelse(Tst  <= 1e-3, "***", TstTxt)
  TstTxt = ifelse(Tst  <= 1e-4, "****", TstTxt)
  TstTxt[is.na(TstTxt)] = "ns"
  return(TstTxt)
}
