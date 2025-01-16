# Fct - Factor
# Phn - Phenotype (ARC)
# Com - Community (ARC)
# Grp - Graphs
# Nod - Nodes 
# Edg - Edges
# Pin = Protein interactions
# Pin - PPI
# C90 - Cox90
# C95 - Cox95
# Kgo - KEGG
# Tmp - Temporal
# Qry - Query
# Frm - Frame
# Age - GenAgeHum
# Mod - GenAgeMod
# a (pre-index) - GenAgeHum-related
# m (pre-index) - GenAgeMod-related
# Htm - Heatmap
# Ngb - Neighbour
# Oth - Othond (Others)
# Dis - Disease
# Cat - Category
# Srt - Sort
# Fcn - Function
# Leg - Legend
# Flt = Flattered (Matrix to vector)

library(rje)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(igraph)

GenPhnCom = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")
AgnPhnHrc_Frm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Showcase/PhnHrc_FrmEnd.csv")
GenTyp="Pin"

TypIntArr = c("EndPin", "EndC95", "EndC90", "EndKgo")
TypIntLen = length(TypIntArr)
SepFct = c("Age", "Mod")
HtmLst =list()
x=2
for(x in 1:TypIntLen){
  print(x)
  TypInt = TypIntArr[x]
  TxtPth = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/",sep="")
  
  Nod = readRDS(paste(TxtPth,"EndNod",".rds",sep=""))
  Edg = readRDS(paste(TxtPth,"EndEdg",".rds",sep=""))
  Grp = readRDS(paste(TxtPth,"EndGrp",".rds",sep=""))
  QryGrp = Grp$GlbGen
  
  PrxQryFrm = read.csv(paste(TxtPth, "Dat/",GenTyp,"/", "GenComCls_PrxFrm.csv",sep="")) %>% as.data.frame()
  PrxQryFrm = PrxQryFrm %>% filter(!(is.na(X)))
  row.names(PrxQryFrm) = PrxQryFrm$X
  PrxQryFrm$X = NULL 
  PrxQryFrm$AgeGen = NULL 
  
  # GET DISTANCE PLOTTING FRAME
  DstQryFrm = (1/PrxQryFrm)-1 # DstFrm
  QryFrm = DstQryFrm 
  MaxVal = QryFrm[QryFrm!=Inf] %>% max()
  QryFrm = abs(QryFrm-MaxVal)
  QryFrm[QryFrm==Inf] = -1
  
  AllFrm = QryFrm %>% as.matrix() %>% round() %>% t()   # AllFrm
  AllGen = colnames(AllFrm)
  
  # Genes of Diseases and GenAge
  DisGen = Edg$GenPhn %>% filter(!(Phn %in% c("HumAge","ModAge"))) %>% pull(Gen) %>% unique() %>% intersect(AllGen) 
  AgeGen = Edg$GenPhn %>% filter(Phn == "HumAge") %>% pull(Gen) %>% unique() %>% intersect(AllGen) 
  ModGen = Edg$GenPhn %>% filter(Phn == "ModAge") %>% pull(Gen) %>% unique() %>% intersect(AllGen)
  
  # GENE-GENE DISTANCES
  DstFrm = QryGrp %>% distances() %>% as.matrix() %>% as.data.frame()
  DstFrm = DstFrm[AllGen, AllGen]
  
  # GENE NEIGHBOURS
  DisNgb = QryGrp %>% ego(order = 1, nodes = DisGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()
  AgeNgb = QryGrp %>% ego(order = 1, nodes = AgeGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()
  ModNgb = QryGrp %>% ego(order = 1, nodes = ModGen, mode = c("all"), mindist = 0) %>% unlist() %>% names() %>% unique()
  
  # GENE DISTANCES TO AGEING AND DISEASES
  AgeDst = DstFrm[AllGen,AgeGen] %>% rowMins()
  ModDst = DstFrm[AllGen,ModGen] %>% rowMins()
  DisDst = DstFrm[AllGen,DisGen] %>% rowMins()
  
  # ------------------------------------------------------------------------------
  
  # GENE MEAN DISTANCE
  TmpGenMen = rowMeans(PrxQryFrm, na.rm=TRUE)
  GenMen = abs(TmpGenMen-max(TmpGenMen))
  DstGenMen = rowMeans(DstQryFrm, na.rm=TRUE)
  
  # GENE ARC-INTERACTORS COUNT
  GenNg1 = DstQryFrm %>% apply(1, function(x){sum(x==1)})
  
  ### ROW SPLITTING ################################################################
  
  # ALL AND EXLCUISE NEIGHBOURS
  aAllNgb = c(DisNgb, AgeNgb) %>% unique()
  #aExcNgb = aAllNgb %>% setdiff(AgeGen) %>% setdiff(DisGen)
  #mAllNgb = c(DisNgb, ModNgb) %>% unique()
  #mExcNgb = mAllNgb %>% setdiff(ModGen) %>% setdiff(DisGen)
  
  #aAllGenAlt = DisNgb %>% c(AgeGen) %>% c(DisGen) %>% unique()
  #mAllGenAlt = DisNgb %>% c(ModGen) %>% c(DisGen) %>% unique()
  #aOthGenAlt = setdiff(AllGen, aAllGenAlt)
  #mOthGenAlt = setdiff(AllGen, mAllGenAlt)
  aDisAgeNgbGen = DisNgb %>% c(AgeGen) %>% c(DisGen) %>% unique()
  mDisAgeNgbGen = DisNgb %>% c(ModGen) %>% c(DisGen) %>% unique()
  aOthGen = setdiff(AllGen, aDisAgeNgbGen)
  mOthGen = setdiff(AllGen, mDisAgeNgbGen)
  
  # DIFFERENTIAL DATA
  AllColNam = AllFrm %>% colnames()
  DisGen_NonAge = setdiff(DisGen, AgeGen)
  DisGen_NonMod = setdiff(DisGen, ModGen)
  AgeGen_NonDis = setdiff(AgeGen, DisGen)
  ModGen_NonDis = setdiff(ModGen, DisGen)
  
  AllColFct = AllColNam
  
  
  if("Age" %in% SepFct){
    AllColFct = ifelse(AllColFct %in% AgeGen, "Age", AllColFct)
    AllColFct = ifelse(AllColFct %in% DisGen_NonAge, "Dis", AllColFct)
    
    AllColFct = ifelse((AllColFct %in% aAllNgb), "Ngb", AllColFct)
    AllColFct = ifelse(!(AllColFct %in% c("Dis", "Age", "Ngb")), "Oth", AllColFct)
    
    AllColFct = as.factor(AllColFct)
  }
  if("Mod" %in% SepFct){
    AllColFct = ifelse(AllColFct %in% ModGen, "Mod", AllColFct)
    AllColFct = ifelse(AllColFct %in% DisGen_NonMod, "Dis", AllColFct)
    
    AllColFct = ifelse((AllColFct %in% aAllNgb), "Ngb", AllColFct)
    AllColFct = ifelse(!(AllColFct %in% c("Dis", "Mod", "Ngb")), "Oth", AllColFct)
    
    AllColFct = as.factor(AllColFct)
  }
  
  #################################
  
  GenTypFrm = data.frame(Gen=AgeGen, GenTyp="AgeGen") %>% 
    rbind( data.frame(Gen=ModGen, GenTyp="ModGen") ) %>%
    rbind( data.frame(Gen=DisGen, GenTyp="DisGen") ) %>%
    rbind( data.frame(Gen=DisGen, GenTyp="PhnGen") ) %>%
    rbind( data.frame(Gen=AgeNgb, GenTyp="AgeNgb") ) %>%
    rbind( data.frame(Gen=ModNgb, GenTyp="ModNgb") ) %>%
    rbind( data.frame(Gen=DisNgb, GenTyp="DisNgb") ) %>%
    rbind( data.frame(Gen=DisNgb, GenTyp="PhnNgb") ) %>%
    rbind( data.frame(Gen=AllGen, GenTyp="AllGen") ) %>% 
    #rbind( data.frame(Gen=aOthGenAlt, GenTyp="aOthGenAlt") ) %>%
    #rbind( data.frame(Gen=mOthGenAlt, GenTyp="mOthGenAlt") ) 
    rbind( data.frame(Gen=aOthGen, GenTyp="aOthGen") ) %>%
    rbind( data.frame(Gen=mOthGen, GenTyp="mOthGen") ) 
    
    
  # CATEGORICAL FRAME
  GenCatFrm = data.frame("Gen"=AllColNam, "Cat"=AllColFct)
  
  ###  COLUMNS SPLITTING ###########################################################
  
  PhnNam = AllFrm %>% row.names()
  
  AllRowFct = PhnNam %>% as.factor
  AllPhnFct = AgnPhnHrc_Frm %>% filter(RutMng %in% PhnNam) %>% pull(Mng) %>% as.factor() 
  
  
  ## ANNOTATION SETS ##########################################################
  
  AgeFrm = AllFrm[,AgeGen] %>% as.matrix() %>% round() %>% t() 
  ModFrm = AllFrm[,ModGen] %>% as.matrix() %>% round() %>% t()  
  DisFrm = AllFrm[,DisGen] %>% as.matrix() %>% round() %>% t()
  NgbFrm = AllFrm[,DisNgb] %>% as.matrix() %>% round() %>% t()
  
  setdiff(DisNgb,colnames(AllFrm))
  
  ### GENE-GENE_LABEL FRAMES FOR PLOTTING #####################################
  
  PltAgeFrm = AgeFrm
  PltModFrm = ModFrm
  PltDisFrm = DisFrm
  PltNgbFrm = NgbFrm
  
  PltAgeGenOrg = row.names(PltAgeFrm)
  PltModGenOrg = row.names(PltModFrm)
  PltDisGenOrg = row.names(PltDisFrm)
  PltNgbGenOrg = row.names(PltNgbFrm)
  
  row.names(PltAgeFrm) = paste("Age",1:nrow(PltAgeFrm),sep="")
  row.names(PltModFrm) = paste("Mod",1:nrow(PltModFrm),sep="")
  row.names(PltDisFrm) = paste("Dis",1:nrow(PltDisFrm),sep="")
  row.names(PltNgbFrm) = paste("Ngb",1:nrow(PltNgbFrm),sep="")
  
  PltAllFrm = PltAgeFrm %>% rbind(PltModFrm) %>% rbind(PltDisFrm) %>% rbind(PltNgbFrm) %>% t()
  colnames(PltAllFrm)[1:(nrow(PltAgeFrm)+nrow(PltModFrm)+1)]
  
  PltAllColFct = c(rep("Age",nrow(PltAgeFrm)), rep("Mod",nrow(PltModFrm)), rep("Dis",nrow(PltDisFrm)), rep("Ngb",nrow(PltNgbFrm))) %>% as.factor()
  
  PltAgeGenNew = row.names(PltAgeFrm)
  PltModGenNew = row.names(PltModFrm)
  PltDisGenNew = row.names(PltDisFrm)
  PltNgbGenNew = row.names(PltNgbFrm)
  
  # GENE AND LABEL RELATIONSHIP FRAMES
  PltAgeGenFrm = data.frame(Gen = PltAgeGenOrg, New = PltAgeGenNew)
  PltModGenFrm = data.frame(Gen = PltModGenOrg, New = PltModGenNew)
  PltDisGenFrm = data.frame(Gen = PltDisGenOrg, New = PltDisGenNew)
  PltNgbGenFrm = data.frame(Gen = PltNgbGenOrg, New = PltNgbGenNew)
  row.names(PltAgeGenFrm) = PltAgeGenOrg
  row.names(PltModGenFrm) = PltModGenOrg
  row.names(PltDisGenFrm) = PltDisGenOrg
  row.names(PltNgbGenFrm) = PltNgbGenOrg
  
  # BIND FRAMES
  PltAllGenFrm = PltAgeGenFrm %>% rbind(PltModGenFrm) %>% rbind(PltDisGenFrm) %>% rbind(PltNgbGenFrm)
  
  ### GENE_LABEL FRAME WITI DISTANCES ###########################################
  
  PltAgeAllFrm = PltAllFrm[,PltAllColFct%in%"Age"]
  PltModAllFrm = PltAllFrm[,PltAllColFct%in%"Mod"]
  PltDisAllFrm = PltAllFrm[,PltAllColFct%in%"Dis"]
  PltNgbAllFrm = PltAllFrm[,PltAllColFct%in%"Ngb"]
  
  NumPltGen = ncol(PltAgeAllFrm)
  
  PltNewDisAllFrm = PltDisAllFrm[, sample(1:ncol(PltDisAllFrm), NumPltGen, replace = TRUE) ]
  PltNewAgeAllFrm = PltAgeAllFrm
  PltNewModAllFrm = PltModAllFrm[, sample(1:ncol(PltModAllFrm), NumPltGen, replace = TRUE) ]
  PltNewNgbAllFrm = PltNgbAllFrm[, sample(1:ncol(PltNgbAllFrm), NumPltGen, replace = TRUE) ]
  
  # ------------------------------------------------------------------------------
  
  GenTypCol = CatFrmColFcn(GenCatFrm, Pal="Pastel 1") 
  TmpGenTypCol = GenTypCol %>% unique()
  
  PltAllFrm = PltNewDisAllFrm %>% cbind(PltNewAgeAllFrm) %>% cbind(PltNewModAllFrm) %>% cbind(PltNewNgbAllFrm)
  PltAllColFct = c(rep("Dis",NumPltGen), rep("Age",NumPltGen), rep("Mod",NumPltGen), rep("Ngb",NumPltGen) )
  PltGenTypCol = c( rep(TmpGenTypCol[1],NumPltGen) , rep(TmpGenTypCol[2],NumPltGen) , rep(TmpGenTypCol[3],NumPltGen) )
  PltGen = colnames(PltAllFrm)
  
  ### LABELLING GENES BY GENE TYPE ###############################################
  
  PltAllGenFrm = PltAgeGenFrm %>% rbind(PltModGenFrm) %>% rbind(PltDisGenFrm) %>% rbind(PltNgbGenFrm)
  
  PltDisMen = GenMen[PltDisGenFrm$Gen]
  PltAgeMen = GenMen[PltAgeGenFrm$Gen]
  PltModMen = GenMen[PltModGenFrm$Gen]
  PltNgbMen = GenMen[PltNgbGenFrm$Gen]
  
  PltDstDisMen = DstGenMen[PltDisGenFrm$Gen]
  PltDstAgeMen = DstGenMen[PltAgeGenFrm$Gen]
  PltDstModMen = DstGenMen[PltModGenFrm$Gen]
  PltDstNgbMen = DstGenMen[PltNgbGenFrm$Gen]
  
  PltDisNg1 = GenNg1[PltDisGenFrm$Gen]
  PltAgeNg1 = GenNg1[PltAgeGenFrm$Gen]
  PltModNg1 = GenNg1[PltModGenFrm$Gen]
  PltNgbNg1 = GenNg1[PltNgbGenFrm$Gen]
  
  names(PltDisMen) = PltDisGenFrm$New
  names(PltAgeMen) = PltAgeGenFrm$New
  names(PltModMen) = PltModGenFrm$New
  names(PltNgbMen) = PltNgbGenFrm$New
  
  names(PltDstDisMen) = PltDisGenFrm$New
  names(PltDstAgeMen) = PltAgeGenFrm$New
  names(PltDstModMen) = PltModGenFrm$New
  names(PltDstNgbMen) = PltNgbGenFrm$New
  
  names(PltDisNg1) = PltDisGenFrm$New
  names(PltAgeNg1) = PltAgeGenFrm$New
  names(PltModNg1) = PltModGenFrm$New
  names(PltNgbNg1) = PltNgbGenFrm$New
  
  PltPrxGenMen = c(PltDisMen, PltAgeMen, PltModMen, PltNgbMen) #PltPrxGenMen
  PltDstGenMen = c(PltDstDisMen, PltDstAgeMen, PltDstModMen, PltDstNgbMen) #PltDstGenMen
  PltGenNg1 = c(PltDisNg1, PltAgeNg1, PltModNg1, PltNgbNg1)
  
  ### Fin cabos perdidos ############################################################## 
  
  PltGenNg1 = PltGenNg1[PltGen]
  PltGenMen = PltPrxGenMen[PltGen]
  PltDstGenMen = PltDstGenMen[PltGen]
  UseAllFrm = PltAllFrm
  
  ### COLOURS ####################################################################
  
  NseFrm = (PltAllFrm + runif(length(PltAllFrm), min = 0, max = 1e-10)) %>% t()
  AllDstFltCol = CatHtmColFcn(UseAllFrm, Pal1 = c("blue", "white"), Pal2 = c("white", "red"), MidNum = 2)
  AllDstFltCol[AllDstFltCol=="#000000"] = "#464646"
  
  # CATEGORICAL
  GenTypCol = CatFrmColFcn(GenCatFrm, Pal="Pastel 1")
  ComTypCol = CatArrColFcn(AllRowFct, Pal="Dynamic")
  
  ### CONTINUOUS COLORS ##########################################################
  
  Col_Ng1 = colorRamp2(c(min(PltGenNg1, na.rm=TRUE), max(PltGenNg1, na.rm=TRUE)/2, max(PltGenNg1, na.rm=TRUE)), c("blue", "white", "red"))
  Col_Men = colorRamp2(c(min(PltPrxGenMen, na.rm=TRUE), max(PltPrxGenMen, na.rm=TRUE)/2, max(PltPrxGenMen, na.rm=TRUE)), c("blue", "white", "red"))
  Col_Dst = colorRamp2(c(min(PltPrxGenMen, na.rm=TRUE), max(PltPrxGenMen, na.rm=TRUE)/2, max(PltPrxGenMen, na.rm=TRUE)), c("blue", "white", "red"))
  
  
  ### DISTANCES ##################################################################
  
  PltGenMenDst = PltDstGenMen
  NanPltGenMenDst = PltGenMenDst
  NanPltGenMenDst[NanPltGenMenDst==Inf] = NA
  MaxValGenMenDst = max(NanPltGenMenDst, na.rm=TRUE)
  MaxPltGenMenDst = NanPltGenMenDst
  MaxPltGenMenDst[is.na(MaxPltGenMenDst)] = MaxValGenMenDst
  
  ### HEATMAP TITLES #############################################################
  
  if(TypInt=="EndKgo"){
    TypIntTxt = "KEGG-based"
  }
  if(TypInt=="EndPin"){
    TypIntTxt = "PPI-based"
  }
  if(TypInt=="EndC90"){
    TypIntTxt = "COX.90-based"
  }
  if(TypInt=="EndC95"){
    TypIntTxt = "COX.95-based"
  }
  
  DisTxt = paste("Diseases\n",length(DisGen)," Genes",sep="")
  AgeTxt = paste("GenAge.Hum\n",length(AgeGen)," Genes",sep="")
  ModTxt = paste("GenAge.Mod\n",length(ModGen)," Genes",sep="")
  NgbTxt = paste("Neighbours\n",length(DisNgb)," Genes",sep="")
  
  
  ### HEATMAP ####################################################################
  
  PltAllColFctFct = factor(PltAllColFct, levels = c("Dis", "Mod", "Age", "Ngb"))
  
  cn = row.names(UseAllFrm)
  FntSiz = 8
  TxtDeg = 45
  PltAllFrm = as.matrix(UseAllFrm + runif(length(UseAllFrm), min = 0, max = 1e-10))
  
  Htm = Heatmap(t(UseAllFrm), name = "Distance", 
                column_title = paste("Distance to ARCs (",TypIntTxt,")",sep=""), 
                
                row_title = c(ModTxt, AgeTxt, DisTxt, NgbTxt),
                show_row_names = FALSE, show_column_names = FALSE,
                row_title_gp = gpar(fontsize = 9),
                row_split = PltAllColFctFct,
                border = TRUE,
                col = AllDstFltCol,
                
                bottom_annotation = HeatmapAnnotation(
                  text = anno_text(row.names(UseAllFrm), rot = TxtDeg, gp = gpar(fontsize = FntSiz), offset = unit(1, "npc"), just = "right"),
                  annotation_height = max_text_width(AllPhnFct) * sin(TxtDeg*pi/180) * FntSiz/11
                ),
                
                heatmap_legend_param = list(
                  title_gp = gpar(fontsize = 9),   
                  labels_gp = gpar(fontsize = 8),  
                  at = SrtLegFcn(AllDstFltCol) %>% as.numeric() %>% sort() %>% rev() %>% as.character(),
                  #labels = TruLeg(AllDstFltCol) %>% as.numeric() %>% sort() %>% as.character() %>% DirIndFcn() 
                  labels = MaxDifFcn(AllDstFltCol, MaxVal) %>% as.numeric() %>% sort() %>% as.character() %>% DirIndFcn() 
                )
                
                
                
                
  ) +
  
  
  Heatmap(PltGenNg1,
          width = unit(5, "mm"),
          name = "iARC-Interactions",
          border=TRUE,
          col = circlize::colorRamp2( c( min(PltGenNg1), ( max(PltGenNg1) + min(PltGenNg1) )/3 , max(PltGenNg1) ), c("#1E8AC6" , "white", "#F8696B") ),
          
          heatmap_legend_param = list(
            title_gp = gpar(fontsize = 9),   
            labels_gp = gpar(fontsize = 8)),  
          
          show_row_names = FALSE, show_column_names = FALSE,
          bottom_annotation = HeatmapAnnotation(
            text = anno_text("iARC-Interactions", rot = TxtDeg, gp = gpar(fontsize = 8), offset = unit(1, "npc"), just = "right"),
            annotation_height = max_text_width(cn) * 0.7
          )
  ) + 
    Heatmap(as.numeric(MaxPltGenMenDst),
            width = unit(5, "mm"),
            name = "Mean distance",
            border=TRUE,
            col = circlize::colorRamp2( c(min(MaxPltGenMenDst) , 2 , max(MaxPltGenMenDst) ), c("#F8696B" , "#FED280", "#1E8AC6") ) ,
            heatmap_legend_param = list(
              title_gp = gpar(fontsize = 9),   
              labels_gp = gpar(fontsize = 8)),  
            show_row_names = FALSE, show_column_names = FALSE,
            bottom_annotation = HeatmapAnnotation(
              text = anno_text("Mean distance", rot = TxtDeg, gp = gpar(fontsize = 8), offset = unit(1, "npc"), just = "right"),
              annotation_height = max_text_width(cn) * 0.7
            )
  )
  
  
  # DRAW HEATMAP
  HtmLst[[TypInt]] = draw(Htm, padding = unit(c(-10, 2, 1.5, 1), "mm")) 

  
  
  ### TESTS ######################################################################
  
  GenTypNew = AllColFct
  
  GenAgeDst = AgeDst
  GenAgeDst[GenAgeDst==Inf] = NA
  GenAgeDst = GenAgeDst %>% as.numeric()
  
  GenModDst = ModDst
  GenModDst[GenModDst==Inf] = NA
  GenModDst = GenModDst %>% as.numeric()
  
  PhnDst = DisDst
  PhnDst[PhnDst==Inf] = NA
  PhnDst = PhnDst %>% as.numeric()
  
  GenAvg = abs(GenMen - max(GenMen,na.rm=TRUE))
  GenAge = ifelse(GenAgeDst==0, 1, 0)
  GenMod = ifelse(GenModDst==0, 1, 0)
  GenDis = 1*(PhnDst==0)
  GenNgb = ifelse(GenNg1>=1, 1, 0)
  
  
  # CREATE COMPARISON FRAME
  GenPhn = GenDis
  PhnGen = GenDis
  NewCmpFrm = data.frame(Gen=names(GenNg1), 
                         AgeDst=GenAgeDst, ModDst=GenModDst, PhnDst, 
                         GenNg1, 
                         GenAge, GenMod, GenDis, 
                         GenNgb, 
                         GenAvg, GenMen, TypInt)
  
  # MERGE ALL
  sGenTypCmpFrm = merge(GenTypFrm, NewCmpFrm)
  
  if(x == 1){
    GenTypCmpFrm = sGenTypCmpFrm
  } else{
    GenTypCmpFrm = rbind(GenTypCmpFrm, sGenTypCmpFrm)
  }

}

GenTypCmpFrm$TypInt %>% table()

saveRDS(GenTypCmpFrm,"C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenTypCmpFrm.rds")
saveRDS(HtmLst,"C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/HtmLst.rds")

################################################################################
# FUNCTIONS
################################################################################

# COLOR FOR CATEGORICAL ARRAY
CatFrmColFcn = function(GenCatFrm, Pal="Pastel 1"){
  
  # UNIQUE DISTANCES AND LENGTH
  GenCatTbl = GenCatFrm$Cat %>% unique() %>% sort() %>% as.character()
  GenCatLen = length(GenCatTbl)
  
  # DEFINE PALETTE
  GenColArr = hcl.colors(GenCatLen, palette = Pal)
  
  # PALETTE-DST TABLE
  GenColFrm = data.frame('Col'=GenColArr, 'Cat'=GenCatTbl)
  row.names(GenColFrm) = GenColFrm$Cat 
  GenColNamArr = GenColFrm$Col
  names(GenColNamArr) = GenColFrm$Cat %>% as.character()
  
  # COL-DST FRAME
  GenCatColFrm = GenColFrm[as.character(GenCatFrm$Cat),]
  
  # COLOR ARRAY WITH DISTANCE NAMES
  GenCatFltCol = GenCatColFrm$Col
  names(GenCatFltCol) = GenCatColFrm$Cat %>% as.vector() %>% as.character()
  
  return(GenCatFltCol)
  
}




# COLOR FOR CATEGORICAL ARRAY
CatArrColFcn = function(GenCatArr, Pal="Pastel 1"){
  
  # UNIQUE DISTANCES AND LENGTH
  GenCatTbl = GenCatArr %>% unique() %>% sort() %>% as.character()
  GenCatLen = length(GenCatTbl)
  
  # DEFINE PALETTE
  GenColArr = hcl.colors(GenCatLen, palette = Pal)
  
  # PALETTE-DST TABLE
  GenColFrm = data.frame('Col'=GenColArr, 'Cat'=GenCatTbl)
  row.names(GenColFrm) = GenColFrm$Cat
  GenColNamArr = GenColFrm$Col
  names(GenColNamArr) = GenColFrm$Cat %>% as.character()
  
  # COL-DST FRAME
  GenCatColFrm = GenColFrm[as.character(GenCatArr),]
  
  # COLOR ARRAY WITH DISTANCE NAMES
  GenCatFltCol = GenCatColFrm$Col
  names(GenCatFltCol) = GenCatColFrm$Cat %>% as.vector() %>% as.character()
  
  return(GenCatFltCol)
  
}


CatHtmColFcn = function(AllFrm, Pal1 = c("green", "white"), Pal2 = c("white", "red"), MidNum = 3){
  
  MidNum = MidNum + 1
  
  # DEFINE PALETTE
  ColPal1 <- colorRampPalette(Pal1)
  ColPal2 <- colorRampPalette(Pal2)
  
  # CREATE DISTANCE TABLE
  AllDstTbl = AllFrm %>% as.matrix() %>% as.vector() %>% unique() %>% sort()
  AllDstLen = length(AllDstTbl)
  
  # CREATE COLOR FOR TABLE
  if(AllDstTbl[1] == -1){
    DstColArr = c( "#000000", ColPal1(AllDstLen-MidNum) , ColPal2(MidNum+1)[2:(MidNum) ])
  } else if(AllDstTbl[VecDstLen] == Inf){
    DstColArr = c( ColPal1(AllDstLen-MidNum) , ColPal2(MidNum+1)[2:(MidNum), "#000000" ])
  } else{
    DstColArr = c( ColPal1(AllDstLen-MidNum) , ColPal2(MidNum+1+1)[2:(MidNum+1)])
  }
  
  # PALETTE-DST TABLE
  DstColFrm = data.frame('Col'=DstColArr, 'Dst'=AllDstTbl)
  row.names(DstColFrm) = DstColFrm$Dst %>% round()
  DstColNamArr = DstColFrm$Col
  names(DstColNamArr) = DstColFrm$Dst %>% round() %>% as.character()
  
  # FLATTERN DATAFRAME
  AllFltCol = AllFrm %>% as.matrix() %>% as.vector() %>% round() %>% as.character()
  
  # COL-DST FRAME
  AllFltColFrm = DstColFrm[AllFltCol,]
  
  # COLOR ARRAY WITH DISTANCE NAMES
  AllDstFltCol = AllFltColFrm$Col
  names(AllDstFltCol) = AllFltColFrm$Dst %>% round() %>% as.vector() %>% as.character()
  
  return(AllDstFltCol)
  
}


SrtLegFcn = function(DstFltCol){
  SrtNam = DstFltCol %>% names() %>% unique() %>% sort() %>% rev()
  return(SrtNam)
}


DirIndFcn = function(Dst){
  Dst[Dst=="0"] = "0 (Direct association)"
  Dst[Dst=="1"] = "1 (iARC-Interaction)"
  return(Dst) 
}


MaxDifFcn = function(Dst,MaxVal){
  SrtDst = Dst %>% names() %>% unique() %>% sort() %>% rev()
  SrtDstNum = SrtDst %>% as.numeric()
  SrtDstNum[SrtDstNum==-1]=Inf
  SrtDstNumDif = abs(SrtDstNum-MaxVal) %>% round() %>% as.character()
  return(SrtDstNumDif)
}
