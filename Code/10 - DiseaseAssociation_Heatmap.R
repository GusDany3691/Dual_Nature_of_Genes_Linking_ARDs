library(dplyr)
library(reshape2)
library(dplyr)
library(igraph)
library(rje)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(devtools)

## INITIAL PARAMETERS #########################################################################################

# SAVING PARAMETERS
#TypInt = "EndPin"
#VerTxt = "_V2_DelHmc"

# SAVING LINK AND FURTHER LOADS
#TxtPth = paste("G:/UK_Biobank/Data/Generated/Sorted/MPA/Connectivity_analysis/DLHCP_dAcA_ThrAll/",TypInt,"/",sep="")

#TxtPth = paste("G:/UK_Biobank/Data/Generated/Sorted/MPA/Connectivity_analysis/DLHCP_dAcA_ThrAll/",TypInt,"/",sep="")
NorFrmLst = list()

### RETRIEVE DATA #################################################

GenPhnComFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")
GenPhnAgeFrm = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnAgeFrm.rds')
GenAgeHumArr = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeHum_Genes.rds")
GenAgeModArr = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeMod_Genes.rds")
#GenAgePhnFrm$Phn %>% unique()

#NodAll = readRDS(paste(TxtPth,"EndNod",VerTxt,".rds",sep=""))
#EdgAll = readRDS(paste(TxtPth,"EndEdg",VerTxt,".rds",sep=""))

#EdgAll$GenPhn$Phn %>% unique() %>% sort()
#EdgAll$GenPhn = EdgAll$GenPhn %>% filter(!(Phn %in% c("C1007",  "C1010",  "C1028",  "C1046",  "C1052",  "C1067",  "C1068",  "C1069")))

#AgeModGenExl = read.csv("G:/UK_Biobank2/The Project/Databases/DataSetsHagr/Final/OMA/AgeingFinal.csv")
#AgeModGen = AgeModGenExl$AgeingGenes
#AgeModFrm = data.frame(Gen=AgeModGen, Phn="AgeMod")

PhnArr = GenPhnAgeFrm$Phn %>% unique()
#NodAll$Phn = c(NodAll$Phn, "ModAge") %>% unique()
EdgAll$GenPhn = rbind(EdgAll$GenPhn,AgeModFrm)
#HumAge = EdgAll$GenPhn %>% filter(Phn == "HumAge") %>% pull(Gen)

NodAll$GenPhn = c(NodAll$GenPhn, AgeModGen) %>% unique()

#GenPhnCom = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")
#GenPhnCom = GenPhnCom %>% filter(Mng != 'cancers')
ComArr = GenPhnComFrm$Com %>% unique()
#ComArr = GenPhnCom$Com %>% unique()


# GenAge Humans Association of Communities
YesHumCom = GenPhnComFrm %>% filter(Gen %in% AgeHumGen) %>% pull(Com) %>% unique()
NonHumCom = setdiff(ComArr,YesHumCom)
YesHumComFrm = data.frame(Dst = 0, Com = YesHumCom)
NonHumComFrm = data.frame(Dst = 1, Com = NonHumCom)
HumComFrm = rbind(YesHumComFrm,NonHumComFrm)
row.names(HumComFrm) = HumComFrm$Com


# GenAge Models Association of Communities
YesModCom = GenPhnComFrm %>% filter(Gen %in% AgeModGen) %>% pull(Com) %>% unique()
NonModCom = setdiff(ComArr,YesModCom)
YesModComFrm = data.frame(Dst = 0, Com = YesModCom)
NonModComFrm = data.frame(Dst = 1, Com = NonModCom)
ModComFrm = rbind(YesModComFrm,NonModComFrm)
row.names(ModComFrm) = ModComFrm$Com


#AgnPhnHrc_Frm = read.csv("G:/UK_Biobank/Data/Generated/Sorted/Phenotypes/PhnHrc_Frm.csv")
AgnPhnHrc_Frm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Showcase/PhnHrc_FrmEnd.csv")

#AgnPhnHrc_Frm = AgnPhnHrc_Frm %>% filter(Nod %in% NodAll$Phn, Acl %in% c(1,2))
AgnPhnHrc_Frm = AgnPhnHrc_Frm %>% filter(Nod %in% PhnArr, Acl %in% c(1,2))
#SmlAgnPhnHrc_Frm = AgnPhnHrc_Frm %>% dplyr::select(Nod, Mng, RutMng) %>% dplyr::rename(Phn=Nod)



#GenPhnInf = inner_join(EdgAll$GenPhn, SmlAgnPhnHrc_Frm, by="Phn")
#GenPhnInf = inner_join(GenPhnAgeFrm, SmlAgnPhnHrc_Frm, by="Phn")


#GenAge = EdgAll$GenPhn %>% filter(Phn == "HumAge") %>% pull(Gen)
#GenAgeHum = GenPhnAgeFrm %>% filter(Phn == "HumAge") %>% pull(Gen)
#GenArr = GenPhnInf$Gen %>% unique()

DisGenArr = GenPhnComFrm$Gen %>% unique()
DisGenLen = length(DisGenArr)
#GenComFrm = data.frame()
ComNumFrm = data.frame()

ComArr = AgnPhnHrc_Frm %>% filter(Type=="D") %>% pull(RutMng) %>% unique()
ComAssFrm = data.frame(Cnt=ComArr)
row.names(ComAssFrm) = ComAssFrm$Cnt 
ComAssFrm$Cnt = 0
ComAssFrm = t(ComAssFrm)
sComAssFrm = ComAssFrm
RowComAssFrm = ComAssFrm
for(i in 1:DisGenLen){
  sComAssFrm = RowComAssFrm
  print(i)
  GenQry = DisGenArr[i]#GenArr[i]
  #ComArr = GenPhnInf %>% filter(Gen == GenQry) %>% pull(ComMng)
  ComArr = GenPhnComFrm %>% filter(Gen == GenQry) %>% pull(Com)
  NumPhn = ComArr %>% length()
  NumCom = ComArr %>% unique() %>% length()
  ComTbl = table(ComArr)
  sComAssFrm[,names(ComTbl)] = as.numeric(ComTbl)
  row.names(sComAssFrm) = GenQry
  if(i == 1){
    ComAssFrm = sComAssFrm
  } else{
    ComAssFrm = rbind(ComAssFrm, sComAssFrm)
  }
  sComNumFrm = data.frame(Gen=GenQry, NumPhn, NumCom)
  ComNumFrm = rbind(ComNumFrm, sComNumFrm)
  
}


#ComFrm = data.frame(ComFrm)
#ComFrm$Gen = row.names(ComFrm)
ComAssFrm = data.frame(ComAssFrm)
ComAssFrm$Gen = row.names(ComAssFrm)
ComAssNumFrm = merge(ComAssFrm, ComNumFrm, by="Gen") %>% arrange(desc(NumCom))
#TmpComFrm  = TmpComFrm %>% arrange(desc(NumCom))
#BigComFrm = TmpComFrm[,!(colnames(TmpComFrm) %in% c("dietary.restriction", "human.ageing", "cell.age"))]
#BigComFrm  = TmpComFrm %>% arrange(desc(NumCom))
#ComFrm = BigComFrm

#ComAssNumFrm %>% mutate(GenAgeHumTyp = ifelse(Gen %in% GenAgeHumArr,))


##########################################################################################################
### PLOT #################################################################################################
##########################################################################################################

#PltNgbFrm = readRDS( paste(TxtPth,"PltNgbFrm",VerTxt,".rds",sep="") ) #LOOK FOR THE FILE: Datasets_plotting
#ComAgeNgbAdj = readRDS( paste(TxtPth,"ComAgeNgbAdj",VerTxt,".rds",sep="") ) #LOOK FOR THE FILE: Datasets_plotting

#ComAgeNgbAdj = ComAgeNgbAdj %>% filter(!Com %in% c("gastrointestinal cancer", "genital tract cancer", "breast cancer", "skin cancer"))
#CanComAgeNgbAdj = data.frame(Dst=0, Com="cancer")
#ComAgeNgbAdj = rbind(ComAgeNgbAdj,CanComAgeNgbAdj)


#NumAgeCol = PltNgbFrm$GenAgeDst %>% unique() %>% length()
#AgeNgbArr = PltNgbFrm$GenAgeDst + 1
#AgeNgbIdxArr = AgeNgbArr
#AgeNgbIdxArr[AgeNgbIdxArr == Inf] = NumAgeCol


#GrnDrgScrFrm = readRDS("G:/UK_Biobank/Data/Generated/Sorted/MPA/Connectivity_analysis/DLHCP_dAcA_ThrAll/Pin/Drg/GrnDrgScrFrm.rds")
#MpaGen = GrnDrgScrFrm %>% filter(Mpa == TRUE) %>% pull(Gen) %>% unique()


#MpaRutFrm = BigRutFrm %>% filter(Gen %in% MpaGen, NumRut>=1)
#SmlRutFrm = MpaRutFrm 
#SmlGen = SmlRutFrm$Gen

#AgeDifGen = setdiff(GenAgeHum,SmlGen)
#AgeDifLen = length(AgeDifGen)
#SmlRutCol = colnames(SmlRutFrm)
#SmlRutColLen = length(SmlRutCol)

#AgeDiFrm = matrix(0,AgeDifLen,SmlRutColLen ) %>% as.data.frame()
#colnames(AgeDiFrm) = SmlRutCol
#AgeDiFrm$Gen = AgeDifGen


#QryFrm = SmlRutFrm[,!(colnames(SmlRutFrm) %in% c("NumMng", "NumRut"))]
QryFrm = ComAssNumFrm[,!(colnames(ComAssNumFrm) %in% c("NumPhn", "NumCom"))]
row.names(QryFrm) = QryFrm$Gen
QryFrm$Gen = NULL
QryFrm = QryFrm[rowSums(QryFrm)>0,]
tQryFrm = t(QryFrm)


#MorSmlRutFrm = SmlRutFrm %>% left_join(PltNgbFrm, by="Gen")

# COLUMN-SCORE FRAMES
FrmMng = ComAssNumFrm$NumPhn #MorSmlRutFrm$NumMng
FrmRut = ComAssNumFrm$NumCom #MorSmlRutFrm$NumRut
#FrmAge = #MorSmlRutFrm$GenAgeDst#GenAgeTyp
#FrmAge[is.na(FrmAge)] = 0
#FrmAge_NonInf = FrmAge + 1
#FrmAge_NonInf[FrmAge_NonInf == Inf] = NumAgeCol

MaxMng = max(FrmMng)
MedMng = max(FrmMng)/2
MinMng = min(FrmMng)

MaxRut = max(FrmRut)
MedRut = max(FrmRut)/2
MinRut = min(FrmRut)

MaxFrm = max(QryFrm)
MinFrm = min(QryFrm)

WhtRedPal <- colorRampPalette(c("white", "red"))   # Apply colorRampPalette
WhtBluPal <- colorRampPalette(c("white", "blue"))   # Apply colorRampPalette
WhtBlkPal <- colorRampPalette(c("white", "black"))   # Apply colorRampPalette

WhtGrnPal <- colorRampPalette(c("green", "white", "orange"))   # Apply colorRampPalette

FltQryFrm = reshape(cbind(QryFrm, id=1, time=rownames(QryFrm)),direction="wide", sep="_", new.row.names=0)[-1]
FltQryVal = as.numeric(FltQryFrm)


RutColArr = colorRamp2(c(0, MaxRut), c("white", "blue"))
MngColArr = colorRamp2(c(0, MedMng, MaxMng), c("green", "white", "orange"))
FrmColArr = colorRamp2(c(0, MaxFrm), c("white", "red"))
#AgeColArr = colorRamp2(c(0, NumAgeCol), c("Black", "white"))


FrmCol = WhtRedPal(MaxFrm-MinFrm+1)
RutCol = WhtBluPal(MaxRut-MinRut+1+1)[(1+1):(MaxRut-MinRut+1+1)]
#AgeCol = c("#000000", hcl.colors(NumAgeCol-2, palette = "Pastel 1"), "#FFFFFF")


MngCol = WhtGrnPal(MaxMng-MinMng+1)
FrmColArr = FrmCol[FltQryVal-MinFrm+1]
RutColArr = RutCol[FrmRut-MinRut+1]
MngColArr = MngCol[FrmMng-MinMng+1]
#AgeColArr = AgeCol[FrmAge_NonInf]


names(FrmColArr) = FltQryVal
names(RutColArr) = FrmRut
names(MngColArr) = FrmMng
#names(AgeColArr) = FrmAge

RutLeg = FrmRut %>% unique() %>% sort(decreasing = TRUE)
MngLeg = FrmMng %>% unique() %>% sort(decreasing = TRUE)
FrmLeg = FltQryVal %>% unique() %>% sort(decreasing = TRUE)
#AgeLeg = c("GenAge", "GenAge_Ngb1", "GenAge_Ngb2", "GenAge_Ngb3", "GenAge_Ngb4", "Genge_Ngb5", "GenAge_Unrecahable")

#AgeColArr[]


#FrmAgeTrf = FrmAge
#FrmAgeTrf[FrmAgeTrf == Inf] = 6
#FrmAgeTrf = abs(FrmAgeTrf-6)/6

#NewFrmAge = ifelse(FrmAge==0,"GenAge","Diseases")
#NewAgeColArr = AgeColArr
#names(NewAgeColArr) = NewFrmAge
#NewAgeColArr[NewFrmAge=="GenAge"] = "#FFC5D0" 
#NewAgeColArr[NewFrmAge=="Diseases"] = "#D4D8A7" 


### COLUMNS ANNOTATION ###########################################################################

NewFrmRut = FrmRut
NewRutColArr = RutColArr
NewMngColArr = MngColArr
NewRutLeg = RutLeg
NewMngLeg = MngLeg

NewFrmRut %>% length()
FrmRut %>% length()

NewRutColArr %>% length()
RutColArr %>% length()

NewMngColArr %>% length()
MngColArr %>% length()

NewRutLeg %>% length()
RutLeg %>% length()

NewMngLeg %>% length()
MngLeg %>% length()






ha = HeatmapAnnotation(
  #Total_ARCs = FrmRut,
  Total_ARCs = FrmRut,#data.frame(Mean = PhnMen),
  Total_ARDs = FrmMng,
  
  col = list(Total_ARCs = RutColArr, Total_ARDs = MngColArr),
  
  annotation_height = unit(4, "mm"), 

  annotation_name_gp= gpar(fontsize = 8, fontface = "bold"),
  
  annotation_legend_param = list(


                                 Total_ARCs = list(title = "Total_ARCs",
                                                          at = RutLeg,
                                                          title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                                          border = "black", 
                                                          labels_gp = gpar(col = "black", fontsize = 8), 
                                                          legend_height = unit(6, "cm"),
                                                          direction = "horizontal"
                                                          ),
                                 
                                 Total_ARDs = list(title = "Total_ARDs", 
                                                         at = MngLeg,
                                                         title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                                         border = "black", 
                                                         labels_gp = gpar(col = "black", fontsize = 8), 
                                                         legend_height = unit(6, "cm"),
                                                         legend_direction = "horizontal"
                                                         )#,
                                 
                                 
                                 )
)


### ROW ANNOTATION ##############################################################################

# PREPARE GEN_AGE DST DATA
RowNam = row.names(tQryFrm)

#row.names(ComAgeNgbAdj) = ComAgeNgbAdj$Com


RowNam = RowNam %>% str_replace("[.]", "/") %>% str_replace("eye.", "eye/") %>% str_replace("[.]", " ") %>% str_replace("genital/tract", "genital tract") # %>% str_replace("/cancer", " cancer")
#EndComAgeNgbAdj = ComAgeNgbAdj[RowNam,]

### MODEL AGE ###################################################################################

#ModComFrm = ModComFrm[RowNam,]
#HumComArr = HumComFrm$Dst
#ModComArr = ModComFrm$Dst

HumComArr = HumComFrm[RowNam,]
ModComArr = ModComFrm[RowNam,]

# ADJUST GEN_AGE_ DST COLORS

ModComFrm

#ComDstLen = EndComAgeNgbAdj$Dst %>% unique() %>% length()
#ComDstArr = EndComAgeNgbAdj$Dst
#ComDstColVal = EndComAgeNgbAdj$Dst %>% unique()

#ComDstLen = ModComFrm$Dst %>% unique() %>% length()
#ComDstArr = ModComFrm$Dst

#ComDstColVal = ModComFrm$Dst %>% unique()
#ComDstColNam = hcl.colors(ComDstLen, palette = "Dynamic")#c("#000000", hcl.colors(ComDstLen-1, palette = "Pastel 1"))

ComDstColNam = hcl.colors(2, palette = "Dynamic")#c("#000000", hcl.colors(ComDstLen-1, palette = "Pastel 1"))
#ComDstColArr = ComDstColNam[ComDstArr+1]
#names(ComDstColArr) = ComDstArr



# ADJUST GEN_AGE_ DST COLORS
ComGenSumArr = (tQryFrm > 0) %>% rowSums() %>% as.data.frame()
colnames(ComGenSumArr) = "NumGen"


MaxScr = max(ComGenSumArr$NumGen)
MedScr = MaxScr/2

#ModComArr
#ModPltComDstArr = ifelse(ModComArr == 0, "Ovelapping", "Non.Ovelapping")
#ModComDstColArr = ComDstColNam[ModComArr+1]
#names(ModComDstColArr) = ModComArr
#ModPltComDstColArr = ModComDstColArr
#names(ModPltComDstColArr) = ifelse(ModComArr=="0","Ovelapping","Non.Ovelapping")
ModPltComDstArr = ifelse(ModComArr$Dst == 0, "Ovelapping", "Non.Ovelapping")
ModComDstColArr = ComDstColNam[ModComArr$Dst+1]
names(ModComDstColArr) = ModComArr
ModPltComDstColArr = ModComDstColArr
names(ModPltComDstColArr) = ifelse(ModComArr$Dst=="0","Ovelapping","Non.Ovelapping")


#PltComDstArr = ifelse(ComDstArr == 0, "Ovelapping", "Non.Ovelapping")
#PltComDstColArr = ComDstColArr
#names(PltComDstColArr) = ifelse(names(PltComDstColArr)=="0","Ovelapping","Non.Ovelapping")
HumPltComDstArr = ifelse(HumComArr$Dst == 0, "Ovelapping", "Non.Ovelapping")
HumComDstColArr = ComDstColNam[HumComArr$Dst+1]
names(HumComDstColArr) = HumComArr
HumPltComDstColArr = HumComDstColArr
names(HumPltComDstColArr) = ifelse(HumComArr$Dst=="0","Ovelapping","Non.Ovelapping")



# note how we set the width of this empty annotation
ra = rowAnnotation( Gene_Count = anno_barplot(ComGenSumArr$NumGen),
                    #GenAge.Hum = PltComDstArr,
                    GenAge.Hum = HumPltComDstArr,
                    GenAge.Mod = ModPltComDstArr,
                    
                    #col = list(GenAge.Hum = PltComDstColArr, GenAge.Mod=ModPltComDstColArr),
                    col = list(GenAge.Hum = HumPltComDstColArr, GenAge.Mod=ModPltComDstColArr),
                    
                    annotation_name_gp= gpar(fontsize = 8, fontface = "bold"),
                    show_legend = c(TRUE, TRUE, FALSE),
                    annotation_legend_param = list(
                    GenAge.Hum = list(title = "GenAge (Hum & Mod)",
                                                              title_gp = gpar(fontsize = 8, fontface = "bold"), 
                                                              border = "black", 
                                                              labels_gp = gpar(col = "black", fontsize = 8), 
                                                              legend_height = unit(6, "cm"),
                                                              direction = "horizontal")
          ),
          annotation_name_rot = 90  

)




### HEATHMAP ######################################################################################


#GenAgeTyp=ifelse(colnames(tQryFrm)%in%GenAge,"GenAge","Non-GenAge")

ncol(t(QryFrm))
cn = colnames(QryFrm)
rn = row.names(QryFrm)
LenRow = length(rn)
LenCol = length(cn)

length(FrmColArr) / 9
dim(tQryFrm)

Htm = Heatmap(as.matrix(tQryFrm), name = "Number of diseases", top_annotation = ha, 
              right_annotation = ra,

              col = FrmColArr,
              
              column_title = paste("Genes (",LenRow,")",sep=""), row_title = paste("ARCs (",LenCol,")",sep="" ),
              use_raster = FALSE, # FEDAULT IS TRUE

              show_column_names = FALSE, 
              border="black",

              row_title_gp=gpar(fontsize=10,fontface="bold"),
              column_title_gp=gpar(fontsize=10,fontface="bold"),
              row_names_gp = gpar(fontsize = 8),
              show_row_names=TRUE,
              
              heatmap_legend_param = list(
                title = "ARDs", 
                at = FrmLeg,
                title_gp = gpar(fontsize = 8, fontface = "bold"), 
                border = "black", 
                labels_gp = gpar(col = "black", fontsize = 8), 
                legend_height = unit(6, "cm"),
                
                direction = "horizontal",

                legend_direction = "horizontal"
                )
                                            
              

) 

### PRINT HEATMAP ##############################################################

draw(Htm, 
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom", 
     merge_legend = TRUE, 
     legend_direction = "horizontal")


draw(Htm, 
     show_heatmap_legend = FALSE, 
     show_annotation_legend = FALSE)


draw(Htm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)


draw(Htm,merge_legend = TRUE)

### GEN AGE ASSOCIATIONS #######################################################

GenPhnComFrm %>% select(Gen, Com) %>% unique() %>% filter(Gen %in% GenAgeHumArr) %>% pull(Com) %>% table() %>% as.data.frame()
GenPhnComFrm %>% select(Gen, Com) %>% unique() %>% filter(Gen %in% GenAgeModArr) %>% pull(Com) %>% table() %>% as.data.frame()





