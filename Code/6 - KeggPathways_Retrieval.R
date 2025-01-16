library("KEGGlincs")
library("KEGGgraph")
library("StarBioTrek")
library("dplyr")
library("limma")
library("stringr")
library("plyr")
library('org.Hs.eg.db')

###  DOWNLOAD LIST OF GENES-RELATED KEGG PATHWAYS ##############################                                                                                                                                            
GenCodKegFrm <- getGeneKEGGLinks(species="hsa")                                                                                 
GenCodKegFrm$Symbol <-  as.character(mapIds(org.Hs.eg.db, GenCodKegFrm$GeneID, 'SYMBOL', 'ENTREZID'))#mapIds(org.Hs.eg.db, tab$GeneID, column="SYMBOL", keytype="ENTREZID")
GenCodKegFrm$PathwayID = str_remove(GenCodKegFrm$PathwayID, "path:")
GenCodKegFrm$GeneID =  paste("hsa:", GenCodKegFrm$GeneID ,sep="")
head(GenCodKegFrm)

#GenCodKegFrm %>% saveRDS( 'C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/GenCodKegFrm.rds' )
GenCodKegFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/KEGG/GenCodKegFrm.rds")


### RETRIEVE KEGG PATHWAYS AND CREATE FRAME ####################################

KegArr = GenCodKegFrm$PathwayID %>% unique()

KegFrm = data.frame()
KegLen = length(KegArr)
i=1
tmp<- tempfile()
for(i in 1:KegLen){
  print(i)
  KegCod = KegArr[i]#"hsa04150"
  KegUrl = KEGGgraph::retrieveKGML(KegCod,"hsa",tmp)
  DstFil = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/KEGG/",KegCod,".xml",sep="")
  download.file(url=KegUrl,destfile=DstFil)
  sKegFrm = parseKGML2DataFrame(file=DstFil,expandGenes=T)
  if(nrow(sKegFrm)>0){
    sKegFrm$pathway = KegCod
    KegFrm = rbind(KegFrm, sKegFrm)
  }     
}

#KegFrm %>% saveRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/KEGG/KegFrm.rds")

### CREATE FINAL VERSION OF KEGG FRAME WITH GENE CODES #########################

KegFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/KegFrm.rds") %>% as.data.table()

Keg_from = GenCodKegFrm %>% SelFcn(c('GeneID', 'Symbol')) %>% RenNamFcn(OldCol=c('GeneID', 'Symbol'),NewCol=c("from", "Gen_from")) 
Keg_to   = GenCodKegFrm %>% SelFcn(c('GeneID', 'Symbol')) %>% RenNamFcn(OldCol=c('GeneID', 'Symbol'),NewCol=c("to", "Gen_to")) 

Frm1 = inner_join(KegFrm, Keg_from, by="from") %>% unique()
Frm = inner_join(Frm1, Keg_to, by="to") %>% unique()

GenGenKeg = Frm %>% dplyr::select(Gen_from, Gen_to) %>% unique()

GenGenKeg %>% saveRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/KEGG/GenGenKeg.rds")


################################################################################
# FUNCTIONS
################################################################################

RenNamFcn = function(OrgFrm, OldCol, NewCol){
  ColNam = colnames(OrgFrm)
  LenRep = length(OldCol)
  for(i in 1:LenRep){
    colnames(OrgFrm)[ColNam == OldCol[i]] = NewCol[i]#TrnFrm
  }
  return(OrgFrm)
}

SelFcn = function(Frm, ColArr){
  return(Frm[ColArr])
}

PulFcn = function(Frm, Col){
  return(Frm[[Col]])
}