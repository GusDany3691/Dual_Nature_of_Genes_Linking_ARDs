##########################################################################################################
# GRAPH ##################################################################################################
##########################################################################################################
#Gen = Gene
#Frm = Frame
#Grp = Graph
#Lst = List
#Nod = Node
#Bgd = Biogrid
#Hum = Human
#Ord = Order
#Ngh = Neighbour
#Nco = Number of Connections
#Prw = Pairwise
#Clt = Cluster

library(rlang)
library(dplyr)
library(igraph)
library(ggrepel)
library(data.table)

####################################################################################################
### START ##########################################################################################
####################################################################################################

# OVERALL PARAMETERS ############################################################################### 

#Mpa_Frm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnFrm_g10K_NoMhc.rds")
Mpa_Frm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnAgeFrm.rds") #%>% filter(!(Phn %in% c("HumAge", "ModAge")))

AgnPhnHrc_Frm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Showcase/PhnHrc_FrmEnd.csv")

IntArr = c('EndPin','EndC90','EndC95', 'EndKgo')
IntLen = length(IntArr)

i = 4

#for(AllCnt in 1:IntLen){
for(i in 1:IntLen){

  IntQry = IntArr[i]
  
  ###########################################################################################################
  #### CREATE GRAPH FROM NETWORK DATASETS 
  ###########################################################################################################
  
  
  #readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/KEGG/GenGenKeg.rds")
  
  if(IntQry == "EndPin"){

    IntFrm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Networks/BIOGRID-MV-Physical-4.4.204.tab3.txt", sep="\t") # 04/12 - 164MB
    
    HumIntFrm = IntFrm%>%dplyr::filter(Organism.ID.Interactor.A==9606,Organism.ID.Interactor.B==9606)
    
    WrkGrp = graph.data.frame(data.frame("Nod_A" = HumIntFrm$Official.Symbol.Interactor.A, 
                                         "Nod_B" = HumIntFrm$Official.Symbol.Interactor.B))
    WrkGrp <- delete_vertices(WrkGrp, V(WrkGrp)[name == "NA"])
    WrkGrp = simplify(WrkGrp)
    WrkGrp = as.undirected(WrkGrp)
    
    #WrkGlb_Gen = V(WrkGrp)$name
    #WrkMpa_Gen = intersect(WrkGlb_Gen,Mpa_Frm$Gen)
    
    WrkFrm = as.data.frame(get.edgelist(WrkGrp))# as_data_frame(Wrk)
    colnames(WrkFrm) = c("Nod_A", "Nod_B")
    
  } 
  if(IntQry == "EndC90"){
    
    IntFrm = read.table("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/CoxPrwFrm90.csv")
    
    AllWrkGen = unique(c(Bgd$row, Bgd$col))
    
    DplEnsNam = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/DplEnsNam90.rds")
    
    EnsNam = DplEnsNam[!duplicated(DplEnsNam$ensembl_gene_id),]
    row.names(EnsNam) = EnsNam$ensembl_gene_id
    
    WrkGen_A = EnsNam[IntFrm$row,]$hgnc_symbol
    WrkGen_B = EnsNam[IntFrm$col,]$hgnc_symbol
    
    PreWrkFrm = data.frame(WrkGen_A, WrkGen_B)
    
    PreWrkFrm = PreWrkFrm[PreWrkFrm$WrkGen_A != '', ]
    PreWrkFrm = PreWrkFrm[PreWrkFrm$WrkGen_B != '', ]
    PreWrkFrm = PreWrkFrm[!is.na(PreWrkFrm$WrkGen_A), ]
    PreWrkFrm = PreWrkFrm[!is.na(PreWrkFrm$WrkGen_B), ]

    
    WrkGrp = graph.data.frame(data.frame("Nod_A" = PreWrkFrm$WrkGen_A, 
                                         "Nod_B" = PreWrkFrm$WrkGen_B))
    WrkGrp <- delete_vertices(WrkGrp, V(WrkGrp)[name == "NA"])
    WrkGrp = simplify(WrkGrp)
    WrkGrp = as.undirected(WrkGrp)
    
    WrkGlb_Gen = V(WrkGrp)$name
    WrkMpa_Gen = intersect(WrkGlb_Gen,Mpa_Frm$Gen)
    
    WrkFrm = as.data.frame(get.edgelist(WrkGrp))# as_data_frame(Wrk)
    colnames(WrkFrm) = c("Nod_A", "Nod_B")
    
    dim(WrkFrm)
    
    UnqGen = c(WrkFrm$Nod_A, WrkFrm$Nod_B) %>% unique()
    
    length(UnqGen)
    
  }
  if(IntQry == "EndC95"){
    
    #IntFrm = read.table("G:/UK_Biobank/Data/Retrieved/Co-expression/CoxPrwFrm95.csv")
    IntFrm = read.table("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/CoxPrwFrm95.csv")
    
    AllWrkGen = unique(c(IntFrm$row, IntFrm$col))
    
    DplEnsNam = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/DplEnsNam95.rds")
    
    EnsNam = DplEnsNam[!duplicated(DplEnsNam$ensembl_gene_id),]
    row.names(EnsNam) = EnsNam$ensembl_gene_id
    
    WrkGen_A = EnsNam[IntFrm$row,]$hgnc_symbol
    WrkGen_B = EnsNam[IntFrm$col,]$hgnc_symbol
    
    PreWrkFrm = data.frame(WrkGen_A, WrkGen_B)
    
    PreWrkFrm = PreWrkFrm[PreWrkFrm$WrkGen_A != '', ]
    PreWrkFrm = PreWrkFrm[PreWrkFrm$WrkGen_B != '', ]
    PreWrkFrm = PreWrkFrm[!is.na(PreWrkFrm$WrkGen_A), ]
    PreWrkFrm = PreWrkFrm[!is.na(PreWrkFrm$WrkGen_B), ]
    
    
    WrkGrp = graph.data.frame(data.frame("Nod_A" = PreWrkFrm$WrkGen_A, 
                                         "Nod_B" = PreWrkFrm$WrkGen_B))
    WrkGrp <- delete_vertices(WrkGrp, V(WrkGrp)[name == "NA"])
    WrkGrp = simplify(WrkGrp)
    WrkGrp = as.undirected(WrkGrp)
    
    WrkGlb_Gen = V(WrkGrp)$name
    WrkMpa_Gen = intersect(WrkGlb_Gen,Mpa_Frm$Gen)
    
    WrkFrm = as.data.frame(get.edgelist(WrkGrp))# as_data_frame(Wrk)
    colnames(WrkFrm) = c("Nod_A", "Nod_B")
  }
  
  
  if(IntQry == "EndKgo"){
    
    IntFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/KEGG/GenGenKeg.rds") 
    colnames(IntFrm) = c("Nod_A", "Nod_B")
    
    WrkGrp = graph.data.frame(data.frame("Nod_A" = IntFrm$Nod_A, 
                                         "Nod_B" = IntFrm$Nod_B))
    WrkGrp <- delete_vertices(WrkGrp, V(WrkGrp)[name == "NA"])
    
    WrkGrp = simplify(WrkGrp)
    WrkGrp = as.undirected(WrkGrp)
    
    #WrkGrp
    
    #WrkGlb_Gen = V(WrkGrp)$name
    #WrkMpa_Gen = intersect(WrkGlb_Gen,Mpa_Frm$Gen)
    
    WrkFrm = as.data.frame(get.edgelist(WrkGrp))# as_data_frame(Wrk)
    colnames(WrkFrm) = c("Nod_A", "Nod_B")
    
    
  } 
  
  
  ##############################################################################
  # INTEGRATING NETWORK AND GENE-PHENOTYPE
  ##############################################################################
  
  NodEnd = list()
  
  NodEnd$All = c(WrkFrm$Nod_A, WrkFrm$Nod_B, Mpa_Frm$Gen, Mpa_Frm$Phn) %>% unique()
  NodEnd$Phn = Mpa_Frm$Phn
  NodEnd$Gen = c(WrkFrm$Nod_A, WrkFrm$Nod_B, Mpa_Frm$Gen) %>% unique()
  NodEnd$FulGen = NodEnd$Gen
  
  
  ### EDGES ##################################################################################################
  
  EdgEnd = list()
  
  Mpa_EdgFrm = Mpa_Frm
  colnames(Mpa_EdgFrm) = c("Nod_A", "Nod_B")
  
  EdgEnd$GlbMpa = rbind(Mpa_EdgFrm,WrkFrm)
  EdgEnd$GenPhn = Mpa_Frm
  EdgEnd$GlbPin = WrkFrm
  
  
  
  'NA' %in% WrkFrm$Nod_A
  'NA' %in% WrkFrm$Nod_B
  'NA' %in% c(WrkFrm$Nod_A, WrkFrm$Nod_B)
  
  ### MORE ON NODES ##########################################################################################
  
  NodEnd$GenPhn = EdgEnd$GenPhn$Gen %>% unique()
  
  NodEnd$GenPin = c(EdgEnd$GlbPin$Nod_A, EdgEnd$GlbPin$Nod_B) %>% unique()
  
  NodEnd$GlbPin = NodEnd$GenPin
  NodEnd$GlbGen = c(NodEnd$GenPhn, NodEnd$GenPin) %>% unique()
  
  #NodEnd$GenPlt = EdgEnd$GenPhn$Gen %>% table() %>% as.data.frame() %>% filter(Freq > 1) %>% PulFcn('.') %>% as.character()
  NodEnd$GenPlt = EdgEnd$GenPhn$Gen %>% table() %>% as.data.frame() %>% filter(Freq > 1) %>% pull('.') %>% as.character()
  
  NodEnd$Gen_UncPhn = setdiff(NodEnd$FulGen, NodEnd$GenPhn)
  NodEnd$Gen_UncPin = setdiff(NodEnd$FulGen, NodEnd$GenPin)
  NodEnd$GenPinLoc_All = NodEnd$GenPin[NodEnd$GenPin %in% NodEnd$GenPhn]
  
  NodEnd$GenMpa = c(NodEnd$GenPhn, NodEnd$GlbPin) %>% unique()
  NodEnd$GenDrgMpa = c(NodEnd$GenMpa, NodEnd$GenDrg) %>% unique()
  
  ### MORE EDGES ##############################################################################################
  
  EdgEnd$LocPin = EdgEnd$GlbPin %>% filter( (Nod_A %in% NodEnd$GenPhn) & (Nod_B %in% NodEnd$GenPhn) ) %>% unique()   #ORIGINAL
  EdgEnd$LocMpa = EdgEnd$GlbMpa %>% filter( !(Nod_A %in% NodEnd$Gen_UncPhn) & !(Nod_B %in% NodEnd$Gen_UncPhn) ) %>% unique()
  
  #ONE MORE NODE
  NodEnd$GenPinLoc_Inn = c(EdgEnd$LocPin$Nod_A, EdgEnd$LocPin$Nod_B) %>% unique() 
  
  ### GRAPHS ##################################################################################################
  
  GrpEnd = list()
  
  SetEdg = c(EdgEnd$GlbMpa$Nod_A, EdgEnd$GlbMpa$Nod_B) %>% unique()
  SetNod = NodEnd$All %>% unique()
  
  GrpEnd$FulGlbMpa = EdgEnd$GlbMpa %>% graph_from_data_frame(directed = FALSE, vertices = NodEnd$All) 
  
  GrpEnd$FulGlbGen = EdgEnd$GlbPin %>% graph_from_data_frame(directed = FALSE, vertices = data.frame(NodEnd$Gen))

  GrpEnd$GlbMpa = EdgEnd$GlbMpa %>% graph_from_data_frame(directed = FALSE) 
  GrpEnd$GlbGen = EdgEnd$GlbPin %>% graph_from_data_frame(directed = FALSE)
  GrpEnd$GlbGenAll = EdgEnd$GlbPin %>% graph_from_data_frame(directed = FALSE, vertices = NodEnd$Gen)
  
  GrpEnd$LocMpa = EdgEnd$LocMpa %>% graph_from_data_frame(directed = FALSE)
  GrpEnd$LocGen = EdgEnd$LocPin %>% graph_from_data_frame(directed = FALSE)
  GrpEnd$LocGenAll = EdgEnd$LocPin %>% graph_from_data_frame(directed = FALSE, vertices = NodEnd$GenPhn)
  
  #sum(is.na(V(GrpEnd$FulGlbMpa)$name))
  #sum(V(GrpEnd$FulGlbMpa)$name=='NA')
  
  ##############################################################################
  # SAVING 
  ##############################################################################
  
  saveRDS(NodEnd,paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",IntQry,"/EndNod.rds", sep=''))
  saveRDS(EdgEnd,paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",IntQry,"/EndEdg.rds", sep=''))
  saveRDS(GrpEnd,paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",IntQry,"/EndGrp.rds", sep=''))
 
}
