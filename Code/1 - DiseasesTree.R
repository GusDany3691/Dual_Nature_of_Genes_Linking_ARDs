library(igraph)
library(tidyverse)

# CODE DICTIONARY ##############################################################
# Lst = List
# Hrc = Hierarchy
# Vec = Vector
# Prn= Parent
# Src = Source
# Anc = Ancestor
# Lup = Loop
# Rut = Root
# Mng = Meaning
# Vrt = Vertice
# Dpth = Depth
# Acl = Ageing Cluster
# Mdl = Modules
# Inc = Inclusive
# Nti = Not Inclusive
# Cnt = Count

################################################################################


DisCat = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Showcase/PhnHrc_FrmOrg.csv")

GrpDis = graph_from_data_frame(DisCat[,c("Node", "Parent")])

NodVec = V(GrpDis)$name

HrcFrm = data.frame()
for(i in 1:length(NodVec)){
  print(i)
  LupNod = NodVec[i]
  if (LupNod != "Top"){
    HrcLst = AgnAscDis_Fcn(LupNod, GrpDis, DisCat)
    sHrcFrm = HrcLst$NodFrm
    HrcFrm = rbind(HrcFrm, sHrcFrm[nrow(sHrcFrm),])
  }
}

HrcFrm = HrcFrm %>% dplyr::rename(AclSpc = spcAcl, PrntSpc = Parent_Specific)
HrcFrm = unique(HrcFrm)
HrcFrm$Type %>% unique()
HrcFrm$X = NULL

write.csv(HrcFrm, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Showcase/PhnHrc_FrmEnd.csv")


#################################################################################################
### FUNCTIONS ###################################################################################
#################################################################################################


AgnAscDis_Fcn <- function(LupNod, GrpDis, DisCat){
  
  ShrtPth = shortest_paths(GrpDis,
                                  from = LupNod,
                                  to = "Top",
                                  mode = c("out", "all", "in")
  )
  
  AllAncNod = as_ids(ShrtPth$vpath[[1]])
  
  if(!(setequal(AllAncNod,"Top"))){
    RutNodIdx = length(AllAncNod) - 1
    AncFrm = DisCat %>% filter(Node %in% AllAncNod) %>% dplyr::rename(Cod = Coding, Mng = Meaning, Sel = Selectable, Nod = Node, Prnt  = Parent, Acl = Ageing, spcAcl = Specific_Ageing)
    AncFrm$RutMng = DisCat %>% filter(Node == AllAncNod[RutNodIdx]) %>% pull("Meaning")
    AncFrm$RutNod = DisCat %>% filter(Node == AllAncNod[RutNodIdx]) %>% pull("Node")
  
    AllAncNod = rev(AllAncNod)
    AllAncLgt = length(AllAncNod) - 1
    OrdFrm = data.frame("Nod" = AllAncNod, "Dpth" = c(0:AllAncLgt))
    
    NodFrm = merge(AncFrm, OrdFrm, by = "Nod") %>% arrange(Dpth)
    
    #AllAncNod = rev(AllAncNod)
    AgnAncNod = list()
    AgnAncNod$All = AgnAncNod_Fcn(NodFrm, AllAncNod, Mdl = c(1,2))
    AgnAncNod$Cl1 = AgnAncNod_Fcn(NodFrm, AllAncNod, Mdl = c(1))
    AgnAncNod$Cl2 = AgnAncNod_Fcn(NodFrm, AllAncNod, Mdl = c(2))
    
    AllAncMng = c("Top" , NodFrm$Mng)
    
    AgnAncMng = list()
    AgnAncMng$All = AgnAncMng_Fcn(NodFrm, AllAncMng, Mdl = c(1,2))
    AgnAncMng$Cl1 = AgnAncMng_Fcn(NodFrm, AllAncMng, Mdl = c(1))
    AgnAncMng$Cl2 = AgnAncMng_Fcn(NodFrm, AllAncMng, Mdl = c(2))
    
    IncAgnAncAcl = c("Top",NodFrm$Acl)
    AgnAncAcl = IncAgnAncAcl[1:(length(IncAgnAncAcl)-1)]
    
    IncPthAgnAclCnt = list()
    IncPthAgnAclCnt$All = sum(IncAgnAncAcl %in% c(1,2))
    IncPthAgnAclCnt$Cl1 = sum(IncAgnAncAcl %in% c(1))
    IncPthAgnAclCnt$Cl2 = sum(IncAgnAncAcl %in% c(2))
    
    PthAgnAclCnt = list()
    PthAgnAclCnt$All = sum(AgnAncAcl %in% c(1,2))
    PthAgnAclCnt$Cl1 = sum(AgnAncAcl %in% c(1))
    PthAgnAclCnt$Cl2 = sum(AgnAncAcl %in% c(2))
    
  
  #Nodes_List = list("Nodes_Frame"=Nodes_Frame, "Ageing_Ascendants"=Ageing_Asendants, "All_Ascendants"=All_Asendants)
  
  NodLst = list("NodFrm"=NodFrm, "AgnAncNod"=AgnAncNod, "AllAncNod"=AllAncNod, "MskAllAncNod"=AgnAncNod, 
                "AllAncMng" = AllAncMng, "AgnAncMng" = AgnAncMng, "AgnAncAcl" = AgnAncAcl,   # Remove MskAgnAncMng
                "IncAgnAncAcl" = IncAgnAncAcl, "PthAgnAclCnt" = PthAgnAclCnt, "IncPthAgnAclCnt" = IncPthAgnAclCnt)
  }
  return(NodLst)
}



AgnAncNod_Fcn = function(NodFrm, AllAncNod, Mdl){  
  AgnAncNod_Lst = list()
  AgnAncNod = NodFrm %>% filter(Acl %in% Mdl) %>% pull("Nod") %>% as.character() #%>% rev()
  MskAgnAncNod = ifelse(AllAncNod %in% AgnAncNod,AllAncNod,"X")
  MskAgnAncNod[1] = "Top"
  AgnAncNod_Lst$AgnAncNod = AgnAncNod
  AgnAncNod_Lst$MskAgnAncNod = MskAgnAncNod
  return(AgnAncNod_Lst)
}


AgnAncMng_Fcn = function(NodFrm, AllAncMng, Mdl){  
  AgnAncMng_Lst = list()
  AgnAncMng = NodFrm %>% filter(Acl %in% Mdl) %>% pull("Mng") %>% as.character()
  AgnAncMng = c("Top", AgnAncMng)
  MskAgnAncMng = ifelse(AllAncMng %in% AgnAncMng,AllAncMng,"X")
  MskAgnAncMng[1] = "Top"
  AgnAncMng_Lst$AgnAncMng = AgnAncMng
  AgnAncMng_Lst$MskAgnAncMng = MskAgnAncMng
return(AgnAncMng_Lst)
}


Agieng_Asendants_Participant <- function(Nodes, Dieases_Graph){
  
  N_Nodes = length(Nodes)
  Ageing_Nodes_List = list()
  for(i in 1:N_Nodes){
    Node = Nodes[i]
    Nodes_List = Node_Asendants(Node, Dieases_Graph)
    Ageing_Nodes_List[[i]] = Nodes_List$Ageing_Nodes
  }
  Agein_Related_Nodes = unique(unlist(Ageing_Nodes_List))
  return(Agein_Related_Nodes)
}



