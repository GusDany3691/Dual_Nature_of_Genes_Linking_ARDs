library(igraph)

TypIntArr = c("EndPin", "EndC95", "EndC90", "EndKgo")

### LOADING DATA ###############################################################

GenPhnCom = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")
GenAge = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeHum_Genes.rds")
GenComAge = data.frame("Gen"=GenAge,"Com"="GenAge_Hum")
GenComFrm = data.frame(GenPhnComFrm$Gen, GenPhnComFrm$Com) %>% unique() 
colnames(GenComFrm) = c("Gen", "Com")
GenComFrm = GenComFrm%>% rbind(GenComAge)


AgeModGen = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeMod_Genes.rds")
AgeModFrm = data.frame(Gen=AgeModGen, Com="GenAge_Mod")
GenComFrm = GenComFrm %>% rbind(AgeModFrm)

HumGen = GenComFrm %>% filter(Com == "GenAge_Hum") %>% pull(Gen)
ModGen = GenComFrm %>% filter(Com == "GenAge_Mod") %>% pull(Gen)
ImmGen = GenComFrm %>% filter(Com == "immunological/systemic disorders") %>% pull(Gen)
ArcGen = GenComFrm %>% filter(!(Com %in% c("GenAge_Hum", "GenAge_Mod", "cancer"))) %>% pull(Gen)
AgeGen = intersect(HumGen, ModGen)
ArcAgeGen = GenComFrm$Gen %>% unique()


## NEXT #################################################

TypIntLen = length(TypIntArr)
CntFrm = data.frame()
for(x in 1:TypIntLen){
  
  TypInt = TypIntArr[x]
  TxtPth = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/",sep="")
  NorFrmLst = list()
  
  ### RETRIEVE DATA #################################################
  
  Grp = readRDS(paste(TxtPth,"EndGrp.rds",sep=""))
  
  IntFrm = calculate_overlaps(Grp$GlbGen, ArcGen)
  colnames(IntFrm) = c("Gen", "Deg", "Ovl", "Rat")
  IntGrp = IntFrm$Rat
  names(IntGrp) = IntFrm$Gen
  
  ### COMPUTE CENTALITIES ###########################################
  paste(x,".Deg",sep="") %>% print()
  DegGrp = igraph::degree(Grp$GlbGen) %>% sort() %>% rev()
  
  paste(x,".Btw",sep="") %>% print()
  BtwGrp = igraph::betweenness(Grp$GlbGen) %>% sort() %>% rev()
  
  paste(x,".Cls",sep="") %>% print()
  ClsGrp = igraph::closeness(Grp$GlbGen) %>% sort() %>% rev()
  
  paste(x,".Trn",sep="") %>% print()
  TrnGrp = full_transitivity(Grp$GlbGen) %>% sort() %>% rev()
  
  paste(x,".Eig",sep="") %>% print()
  EigGrp = igraph::eigen_centrality(Grp$GlbGen)$vector %>% sort() %>% rev() 


  ### MAP CENTRALITIES ##############################################
  
  GrpGen = V(Grp$GlbGen)$name
  HumGenGrp = intersect(HumGen, GrpGen)
  ModGenGrp = intersect(ModGen, GrpGen)
  ImmGenGrp = intersect(ImmGen, GrpGen)
  ArcGenGrp = intersect(ArcGen, GrpGen)
  AgeGenGrp = intersect(AgeGen, ModGen)
  OthGenGrp = setdiff(GrpGen, ArcAgeGen)

  # CREATE CENTRALITY FRAMES
  DegHumFrm = data.frame(Com="Hum", Scr=DegGrp[HumGenGrp], Cnt="Deg", Ntw=TypInt)
  DegModFrm = data.frame(Com="Mod", Scr=DegGrp[ModGenGrp], Cnt="Deg", Ntw=TypInt)
  DegImmFrm = data.frame(Com="Imm", Scr=DegGrp[ImmGenGrp], Cnt="Deg", Ntw=TypInt)
  DegArcFrm = data.frame(Com="Arc", Scr=DegGrp[ArcGenGrp], Cnt="Deg", Ntw=TypInt)
  DegAgeFrm = data.frame(Com="Age", Scr=DegGrp[AgeGenGrp], Cnt="Deg", Ntw=TypInt)
  DegOthFrm = data.frame(Com="Oth", Scr=DegGrp[OthGenGrp], Cnt="Deg", Ntw=TypInt)
  
  BtwHumFrm = data.frame(Com="Hum", Scr=BtwGrp[HumGenGrp], Cnt="Btw", Ntw=TypInt)
  BtwModFrm = data.frame(Com="Mod", Scr=BtwGrp[ModGenGrp], Cnt="Btw", Ntw=TypInt)
  BtwImmFrm = data.frame(Com="Imm", Scr=BtwGrp[ImmGenGrp], Cnt="Btw", Ntw=TypInt)
  BtwArcFrm = data.frame(Com="Arc", Scr=BtwGrp[ArcGenGrp], Cnt="Btw", Ntw=TypInt)
  BtwAgeFrm = data.frame(Com="Age", Scr=BtwGrp[AgeGenGrp], Cnt="Btw", Ntw=TypInt)
  BtwOthFrm = data.frame(Com="Oth", Scr=BtwGrp[OthGenGrp], Cnt="Btw", Ntw=TypInt)
  
  ClsHumFrm = data.frame(Com="Hum", Scr=ClsGrp[HumGenGrp], Cnt="Cls", Ntw=TypInt)
  ClsModFrm = data.frame(Com="Mod", Scr=ClsGrp[ModGenGrp], Cnt="Cls", Ntw=TypInt)
  ClsImmFrm = data.frame(Com="Imm", Scr=ClsGrp[ImmGenGrp], Cnt="Cls", Ntw=TypInt)
  ClsArcFrm = data.frame(Com="Arc", Scr=ClsGrp[ArcGenGrp], Cnt="Cls", Ntw=TypInt)
  ClsAgeFrm = data.frame(Com="Age", Scr=ClsGrp[AgeGenGrp], Cnt="Cls", Ntw=TypInt)
  ClsOthFrm = data.frame(Com="Oth", Scr=ClsGrp[OthGenGrp], Cnt="Cls", Ntw=TypInt)
  
  TrnHumFrm = data.frame(Com="Hum", Scr=TrnGrp[HumGenGrp], Cnt="Trn", Ntw=TypInt)
  TrnModFrm = data.frame(Com="Mod", Scr=TrnGrp[ModGenGrp], Cnt="Trn", Ntw=TypInt)
  TrnImmFrm = data.frame(Com="Imm", Scr=TrnGrp[ImmGenGrp], Cnt="Trn", Ntw=TypInt)
  TrnArcFrm = data.frame(Com="Arc", Scr=TrnGrp[ArcGenGrp], Cnt="Trn", Ntw=TypInt)
  TrnAgeFrm = data.frame(Com="Age", Scr=TrnGrp[AgeGenGrp], Cnt="Trn", Ntw=TypInt)
  TrnOthFrm = data.frame(Com="Oth", Scr=TrnGrp[OthGenGrp], Cnt="Trn", Ntw=TypInt) 
  
  EigHumFrm = data.frame(Com="Hum", Scr=EigGrp[HumGenGrp], Cnt="Eig", Ntw=TypInt)
  EigModFrm = data.frame(Com="Mod", Scr=EigGrp[ModGenGrp], Cnt="Eig", Ntw=TypInt)
  EigImmFrm = data.frame(Com="Imm", Scr=EigGrp[ImmGenGrp], Cnt="Eig", Ntw=TypInt)
  EigArcFrm = data.frame(Com="Arc", Scr=EigGrp[ArcGenGrp], Cnt="Eig", Ntw=TypInt)
  EigAgeFrm = data.frame(Com="Age", Scr=EigGrp[AgeGenGrp], Cnt="Eig", Ntw=TypInt)
  EigOthFrm = data.frame(Com="Oth", Scr=EigGrp[OthGenGrp], Cnt="Eig", Ntw=TypInt) 
  
  IntHumFrm = data.frame(Com="Hum", Scr=IntGrp[HumGenGrp], Cnt="Int", Ntw=TypInt)
  IntModFrm = data.frame(Com="Mod", Scr=IntGrp[ModGenGrp], Cnt="Int", Ntw=TypInt)
  IntImmFrm = data.frame(Com="Imm", Scr=IntGrp[ImmGenGrp], Cnt="Int", Ntw=TypInt)
  IntArcFrm = data.frame(Com="Arc", Scr=IntGrp[ArcGenGrp], Cnt="Int", Ntw=TypInt)
  IntAgeFrm = data.frame(Com="Age", Scr=IntGrp[AgeGenGrp], Cnt="Int", Ntw=TypInt)
  IntOthFrm = data.frame(Com="Oth", Scr=IntGrp[OthGenGrp], Cnt="Int", Ntw=TypInt)
  
  
  paste(x,".BndDeg",sep="") %>% print()
  DegFrm = DegHumFrm %>% rbind(DegModFrm) %>% rbind(DegImmFrm) %>% rbind(DegArcFrm) %>% rbind(DegAgeFrm) %>% rbind(DegOthFrm)
  
  paste(x,".BndBtw",sep="") %>% print()
  BtwFrm = BtwHumFrm %>% rbind(BtwModFrm) %>% rbind(BtwImmFrm) %>% rbind(BtwArcFrm) %>% rbind(BtwAgeFrm) %>% rbind(BtwOthFrm)
  
  paste(x,".BndCls",sep="") %>% print()
  ClsFrm = ClsHumFrm %>% rbind(ClsModFrm) %>% rbind(ClsImmFrm) %>% rbind(ClsArcFrm) %>% rbind(ClsAgeFrm) %>% rbind(ClsOthFrm)
  
  paste(x,".BndTrn",sep="") %>% print()
  TrnFrm = TrnHumFrm %>% rbind(TrnModFrm) %>% rbind(TrnImmFrm) %>% rbind(TrnArcFrm) %>% rbind(TrnAgeFrm) %>% rbind(TrnOthFrm)
  
  paste(x,".BndEig",sep="") %>% print()
  EigFrm = EigHumFrm %>% rbind(EigModFrm) %>% rbind(EigImmFrm) %>% rbind(EigArcFrm) %>% rbind(EigAgeFrm) %>% rbind(EigOthFrm)
  
  paste(x,".BndInt",sep="") %>% print()
  IntFrm = IntHumFrm %>% rbind(IntModFrm) %>% rbind(IntImmFrm) %>% rbind(IntArcFrm) %>% rbind(IntAgeFrm) %>% rbind(IntOthFrm)

  paste(x,".sBndAll",sep="") %>% print()
  sCntFrm = DegFrm %>% rbind(BtwFrm) %>% rbind(ClsFrm) %>% rbind(TrnFrm) %>% rbind(EigFrm) %>% rbind(IntFrm)
  
  paste(x,".BndAll",sep="") %>% print()
  CntFrm = rbind(CntFrm, sCntFrm)
  
}


### TABLE ##########################################################


rnd = 2
TrmCntFrm = transform_data(CntFrm)
TrmCntFrm$Btw = round(TrmCntFrm$Btw)
TrmCntFrm$Deg = round(TrmCntFrm$Deg)

TrmCntFrm$Cls = paste(round(TrmCntFrm$Cls*100), "%", sep="")
TrmCntFrm$Eig = paste(round(TrmCntFrm$Eig*100), "%", sep="")
TrmCntFrm$Trn = paste(round(TrmCntFrm$Trn*100), "%", sep="")
TrmCntFrm$Int = paste(round(TrmCntFrm$Int*100), "%", sep="")

TrmCntFrm = TrmCntFrm %>% select(Ntw,Com,Deg,Btw,Cls,Eig,Trn,Int)
TrmCntFrm = TrmCntFrm %>% filter(Com!="Age")
TrmCntFrm$Com[TrmCntFrm$Com=="Hum"] = "GenAge.Hum"
TrmCntFrm$Com[TrmCntFrm$Com=="Mod"] = "GenAge.Mod"
TrmCntFrm$Com[TrmCntFrm$Com=="Arc"] = "Diseases"
TrmCntFrm$Com[TrmCntFrm$Com=="Imm"] = "Immune Disorders"
TrmCntFrm$Com[TrmCntFrm$Com=="Oth"] = "Others"
TrmCntFrm$Eig = NULL
colnames(TrmCntFrm) = c("Network", "Genes", "Degree Centrality", "Betweenness Centrality", "Closeness Centrality",
                        "Clustering Coefficient", "Percentage of Neighbour Genes related with Diseases") 

NtwArr = TrmCntFrm$Network %>% unique()
NtwNum = c(3,4,2,1) 
NtwFrm = data.frame(Network=NtwArr,NtwNum)

TrmCntFrm = TrmCntFrm %>% 
  merge(NtwFrm) %>%
  group_by(Network) %>%
  arrange(NtwNum, Genes) %>%
  ungroup() 
TrmCntFrm$NtwNum = NULL

write.csv(TrmCntFrm,"C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/Centrality_Table.csv")


################################################################################
### FUNCTIONS
################################################################################

calculate_overlaps <- function(graph, gene_set) {
  # Initialize a data frame to store the results
  result_df <- data.frame(Gene = V(graph)$name, Degree = igraph::degree(graph,mode="out"), Number_of_overlaps = NA)
  #result_df <- data.frame(Gene = V(graph)$name, Number_of_overlaps = NA)
  
  # Loop through each gene (node) in the graph
  NumNod = vcount(graph)
  i=3
  for (i in 1:NumNod) {
    
    paste(i,NumNod) %>% print()
    
    # Get the first-order neighbors of the gene
    neighbors <- neighbors(graph, i, mode="out")
    
    # Calculate the overlap with the gene set (excluding the gene itself)
    #overlap <- length(intersect(neighbors$name, gene_set)) - (V(graph)$name[i] %in% gene_set)
    overlap <- length(intersect(neighbors$name, gene_set)) 
    
    # Store the result in the data frame
    result_df$Number_of_overlaps[i] <- overlap
  }
  
  result_df$Rat = result_df$Number_of_overlaps/result_df$Degree
  
  return(result_df)
}


full_transitivity <- function(graph) {
  # Check if the input is a graph object
  if (!inherits(graph, "igraph")) {
    stop("The input is not an igraph graph object.")
  }
  
  # Calculate local transitivity for non-isolated vertices
  transitivity_scores <- transitivity(graph, type = "local")
  
  # Get all vertex names
  all_vertex_names <- V(graph)$name
  
  # Initialize a named vector with NA values for all vertices
  full_scores <- setNames(rep(NA, length(all_vertex_names)), all_vertex_names)
  
  # Identify non-isolated vertices (vertices with at least one edge)
  non_isolated <- which(igraph::degree(graph) > 0)
  
  # Update the full_scores with the transitivity scores for non-isolated vertices
  full_scores[names(full_scores) %in% all_vertex_names[non_isolated]] <- transitivity_scores
  
  return(full_scores)
}
  

transform_data <- function(df) {
  # Reshape the data and compute the mean score for each combination of Ntw, Com, and Cnt
  result_df <- df %>%
    group_by(Ntw, Com, Cnt) %>%
    summarize(Mean_Scr = mean(Scr, na.rm = TRUE), .groups = 'drop') %>%
    spread(key = Cnt, value = Mean_Scr)
  
  return(result_df)
}
  