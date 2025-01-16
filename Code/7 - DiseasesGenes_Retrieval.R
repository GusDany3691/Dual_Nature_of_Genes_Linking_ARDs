library(tidyverse)
library(data.table)
library(VariantAnnotation)
library(GenomicRanges)
library(tidyverse)
library(qqman)
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

Tss = FALSE

##### GENES FRAME ####################################################################################################################

genemap = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/genemap.rds") #OK
genes.granges = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/GR_nTss_H.rds") #PSEUDO OK
GenesFrame = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/HumanGenesFrameComplete.rds") #OK
DiseasesCategories = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Showcase/DiseasesCategories.csv") #UK BIOBANK
colnames(DiseasesCategories)

colnames(GenesFrame)[colnames(GenesFrame) == "hgnc_symbol"] ="Genes"

##############################################################################################################################
### OTHER ####################################################################################################################
##############################################################################################################################

Nodes = DiseasesCategories[DiseasesCategories$Ageing %in% c(1,2,3,4),]$Node
Names = DiseasesCategories [DiseasesCategories$Ageing %in% c(1,2,3,4),]$Meaning
Nodes = as.character(Nodes)

Interval_Gene = list()
sgpList_d10K = list()
sgpList_g10K = list()

N_Nodes = length(Nodes)

for(i in 1:N_Nodes){
  
  LupTimeStart <- Sys.time()
  
  paste("### PHENOTYPE: ", i, "/", N_Nodes," #############################################", sep="") %>% print()
  noquote("")
  
  Node = Nodes[i]
  Name = Names[i]
  
  # GET PHENOTYPE DATA
  
  FilePath = paste("G:/UK_Biobank2/Data/Retrieved/Ukbb/bolt/a", Node, ".imp.stats", sep = "")
  PathExists = FilePath %>% file.exists()
  
  if(PathExists){
    
    # GET PHENOTYPE DATA
    gwasRes <- read_tsv( FilePath )
    Name_short = str_remove(Name, "Pilling_2017_UKB_GWAS_")
    PgwasRes = gwasRes
    gwasRes$P_BOLT_LMM = as.numeric(gwasRes$P_BOLT_LMM)
    gwasRes = AdaptGwasFrame(gwasRes)
    sgwasRes = gwasRes[gwasRes$pvalue < 5E-8,]
    
    sgpList_g10K[[as.character(Node)]] = summary2genes_intervals(sgwasRes, Node, Name_short, thr = 10000, GenesFrame, genes.granges)
    
  }
}

GenPhnFrm_g10K = FrameFromList(sgpList_g5K) %>% dplyr::select(PhenoCode, Gene) %>% unique()
colnames(GenPhnFrm_g10K) = c("Phn","Gen")
GenPhnFrm_g10K$Phn = paste("D",GenPhnFrm_g10K$Phn,sep="")

MhcGenArr = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Ranges/MhcGenArr.rds")
GenPhnFrm_g10K_NoMhc = GenPhnFrm_g10K %>% filter(!(Gen %in% MhcGenArr))

# INTEGRATING GENAGE ###########################################################

GenAgeHum = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeHum_Genes.rds")
GenAgeMod = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/HAGR/GenAgeMod_Genes.rds")

GenAgeHumFrm = data.frame(Phn='HumAge',Gen=GenAgeHum)
GenAgeModFrm = data.frame(Phn='ModAge',Gen=GenAgeMod)

PhnAgeFrm = GenPhnFrm_g10K_NoMhc %>% rbind(GenAgeHumFrm) %>% rbind(GenAgeModFrm)

saveRDS(PhnAgeFrm,"C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnAgeFrm.rds")

HrcFrm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Showcase/PhnHrc_FrmEnd.csv")
HrcFrm$X = NULL
PhnComFrm = HrcFrm %>% dplyr::select(Nod,Mng,RutMng)
colnames(PhnComFrm) = c('Phn','Mng','Com')

GenPhnComFrm = merge(PhnAgeFrm,PhnComFrm,by='Phn') %>% 
  dplyr::select('Phn','Mng','Com', 'Gen') %>% 
  filter(!(Com %in% c("human ageing", "model ageing homologs")))

saveRDS(GenPhnComFrm,"C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnComFrm.rds")



###################################################################################################################################
##### FUNCTIONS ###################################################################################################################
###################################################################################################################################

gwas2GRanges <- function(gwasRes, SNP = "RefSNP_id", start = "BP", chr = "CHR", cols2retain = c('ALLELE0', 'ALLELE1'), genome='hg19' ){
  library(tidyverse)
  library(GenomicRanges)
  gwasRes <- gwasRes %>%
    dplyr::rename(RefSNP_id = SNP,
                  start = start,
                  chr = chr) %>%
    dplyr::mutate(end=`start`)
  gwasRes <- gwasRes %>%
    dplyr::select( RefSNP_id, chr, start, end, cols2retain)
  snpinf <- makeGRangesFromDataFrame(gwasRes)
  genome(snpinf)=genome
  for( colnm in c('RefSNP_id', cols2retain)){
    values(snpinf)[[colnm]]=gwasRes[[colnm]]
  }
  return(snpinf)
}



gwas2GRanges_chr <- function(gwasRes, SNP = "RefSNP_id", start = "BP", chr = "CHR", cols2retain = c('ALLELE0', 'ALLELE1'), genome='hg19' ){
  library(tidyverse)
  library(GenomicRanges)
  gwasRes <- gwasRes %>%
    dplyr::rename(RefSNP_id = SNP,
                  start = start,
                  chr = chr) %>%
    dplyr::mutate(end=`start`)
  gwasRes <- gwasRes %>%
    dplyr::select( RefSNP_id, chr, start, end, cols2retain)
  snpinf <- makeGRangesFromDataFrame(gwasRes)
  if (all(substr(seqlevels(snpinf),1,3)!='chr')) {
    seqlevels(snpinf)=paste('chr',seqlevels(snpinf),sep='')
  }
  genome(snpinf)=genome
  for( colnm in c('RefSNP_id', cols2retain)){
    values(snpinf)[[colnm]]=gwasRes[[colnm]]
  }
  return(snpinf)
}


snp2gene_eQTL <- function(gwasRes, eQTLfile, genemap){
  tissue <- strsplit(sapply(strsplit(eQTLfile,'/'),function(x)x[length(x)]),'[.]')[[1]][1]
  eqtl <- read_tsv(eQTLfile) %>%
    inner_join(gwasRes) %>%
    dplyr::select(SNP,CHR, BP, ALLELE1, ALLELE0, gene_id, slope, pval_beta, tissue)%>%
    unique()
  genemap <- genemap %>%
    unique()%>%
    dplyr::mutate(entrezgene=as.character(entrezgene))%>%
    dplyr::rename(gene_id = ensembl_gene_id)
  eqtl <- left_join(eqtl, genemap)
  eqtl$entrezgene[eqtl$entrezgene=='']=NA
  eqtl$hgnc_symbol[eqtl$hgnc_symbol=='']=NA
  eqtl$description[eqtl$description=='']=NA
  eqtl
}



FrameFromList <- function(sgpList){
  sgpFrame = data.frame()
  for (i in 1:length(sgpList)){
    Frame = sgpList[[i]]
    if(dim(Frame)[1] > 0)
      sgpFrame = rbind(sgpFrame, Frame)
  }
  return(sgpFrame)
}



range2GRanges <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = as.numeric(df[,2]), end = as.numeric(df[,3]))
  )
  return(gr)
}


# EXTEND RANGES
extend <- function(x, upstream=0, downstream=0)     
{
  names = names(x)
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
  names(x) = names
  return(x)
}


# ADAPT HEALTHSPAN
AdaptHealthspan <- function(gwasRes)     
{
  names(gwasRes)[names(gwasRes) == 'SNPID'] <- 'SNP'
  names(gwasRes)[names(gwasRes) == 'chr'] <- 'CHR'
  names(gwasRes)[names(gwasRes) == 'pos'] <- 'BP'
  names(gwasRes)[names(gwasRes) == 'EA'] <- 'ALLELE1'
  names(gwasRes)[names(gwasRes) == 'RA'] <- 'ALLELE0'
  names(gwasRes)[names(gwasRes) == 'EAF'] <- 'A1FREQ'
  names(gwasRes)[names(gwasRes) == 'beta'] <- 'BETA'
  names(gwasRes)[names(gwasRes) == 'se'] <- 'SE'
  return(gwasRes)
}



AdaptGwasFrame <- function(gwasRes)     
{
  names = colnames(gwasRes)
  if ("X.log10.p.value." %in% names){
    gwasRes$pvalue = 10^(gwasRes$X.log10.p.value.)
  }
  if ("P_BOLT_LMM" %in% names){
    gwasRes$pvalue = gwasRes$P_BOLT_LMM
    gwasRes$X.log10.p.value. = -log10(gwasRes$pvalue)
  }
  return(gwasRes)
}



GetSnp2GeneDistance <- function(sgpFrame)     
{
  distance = c()
  for(i in 1:dim(sgpFrame[1])){
    if(sgpFrame[i,"BP"] < sgpFrame[i,"start_position"]){
      distance[i] = paste("-",as.character(sgpFrame[i,"start_position"] - sgpFrame[i,"BP"]))
    }
    else if(sgpFrame[i,"BP"] > sgpFrame[i,"end_position"]){
      distance[i] = paste("+",as.character(sgpFrame[i,"BP"] - sgpFrame[i,"end_position"]))
    }
    else{
      distance[i] = "within"
    }
  }
  sgpFrame$Snp2Gene_distance = distance
  return(sgpFrame)
}



GetSnp2GeneDistance_Apply <- function(SnpGeneFrame, Tss){
  
  
  SnpGeneFrame = SnpGeneFrame %>% dplyr::mutate(
      Snp2Gene_distance = case_when(
        BP < start_position ~ BP - start_position,
        BP > end_position ~ BP - end_position, 
        TRUE ~ 0
        )
      )
  
  if(Tss){
    SnpGeneFrame = SnpGeneFrame %>% dplyr::mutate(Snp2Tss_distance = BP - transcription_start_site)
  }
  
  
  return(SnpGeneFrame)
}  
 


# ADAPT PVALUE APPLY
GetSnp2GeneDistance_Apply_For <- function(SnpGeneFrame, Tss){
  
  Snp2Gene_distance = c()
  
  if(!Tss){
    for(k in 1:dim(SnpGeneFrame)[1]){
      if(SnpGeneFrame$BP[k] < SnpGeneFrame$start_position[k]){
        Snp2Gene_distance[k] = SnpGeneFrame$BP[k] - SnpGeneFrame$start_position[k]
      }else if(SnpGeneFrame$BP[k] > SnpGeneFrame$end_position[k]){
        Snp2Gene_distance[k] = SnpGeneFrame$BP[k] - SnpGeneFrame$end_position[k]
      } else{
        Snp2Gene_distance[k] = 0
      }
    }
  }
  
  SnpGeneFrame$Snp2Gene_distance = Snp2Gene_distance
  
  if(Tss){
    SnpGeneFrame = SnpGeneFrame %>% dplyr::mutate(Snp2Tss_distance = BP - transcription_start_site)
  }
  
  return(SnpGeneFrame)
}  



summary2genes_intervals <- function(sgwasRes, Node, Name, thr, GenesFrame, genes.granges){
  significants = dim(sgwasRes)[1]
  if (significants > 0){
    sgwasRes = AdaptGwasFrame(sgwasRes)
    
    gwas_as_GR <- gwas2GRanges(sgwasRes, SNP = 'SNP',start = 'BP',chr = 'CHR',genome = 'hg19')
    
    r1 = extend(genes.granges, upstream=thr, downstream=thr) #genes.granges
    r2 = gwas_as_GR
    
    names(r2) = r2$RefSNP_id
    
    overlap <- GenomicRanges::findOverlaps(r1, r2)
    snps <- names(r2)[overlap@to]
    genes <- names(r1)[overlap@from]
    snps_genes_frame = data.frame(snps, genes)
    
    colnames(snps_genes_frame) = c("SNP", "Genes")

    SnpGeneFrame = merge(sgwasRes[,c("CHR", "BP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "X.log10.p.value.", "pvalue", "SNP")], snps_genes_frame,by="SNP",all=TRUE)

    SnpGeneFrame = merge(SnpGeneFrame,GenesFrame,by="Genes",all=FALSE)
    
    if(dim(SnpGeneFrame)[1] > 0){

      SnpGeneFrame$PhenoName = Name
      SnpGeneFrame$PhenoCode = Node
      
      
      ChrCond = SnpGeneFrame$chromosome_name == SnpGeneFrame$CHR
      
      SnpGeneFrame$ChrCond = ChrCond
      
      SnpGeneFrame$Interval_any = (SnpGeneFrame$BP >= SnpGeneFrame$start_position - thr) & (SnpGeneFrame$BP <= SnpGeneFrame$end_position + thr)
      SnpGeneFrame$Interval_gene = (SnpGeneFrame$BP >= SnpGeneFrame$start_position) & (SnpGeneFrame$BP <= SnpGeneFrame$end_position)
      SnpGeneFrame$Interval_near = xor(SnpGeneFrame$Interval_gene, SnpGeneFrame$Interval_any)
      
      SnpGeneFrame = GetSnp2GeneDistance(SnpGeneFrame) 
      
      # -------------------------------------------------------------------------------------------------
      
      gwas_as_GR <- gwas2GRanges_chr(sgwasRes, SNP = 'SNP',start = 'BP',chr = 'CHR',genome = 'hg19')
      # FURTHER ASSOCIATIONS
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
      varInfo <- unique(granges(gwas_as_GR))
      allvar <- locateVariants(varInfo, txdb,
                               AllVariants(promoter=PromoterVariants(upstream=thr,
                                                                     downstream=thr)))
      overs <- findOverlaps(gwas_as_GR,allvar)
      
      NumData = dim(as.data.frame(allvar, row.names = NULL))[1]
      
      if (NumData > 0 ) {
        
        for(colnm in colnames(mcols(gwas_as_GR))){
          mydat <-(mcols(gwas_as_GR)[[colnm]])[queryHits(overs)]
          mcols(allvar)[[colnm]] = '*'
          mcols(allvar)[[colnm]][subjectHits(overs)] <- mydat
        }
        geneids <- unique(allvar$GENEID)
        geneids <- geneids[complete.cases(geneids)]
        genemap <- genemap %>% unique() %>%
          mutate(entrezgene=as.character(entrezgene))
        
        mcols(allvar)$RefSNP_id
        mcols(allvar)$LOCATION
        
        SnpDef = unique(data.frame("SNP" = mcols(allvar)$RefSNP_id, "Location" = mcols(allvar)$LOCATION))
        
        GeneDef = unique(data.frame("Genes" = genemap$hgnc_symbol, "Description" = genemap$description))

        SnpGeneFrame = merge(SnpGeneFrame, unique(SnpDef), by = "SNP", all.x = TRUE, all.y = FALSE)
        SnpGeneFrame = merge(SnpGeneFrame, GeneDef, by = "Genes", all.x = TRUE, all.y = FALSE)
        
        
        SnpGeneFrame = SnpGeneFrame[,c("PhenoCode", "PhenoName", "SNP", "pvalue", "X.log10.p.value.", "CHR", "BP", "Location",
                                      "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE",
                                      "Genes", "chromosome_name", "start_position", "end_position", "ChrCond", "Interval_any",
                                      "Interval_gene", "Interval_near", "Snp2Gene_distance", "gene_biotype", "Description")] #Z would be penultime if required
        
        SnpGeneFrame <- SnpGeneFrame %>% 
          dplyr::rename(CHR_SNP = CHR,
                        BP_SNP = BP,
                        Location_SNP = Location,
                        Gene = Genes,
                        CHR_Gene = chromosome_name,
                        Gene_Start = start_position,
                        Gene_End = end_position,
                        CHR_Congruency = ChrCond,
                        Gene_Description = Description)

      } else {
        SnpGeneFrame = data.frame()
      }
    } else {
      SnpGeneFrame = data.frame()
    }
  } else {
    SnpGeneFrame = data.frame()
  }
  
  return(SnpGeneFrame)
}



MinP_Gene = function(G, Expanded_Summary){
  min_p = Expanded_Summary %>% filter(Gen == G) %>% pull(pvalue) %>% min()
  Expanded_Summary %>% filter(Gen == G) %>% mutate(min_pvalue = min_p)
  #return(Expanded_Summary)
}
