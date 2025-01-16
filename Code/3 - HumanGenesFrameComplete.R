library(biomaRt)
library(dplyr)
library("GenomicRanges")

### HUMAN GENES FRAME COMPLETE #################################################
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(
  attributes=c("hgnc_symbol", "ensembl_gene_id" ,"entrezgene_id","chromosome_name", 
               "transcription_start_site", "start_position","end_position", "strand", "gene_biotype"),
  mart = mart)
#attributes=c("hgnc_symbol","chromosome_name", "transcription_start_site", "start_position","end_position", "strand"),
#mart = mart)

genes = genes[genes$chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                           "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                           "20", "21", "22"),]
colnames(genes) = c("Gene", "Ensemble", "Entrez", "chr", "tss", "start_transcript", "end_transcript", "strand", "Function") 

#saveRDS(genes, "E:/UK_Biobank/Data/Retrieved/Biomart/HumanGenesFrameComplete.rds")
saveRDS(genes, "C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/HumanGenesFrameComplete.rds")


### Query for genemap ##########################################################
genemap <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "description"),
  mart = ensembl
)

#saveRDS(genemap, "E:/UK_Biobank/Data/Retrieved/Biomart/genemap.rds")
saveRDS(genemap, "C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/genemap.rds")

### EXCLUDING MHC ##############################################################


#GenRanFrm = readRDS('G:/UK_Biobank2/Data/Generated/Sorted/Ranges/GenesRanges.rds')
GenRanFrm = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/GenesRanges.rds')

#GenRanFrm = readRDS('G:/UK_Biobank2/Data/Generated/Sorted/Ranges/HumanGenesFrame.rds')
GenRanFrm = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/HumanGenesFrame.rds')

HlaGenFrm = GenRanFrm %>% filter(chromosome_name == 6 & start_position >= 28477797 & end_position <= 33448354)

MhcGenArr = HlaGenFrm$hgnc_symbol %>% unique()

#saveRDS(MhcGenArr,"G:/UK_Biobank2/Data/Generated/Sorted/Ranges/MhcGenArr.rds")
saveRDS(MhcGenArr,"C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/MhcGenArr.rds")


### HGNC AND ENSEMBLE RELATIONSHIP #############################################

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# 1. Convert from ensembl.gene to gene.symbol
#genes <- c("ENSG00000150676", "ENSG00000099308", "ENSG00000142676", "ENSG00000180776", "ENSG00000108848", "ENSG00000277370", "ENSG00000103811", "ENSG00000101473")

G_list <- getBM(filters = "ensembl_gene_id", 
                attributes = c("ensembl_gene_id", "external_gene_name"),
                mart = mart)

# saveRDS(G_list, "C:/Users/gusdany/Desktop/Transactions/EnsHgnFrm.rds")
saveRDS(G_list, "C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/EnsHgnFrm.rds")


