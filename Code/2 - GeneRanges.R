library(biomaRt)
Tss = FALSE
Source = "H" # H - HGNC, E - Ensemble

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

colnames(genes) = c("hgnc_symbol", "Ensemble", "Entrez", "chr", "transcription_start_site", "start_transcript", "end_transcript", "strand", "Function") 

if(!Tss){
  genes$transcription_start_site = NULL
  genes = unique(genes)
}

if(Source == "H"){
  genes = genes %>% filter(hgnc_symbol != "")
}

#genes$start = genes$tss
#genes$end = genes$tss

genes$start = genes$start_transcript = genes$start_transcript
genes$end = genes$end_transcript = genes$end_transcript

genes$strand = str_replace(genes$strand, "-1", "-")
genes$strand = str_replace(genes$strand, "1", "+")

genes = genes[order(genes$chr),]

if(Tss){
  GeneCols = genes[,c("hgnc_symbol", "Ensemble", "Entrez", "transcription_start_site", "start_transcript", "end_transcript")]
} else{
  GeneCols = genes[,c("hgnc_symbol", "Ensemble", "Entrez", "start_transcript", "end_transcript")]
}


##### SNPS TO GENE RANGES #################################################################

# GENE
geneinf = makeGRangesFromDataFrame(genes)
genome(geneinf)='hg19'
for( colnm in colnames(GeneCols)){
  values(geneinf)[[colnm]]=GeneCols[[colnm]]
}

if(Source == "H"){
  names(geneinf) = geneinf$hgnc_symbol
}
if(Source == "E"){
  names(geneinf) = geneinf$Ensemble
}


if(Tss){
  TssText = "Tss"
}else{
  TssText = "nTss"
}

#GRtext = paste("E:/UK_Biobank/Data/Retrieved/Biomart/", "GR_", TssText, "_", Source,".rds", sep="")
GRtext = paste("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/", "GR_", TssText, "_", Source,".rds", sep="")
saveRDS(geneinf, GRtext)


