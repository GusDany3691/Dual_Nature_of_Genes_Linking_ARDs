library(tidyverse)
library(data.table)
library(stringr)

Cor - Correlation (e.g., Cor_long, Cor)
Cox - Coexpression (e.g., CoxPrwFrm90, CoxPrwFrm95)
Ens - Ensemble (e.g., EnsHgnFrm, DupEnsNam90)
Nam - Name (e.g., DupEnsNam90, DupEnsNam95)
Dup - Duplicated (e.g., DupEnsNam90, DupEnsNam95)
Frm - Frame (e.g., EnsHgnFrm, FltFltEnsHgnFrm, CoxPrwFrm90)
Hgn - Hgnc (e.g., EnsHgnFrm, FltEnsHgnFrm)
Flt - Flat (i.e., matrix turned into an array, flattered) (e.g., FltEnsHgnFrm, FltFltEnsHgnFrm)
Gen - Gene (e.g., AllGen90, AllGen95)
All - All (e.g., AllGen90, AllGen95)
Abs - Absolute value (e.g., abs_Cor_long)
Pre - Previous (e.g., PreEnsGenCor, inferred from the commented code)
Prw - Pairwise (e.g., CoxPrwFrm90, CoxPrwFrm95)
Phn - Phenotype



################################################################################

Cor = read_tsv("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Networks/human_genes_correlation_matrix.tsv")

################################################################################
# COX NETWORKS 
################################################################################

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# LONG FORMAT OF FRAME
Cor_long <- pivot_longer(Cor, cols = -Ensemble, names_to = "column", values_to = "value")
abs_Cor_long = Cor_long %>% mutate(value = abs(value))

# Coexpression 90
filtered_Cor_long90 = abs_Cor_long %>% filter(value > 0.90) %>% Dupyr::select(Ensemble, column)
colnames(filtered_Cor_long90) = c("row","col")
write.table(filtered_Cor_long90, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/CoxPrwFrm90.csv")

AllGen90 = c(filtered_Cor_long90$row, filtered_Cor_long90$col)

DupEnsNam90 = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                  filters = 'ensembl_gene_id', 
                  values = AllGen90, 
                  mart = ensembl)

DupEnsNam90 %>% saveRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/DupEnsNam90.rds")



# Coexpression 95
filtered_Cor_long95 = abs_Cor_long %>% filter(value > 0.95) %>% Dupyr::select(Ensemble, column)
colnames(filtered_Cor_long95) = c("row","col")
write.table(filtered_Cor_long95, "C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/CoxPrwFrm95.csv")

AllGen95 = c(filtered_Cor_long95$row, filtered_Cor_long95$col)

DupEnsNam95 = getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), 
                    filters = 'ensembl_gene_id', 
                    values = AllGen95, 
                    mart = ensembl)

DupEnsNam95 %>% saveRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/DupEnsNam95.rds")


################################################################################
# COEXPRESSION MATRIX FOR AGEING AND DISEASE-RELATED GENES 
################################################################################

# Extract Ensembl gene IDs from the row and column names of PreEnsGenCor matrix
#gene_ids <- unique(c(rownames(PreEnsGenCor), colnames(PreEnsGenCor)))
gene_ids <- unique(c(rownames(Cor), colnames(Cor)))

EnsHgnFrm = readRDS("/mnt/work/larder/inesl/Gustavo_Work/GeneCorr/Ret/EnsHgnFrm.rds")

row.names(EnsHgnFrm) = EnsHgnFrm$ensembl_gene_id
FltEnsHgnFrm = EnsHgnFrm %>% filter(external_gene_name != "")

NonDupEns = FltEnsHgnFrm$ensembl_gene_id[!duplicated(FltEnsHgnFrm$external_gene_name)]
FltFltEnsHgnFrm = FltEnsHgnFrm[NonDupEns,]

FltFltEnsHgnFrm$external_gene_name %>% table() %>% sort() %>% rev() %>% head(100)


# MATRIX WITH HGNA NAMES -------------------------------------------------------
HgnGenCor = Cor[FltFltEnsHgnFrm$ensembl_gene_id,FltFltEnsHgnFrm$ensembl_gene_id]

EnsRow = row.names(HgnGenCor)
EnsCol = colnames(HgnGenCor)

HgnRow = EnsHgnFrm[EnsRow,"external_gene_name"]
HgnCol = EnsHgnFrm[EnsCol,"external_gene_name"]

# Replace Ensembl gene IDs with gene names in the row and column names of the PreEnsGenCor matrix
rownames(HgnGenCor) <- HgnRow
colnames(HgnGenCor) <- HgnCol

# FILTER FOR DISEASE AND AGEING-RELATED GENES ----------------------------------

GenPhnFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenPhnFrm_g10K_NoMhc.rds")

GenArr = GenPhnFrm$Gen %>% unique()
PhnGenArr = intersect(GenArr, row.names(HgnGenCor))

PhnHgnGenCor = HgnGenCor[PhnGenArr,PhnGenArr]

#write.table(MpaHgnGenCor,"/mnt/work/larder/inesl/Gustavo_Work/GeneCorr/Ret/MpaHgnGenCor.txt")
saveRDS(PhnHgnGenCor,"C:/Users/Usuario/Desktop/Nature/Data/Generated/Networks/Coexpression/PhnHgnGenCor.rds")