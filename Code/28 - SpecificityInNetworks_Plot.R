# Dir - Direct
# Ind - Indirect

library(dplyr)
library(cowplot)
library(ggplot2)
library(rstatix)
library(ggpubr)

TauFrm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Specificity/Tau_gene_V8.csv") 
EnsHgnFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Retrieved/Ranges/EnsHgnFrm.rds")
row.names(EnsHgnFrm) = EnsHgnFrm$ensembl_gene_id
TauFrm$Gen = EnsHgnFrm[TauFrm$gene_id,"external_gene_name"]
TauFrm = TauFrm %>% filter(!is.na(tau)) %>% select(Gen, tau) %>% rename(Tau=tau)

DirIndIntFrm_Pin = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Generated/EndPin/DirIndIntFrm.rds') 

# Pleiotropy
DirFrm = DirIndIntFrm_Pin %>% filter(DirIntCom>0) %>% select(Gen,DirIntCom) %>% rename(Num=DirIntCom) %>%
  mutate(TypInt = "Dir") %>% 
  mutate(TypNum = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "EndAll") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, TypNum, Dts)

# PPI iARC-Interactions
IndFrm_Pin = DirIndIntFrm_Pin %>% 
  filter(IndIntCom>0) %>% select(Gen,IndIntCom) %>% rename(Num=IndIntCom) %>%
  mutate(TypInt = "Ind") %>% 
  mutate(TypNum = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "EndPin") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, TypNum, Dts)

# COX90 iARC-Interactions
IndFrm_C90 = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Generated/Endc90/DirIndIntFrm.rds') %>%
  filter(IndIntCom>0) %>% select(Gen,IndIntCom) %>% rename(Num=IndIntCom) %>%
  mutate(TypInt = "Ind") %>% 
  mutate(TypNum = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "EndC90") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, TypNum, Dts)

# COX95 iARC-Interactions
IndFrm_C95 = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Generated/EndC95/DirIndIntFrm.rds')  %>%
  filter(IndIntCom>0) %>% select(Gen,IndIntCom) %>% rename(Num=IndIntCom) %>%
  mutate(TypInt = "Ind") %>% 
  mutate(TypNum = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "EndC95") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, TypNum, Dts)

# KEGG iARC-Interactions
IndFrm_Kgo = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Generated/EndKgo/DirIndIntFrm.rds') %>%
  filter(IndIntCom>0) %>% select(Gen,IndIntCom) %>% rename(Num=IndIntCom) %>%
  mutate(TypInt = "Ind") %>% 
  mutate(TypNum = ifelse(Num>=4,"4+",Num)) %>% 
  mutate(Dts = "EndKgo") %>%
  merge(TauFrm) %>%
  select(Gen, Tau, TypInt, Num, TypNum, Dts)

TauDirIndFrm = rbind(DirFrm,IndFrm_Pin,IndFrm_C90,IndFrm_C95,IndFrm_Kgo)

TauDirIndFrm %>% filter(Dts=='EndAll') %>% dim()

DtsArr = TauDirIndFrm$Dts %>% unique()
DtsLen = length(DtsArr)

# Plot parameters
Alp = 0.2
TckAng = 45
TstTxtSiz = 3
AnnTxtSiz = 3
YtstJmp = 0.1
TstSep = 0.07
Ylab = "Tau"
MinPnt = 1.05
DotArr = c(1,3.5,1,3,3.5,3.5)
LenAnn = 4

# Plot list and frame
TauDirIndLst = list()
MenFrm = data.frame()
z=1
for(z in 1:DtsLen){
  
  print(z)
  
  DtsQry = DtsArr[z] 
  
  TypIntQry = TauDirIndFrm %>%
    filter(Dts == DtsQry) %>%
    pull(TypInt) %>%
    unique()
  
  if(DtsQry == "EndAll"){
    DtsTxt = "Pleiotropy"
  }
  if(DtsQry == "EndPin"){
    DtsTxt = "PPI"
  }
  if(DtsQry == "EndKgo"){
    DtsTxt = "KEGG"
  }
  if(TypInt == "EndC90"){
    DtsTxt = "COX.90"
  }
  if(DtsQry == "EndC95"){
    DtsTxt = "COX.95"
  }
  
  TstFrm = TauDirIndFrm %>% 
    filter(Dts == DtsQry) %>%
    select(TypNum,Tau,Dts)
  
  pwc = TstFrm %>% pairwise_t_test( Tau ~ TypNum , p.adjust.method = "bonferroni")

  if(TypIntQry == "Dir"){
    TstFrm$TypInt = paste("Tissue specificity\nARC-Pleiotropy", sep="")
    Xlab = "Number of directly\nconnected ARCs"
  } else{
    TstFrm$TypInt = paste("Tissue specificity\niARC-Interactor (",DtsTxt,")", sep="")
    Xlab = "Number of indirectly\nconnected ARCs"
  }
  
  # Statistics
  grouped_means <-  TstFrm %>%
    group_by(TypNum) %>%
    summarise(MenTau = mean(Tau, na.rm = TRUE))
  grouped_means
  
  
  grouped_medians <- TstFrm %>%
    group_by(TypNum) %>%
    summarise(MedTau = median(Tau, na.rm = TRUE))
  
  # Mean Tau Frame
  sMenFrm = grouped_means %>% merge(grouped_medians, by='TypNum')
  sMenFrm$Network = DtsQry
  MenFrm = rbind(MenFrm,sMenFrm)
  
  # Convert the `TypNum` column to a factor with the desired order
  TstFrm$TypNum <- factor(TstFrm$TypNum, levels = c("1", "2", "3", "4+"))
  
  # BOX PLOT 
    p = ggboxplot(TstFrm, x = "TypNum", y = "Tau",
                  fill="#F8776D",
                  size = 0.9,
                  facet.by = "TypInt",
                  combine = TRUE,
                  outlier.shape = NA,
    ) 
  
  
  p = p+geom_dotplot(binaxis='y', stackdir='center',
                     position=position_dodge(1), dotsize=0.004*(DotArr[z]), color="#202020")

  p = p + theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1)) 
  
  p = ggpar(p, xlab =Xlab, ylab = Ylab)+
    font("xlab", size = 8, color = "black")+
    font("ylab", size = 8, color = "black")+
    font("xy.text", size = 8, color = "black") +
    theme(legend.position = "none")
  
  p = p + grids(axis = c("xy"), linetype = "dashed",  color = "grey", size = NULL) #x, y, xy
  
  pwc$y.position = c(MinPnt, MinPnt+(1*YtstJmp), MinPnt, MinPnt+(2*YtstJmp), MinPnt+(3*YtstJmp),MinPnt)
  pwc[[".y."]] = pwc$y.position
  
  # X Separison - Tests
  pwc = pwc %>% mutate(xmin = c(1,1,2,1,2,3)+TstSep )
  pwc = pwc %>% mutate(xmax = c(2,3,3,4,4,4)-TstSep )
  
  p = p + stat_pvalue_manual(pwc, label = "p.adj.signif", tip.length = 0, step.increase = 0.00, label.size = TstTxtSiz) 
  
  p = p +
    annotate("text",
             x = 1:LenAnn,
             y = -0.1,
             label = TstFrm$TypNum %>% table() %>% format_with_commas(),
             col = "black",
             size = AnnTxtSiz,
             vjust = - 1)
  
  p = p + coord_cartesian(ylim = c(-0,1.36))
  
  if(!(z %in% c(2,5))){
    p = p + ylab(NULL)
  }
  
  p = p + scale_x_discrete(labels = c("1" = "1", "2" = "2","3" = "3", "4+" = "4+") )
  
  TauDirIndLst[[DtsQry]] = p

}


plot_grid(TauDirIndLst$EndPin, TauDirIndLst$EndC95, TauDirIndLst$EndC90, 
          TauDirIndLst$EndKgo, TauDirIndLst$EndAll, NULL,
          labels = c("a","b","c","d","e","","","",""), label_size = 12, ncol=3) 


write.csv(TauDirIndFrm,"C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/TauDirIndFrm.csv")
write.csv(MenFrm,"C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/TauDirIndMeanFrm.csv")


################################################################################
# FUNCTIONS
################################################################################

format_with_commas <- function(x) {
  # Check if x is a numeric vector or a matrix
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector or matrix.")
  }
  
  # Format each element in x with commas
  formatted <- format(x, big.mark = ",", scientific = FALSE, justify = "none")
  
  # Remove leading spaces by trimming the formatted strings
  formatted <- trimws(formatted)
  
  # Return the formatted vector
  return(formatted)
}
