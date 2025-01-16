library(rstatix)
library(ggpubr)
library(dplyr)
library(rje)
library(cowplot)

GenTypCmpFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndAll/GenTypCmpFrm.rds")

TypArr = GenTypCmpFrm$TypInt %>% unique()

### Pairwise comparisons ###########################################################################

# QUERY VARIABLES
Xvar = "GenTyp" #AgeDst, GenTyp
Yvar = "GenNg1"#"GenNg1" #"GenNg1", #"GenAvg"


# - GROUPS SORTED ------------------------------------------------------------------------------------
FltElm = c("DisGen", "AgeGen", "ModGen", "DisNgb", "aOthGen") 

xTckLab = c("Diseases", "GenAge.Hum", "GenAge.Mod", "Neighbours", "Others")
FltElmLen = length(FltElm)
FltElmGr1 = data.frame(group1=FltElm, xmin=1:FltElmLen)
row.names(FltElmGr1) = FltElmGr1$group1
FltElmGr2 = data.frame(group2=FltElm, xmax=1:FltElmLen)
row.names(FltElmGr2) = FltElmGr2$group2

### INITIAL VARIABLES -- DisGen, AgeGen, aExcNgb -- GenNgb #####################

WixTst = TRUE
TstTxt = "All" 
AnnNumMen = TRUE #Number of genes

Ymin = 0 

#yLimMax = c(9.7, 11.5, 11.5, 11.5, 11.5, 11.5) 
#YtstMin = c(7.6, 9.0, 9.0, 9.0, 9.0, 9.0)

yLimMax = c(8.7, 11.5, 11.5, 11.5, 11.5, 11.5) 
YtstMin = c(6.6, 9.0, 9.0, 9.0, 9.0, 9.0)

YtstJmp = 0.8*yLimMax/max(yLimMax) 

yLimMin = -c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8)*yLimMax/max(yLimMax)
YannMin = -c(1.9, 1.9, 1.9, 1.9, 1.9, 1.9)*yLimMax/max(yLimMax)

# Parameters for plotting
TstTxtSiz = 3 
AnnTxtSiz = 3
DstGrp = FALSE 
AnnTxtExt = "" 
YtstBia = rep(0,6)
TstSep = 0.08
TckAng = 30
Alp = 0.1

Yscr = "Number of Indirectly\nconnected ARCs"


### PLOTTING LOOP ##############################################################

TypIntArr = c("EndPin", "EndC95", "EndC90", "EndKgo")
GrpLbl = c("PPI", "COX.95", "COX.90", "KEGG")


MenFrm = data.frame()
x=1
TypIntLen = length(TypIntArr)
P = list()
for(x in 1:TypIntLen){
  
  print(x)
  QryTypInt = TypIntArr[x]
  CmpFrm = GenTypCmpFrm %>% filter(GenTyp %in% FltElm) %>% filter(TypInt == QryTypInt)
  
  grouped_means <-  CmpFrm %>%
    group_by(GenTyp) %>%
    summarise(mean_Y = mean(GenNg1, na.rm = TRUE))
  grouped_means
  
  grouped_medians <-  CmpFrm %>%
    group_by(GenTyp) %>%
    summarise(median_Y = median(GenNg1, na.rm = TRUE))
  
  sMenFrm = grouped_means %>% merge(grouped_medians, by='GenTyp')
  sMenFrm$Network = QryTypInt
  MenFrm = rbind(MenFrm,sMenFrm)
  

  # NUMBER OF ELEMENTS
  AgeNum = CmpFrm %>% filter(GenTyp == "AgeGen") %>% nrow()
  ModNum = CmpFrm %>% filter(GenTyp == "ModGen") %>% nrow()
  DisNum = CmpFrm %>% filter(GenTyp == "DisGen") %>% nrow()
  NgbNum = CmpFrm %>% filter(GenTyp == "DisNgb") %>% nrow()
  BeyNum = CmpFrm %>% filter(GenTyp == "aOthGen") %>% nrow()
  NumArr = c(DisNum, AgeNum, ModNum, NgbNum, BeyNum)
  
  
  # FILTER DESIRED FRAME
  TstFrm = CmpFrm %>% select(.data[[Xvar]], .data[[Yvar]], TypInt)
  colnames(TstFrm) = c("X", "Y", "G")
  TstFrm$Y = TstFrm$Y 
  TstFrm = TstFrm %>% filter(!is.na(X))
  TstFrm$G = GrpLbl[x]
  
  pwc <- TstFrm %>% pairwise_t_test( Y ~ X , p.adjust.method = "bonferroni")
  
  GrpLen_Pwc = c(pwc$group1, pwc$group2) %>% unique() %>% length()
  GrpLen_Org = TstFrm$X %>% unique() %>% length()
  
  yMaxLoc = max(TstFrm$Y,na.rm=TRUE)
  
  # IF NOT ALL PWC ROWS ARE PRESENT WHERE THEY SHOULD BE
  if(GrpLen_Pwc < GrpLen_Org){
    
    Xarr = TstFrm$X %>% unique() %>% sort()
    Xprw <- consecutive_pairs(Xarr)
    
    # COMBINATORY
    PrwFrm = abs (combn(Xarr, 2) - max(Xarr) ) %>% as.data.frame()
    SrtPrwFm = PrwFrm[length(PrwFrm):1]
    Xprw = apply(SrtPrwFm, 2, list)
    
    Xlen = length(Xprw)
    pwc = data.frame()
    GrpMaxArr = c()
    
    # COMPUTE ROWWISE TESTS
    for(i in 1:Xlen){
      QryGrpArr = Xprw[[i]] %>% unlist() %>% rev()
      sTstFrm = TstFrm %>% filter(X %in% QryGrpArr)
      GrpMaxArr[i] = sTstFrm %>% pull(Y) %>% max()
      spwc = sTstFrm %>% pairwise_t_test( Y ~ X , p.adjust.method = "bonferroni")
      # IF ROW DOESNT CONTAIN ANYTING
      if(nrow(spwc) == 0){
        npwc = pwc[i-1,]
        GrpLenArr = sTstFrm$X %>% table() %>% as.numeric()
        spwc = npwc %>% mutate(group1 = QryGrpArr[1], group2 = QryGrpArr[2], n1 = GrpLenArr[1], n2 = GrpLenArr[2], p=1, p.signif="ns", p.adj=1, p.adj.signif="ns")
      }
      pwc = rbind(pwc, spwc)
    }
    
  }
  
  
  # WILXON TEXT IF NECCESARY
  if(Yvar %in% c("GenNg0", "GenNg1", "GenNg10") & WixTst==TRUE){
    pwc2 <- pairwise.wilcox.test( TstFrm$Y , TstFrm$X , p.adjust.method = "bonferroni") %>% tidy() %>% mutate(p.adj=p.value, P=p.value, `.y.`="Y") 
    if(nrow(pwc)==0){
      pwc = pwc2
    }
    pwc$p = pwc2$p.value
    pwc$p.adj = pwc2$p.value %>% signif(2)
    CapTxt = "pwc: Wilcox test; p.adjust: Bonferroni"
    #Yscr = "Indirect pleitropy"
    #Yscr = "Number of indirectly\nconnected ARCs"
  }
  
  # - X,Y POSITIONS OF STAT SYMBOLS --------------------------------------------
  pwc = pwc %>% add_xy_position(x = "X")
  PrePwc = pwc
  NaPwx = sum(is.na(PrePwc)) > 0
  

  # IF NA DATA IN X,Y POSITIONS
  if(NaPwx){
    PrePwc$y.position = GrpMaxArr
    PwcLen = nrow(PrePwc)
    for(i in 1:PwcLen){
      PrePwc[i,][['groups']][[1]] = Xprw[[i]] %>% unlist() %>% rev() %>% as.character()
    }
    names(PrePwc$groups) = paste("V",1:PwcLen,sep="")
    
  }
  
  
  PltPwc = PrePwc 
  PltPwc$xmin = FltElmGr1[PltPwc$group1,"xmin"]
  PltPwc$xmax = FltElmGr2[PltPwc$group2,"xmax"]
  PltPwc <- PltPwc[order(PltPwc$xmin, PltPwc$xmax), ]
  PltPwc$y.position = sort(PltPwc$y.position)
  
  
  MinArr = PltPwc$xmin
  MaxArr = PltPwc$xmax
  MinMaxFrm = data.frame(MinArr, MaxArr)
  PltPwc$xmax = rowMaxs(MinMaxFrm)
  PltPwc$xmin = rowMins(MinMaxFrm)
  
  
  yPosArr = sort(PltPwc$y.position)
  FstPltPwc = PltPwc %>% filter(y.position == yPosArr[1]) 
  yPosMin = yPosArr[1]
  
  if(FstPltPwc$xmax-FstPltPwc$xmin != 1){
    NewPos = PltPwc %>% filter(xmax == 1+xmin) %>% pull(y.position) %>% min()
    PltPwc[PltPwc$y.position==yPosMin, "y.position"] = NewPos
  }
  
  
  PltPwc = PltPwc %>% mutate(y.position = ifelse(abs(xmin-xmax)==1,yPosMin,y.position) )
  
  PltPwc = PltPwc %>% mutate(xmin = xmin+0.1, xmax = xmax-0.1)
  PltPwc = PltPwc[order(PltPwc$y.position), ]
  YposTbl = PltPwc$y.position %>% table() %>% names()
  YposTblLen = length(YposTbl)
  MinPnt = YtstMin[x]

  
  for(j in 1:YposTblLen){
    PltPwc[PltPwc$y.position == YposTbl[j], "y.position"] = MinPnt + ((j-1)*YtstJmp[x])
  }

  
  # STATS ASTHERISKS
  PsgnArr = c("ns", "ns", "ns")
  PsgnArr[PltPwc$p > 5e-2] = "ns"
  PsgnArr[5e-2 >= PltPwc$p & PltPwc$p > 1e-2] = "*"
  PsgnArr[1e-2 >= PltPwc$p & PltPwc$p > 1e-3] = "**"
  PsgnArr[1e-3 >= PltPwc$p & PltPwc$p > 1e-4] = "***"
  PsgnArr[1e-4 >= PltPwc$p] = "****"
  
  PadjArr = c("ns", "ns", "ns")
  PadjArr[PltPwc$p.adj > 5e-2] = "ns"
  PadjArr[5e-2 >= PltPwc$p.adj & PltPwc$p.adj > 1e-2] = "*"
  PadjArr[1e-2 >= PltPwc$p.adj & PltPwc$p.adj > 1e-3] = "**"
  PadjArr[1e-3 >= PltPwc$p.adj & PltPwc$p.adj > 1e-4] = "***"
  PadjArr[1e-4 >= PltPwc$p.adj] = "****"
  
  PltPwc$p.signif = PsgnArr
  PltPwc$p.adj.signif = PadjArr
  
  
  if(class(TstFrm$X) == "numeric"){
    ElmOrd = TstFrm$X %>% unique() %>% sort()
    xTckLab = ElmOrd
  } else{
    ElmOrd = FltElm
  }
  
  ElmLen = length(ElmOrd)
  
  
  # - PLOT BOX -----------------------------------------------------------------
  p <- ggboxplot(TstFrm, x = "X", y = "Y",
                 fill='#F8776D',
                 facet.by = "G",
                 combine = TRUE,
                 order = ElmOrd,
                 outlier.size = 0.1
  ) 
  
  
  p = p + coord_cartesian(ylim = c(yLimMin[x],yLimMax[x]))
  
  p = p+stat_pvalue_manual(PltPwc, label = "p.adj.signif", tip.length = 0, step.increase = 0.00, label.size = TstTxtSiz) 
  
  p = p + grids(axis = c("xy"), linetype = "dashed",  color = "grey", size = NULL) #x, y, xy
  
  p = p +
    annotate("text",
             x = 1:ElmLen,
             y = YannMin[x],
             label = NumArr %>% format_with_commas(), 
             col = "black",
             size = AnnTxtSiz,
             vjust = - 1)
  
  p =  p + scale_x_discrete(labels=xTckLab)

  p = p + theme(axis.text.x = element_text(angle = TckAng, vjust = 1.2, hjust=1)) 
  
  p = p + scale_y_continuous(breaks = c(0,2,4,6,8,10,12))  

  # MEAN SYMBOL
  p = p + stat_summary(fun = mean, geom = "point", shape = "*", 
               size = 6, color = "black", position = position_dodge(0.8))
  
  
  p = ggpar(p, xlab ="Gene Type", ylab = Yscr, legend.title = "") +
    font("title", size = 10, color = "black", face = "bold") +
    font("xlab", size = 9, color = "black", face = "bold")+
    font("ylab", size = 9, color = "black", face = "bold")+
    font("xy.text", size = 8, color = "black") +
    theme(legend.position = "none") #+
  
  
  if(x %in% c(1,2)){
    p = p + xlab(NULL)
  }
  
  
  if(x %in% c(2,4,5,6)){
    p = p + ylab(NULL)
  }
  
  
  P[[QryTypInt]] = p
  
}


plot_grid(P$EndPin, P$EndC95, P$EndC90, P$EndKgo, labels = c("a","b","c","d"), label_size = 12, ncol=2)


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

