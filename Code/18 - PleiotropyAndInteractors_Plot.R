library(dplyr)
library(ggplot2)
library(viridis)
library(ggtext)
library(cowplot)


# Logaritmic parameters and colours
numbers <- 1:2125
log_values <- log10(ifelse(numbers == 0, 1, numbers))
normalized <- (log_values - min(log_values)) / (max(log_values) - min(log_values))
color_array <- viridis(length(unique(normalized)))[as.integer(normalized * (length(unique(normalized)) - 1)) + 1]
ColFrm <- data.frame(AllLen = numbers, Col = color_array)


TypIntArr = c("EndPin", "EndC95", "EndC90", "EndKgo")
TypTxtArr = c("PPI", "COX.95", "COX.90", "KEGG")
TypIntLen = length(TypIntArr)

pLst = list()
for(z in 1:TypIntLen){
  
  print(z)
  TypInt = TypIntArr[z]
  TypTxt = TypTxtArr[z]
  
  TxtPth = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/",TypInt,"/",sep="")
  
  NodLst = paste(TxtPth,"EndNod.rds",sep="") %>% readRDS()
  NtwGen = NodLst$GenPin
  
  PltOrg = paste(TxtPth,"DirIndIntFrm.rds",sep="") %>% readRDS() # DirIndIntFrm
  PltOrg = PltOrg %>% filter(Gen %in% NtwGen)
  
  PltOrg$ComCom = paste(PltOrg$DirIntCom,PltOrg$IndIntCom,sep="") 
  ComComArr = PltOrg$ComCom %>% unique()
  ComComLen = length(ComComArr)
  
  cor_value <- cor(PltOrg$DirIntCom, PltOrg$IndIntCom, method = "pearson")  # You can change the method if needed
  
  PieSctFrm = data.frame()
  GenLen = nrow(PltOrg)
  
  for(x in 1:ComComLen){
    ComComQry = ComComArr[x]
    PltOrgQry = PltOrg %>% filter(ComCom == ComComQry)
    AllLen = nrow(PltOrgQry)
    
    sPieSctFrm = data.frame(Dir = unique(PltOrgQry$DirIntCom), 
                            Ind = unique(PltOrgQry$IndIntCom), 
                            AllLen, 
                            AllLog = log10(AllLen),
                            Reg = ComComQry
    )
    PieSctFrm = rbind(PieSctFrm,sPieSctFrm)
    
  }
  
  PieSctFrm$AllPer = PieSctFrm$AllLen / GenLen

  
  ### SECOND PLOT ATTEMPT ########################################################################
  
  MaxSet = max(PieSctFrm$AllLog)
  PieSctFrm$Rad = 0.5*PieSctFrm$AllLog / MaxSet
  
  # - ACTUAL PLOT --------------------------------------------------------------------------------
  
  PieSctFrm$G = TypTxt
  
  PieSctFrm = merge(PieSctFrm, ColFrm, by="AllLen")
  
  
  # PlOT
  p <- ggplot(PieSctFrm, aes(x = Dir, y = Ind)) +
    geom_point(aes(size = 1, fill = AllLen), shape = 21, color = "black", stroke = 1) +
    scale_fill_gradientn(colors = PieSctFrm$Col, trans="log10", name="Gene Count", # I CAN END HERE WITH ")" + COMMENT THE FOLLOWING TWO LINES
                         guide = guide_colourbar(direction = "horizontal")) + theme(legend.position = "bottom") +
    guides(fill = "none")
  
  
  p = p +  scale_x_continuous(breaks = c(0,2,4,6)) +   # Ensure x-axis has ticks at every integer from 1 to 7
    scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
    facet_wrap(~ G) +
    theme(strip.text = element_text(face = "bold", size=10)) # This line makes the facet labels bold
  

  p = p + theme(
    axis.title.x = element_markdown(size = 10), 
    axis.title.y = element_markdown(size = 10)   
  )
  
  
  p = p + labs(x="**ARC-Pleiotropy**<br>(Directly connected ARCs)", y="**iARC-Interactions**<br>(Indirectly connected ARCs)")
  
  p = p + annotate("text", x = 6.5, y = 9, label = sprintf("Corr: %.2f", cor_value), hjust = 1, vjust=0, size=3.5)  # Adjust x and y for positioning
  
  
  if(z %in% c(2,4,6)){
    p = p + theme(axis.title.y = element_blank())
  }
  

  if(z %in% c(1,2)){
    p = p + theme(axis.title.x = element_blank())
  }
  
  p = p + coord_cartesian(ylim = c(-0.5,9.5))
  
  p = p + guides(size = guide_none())
  
  pLst[[TypInt]] = p
  
}


plot_grid(pLst$EndPin, pLst$EndC95, pLst$EndC90, 
          pLst$EndKgo,
          labels = "AUTO", label_size = 12, ncol=2)
  
