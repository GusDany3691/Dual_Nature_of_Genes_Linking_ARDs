library(ggplot2)
library(dplyr)
library(ggrepel)

################################################################################
# LOADING SPECIFICITY AND COEXPRESSION DATA 
################################################################################

CoxFrm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Coexpression/BoxComScrFrm.csv") %>%
  filter(!(Com %in% c("Plt1","Plt2","Plt3","Plt4","Plt5","Plt6")))
TauFrm = read.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/BoxComScrFrm.csv") %>%
  filter(!(Com %in% c("Plt1","Plt2","Plt3","Plt4","Plt5","Plt6")))

ComArr = CoxFrm$Com %>% unique()
ComLen = length(ComArr)
StsCox = data.frame()
StsTau = data.frame()
i=1
for(i in 1:ComLen){
  print(i)
  ComQry = ComArr[i]
  CoxArr = CoxFrm %>% filter(Com %in% ComQry) %>% pull(Scr)
  TauArr = TauFrm %>% filter(Com %in% ComQry) %>% pull(Scr)
  sStsCox = calculateStats(CoxArr, ComQry)
  sStsTau = calculateStats(TauArr, ComQry)
  StsCox = rbind(StsCox,sStsCox)
  StsTau = rbind(StsTau,sStsTau)
}

colnames(StsCox) = c("Com", "Men.Cox", "Med.Cox", "Q1.Cox", "Q3.Cox", "LowCnf.Cox", "UppCnf.Cox")
colnames(StsTau) = c("Com", "Men.Tau", "Med.Tau", "Q1.Tau", "Q3.Tau", "LowCnf.Tau", "UppCnf.Tau")


ColFrm = data.frame(Com = c("immunological/systemic disorders", "cardiovascular", "endocrine/diabetes", "gastrointestinal/abdominal", "musculoskeletal/trauma",       
                            "haematology/dermatology", "neurology/eye/psychiatry", "renal/urology", "GenAge.Hum", "GenAge.Mod"),
                    Col = c("Immunological", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "GenAge", "GenAge"))

StsCoxTau = merge(StsCox, StsTau, by="Com") %>% merge(ColFrm, by="Com")


################################################################################
# PLOT 
################################################################################

RedCol = "#F8696B"
ModCol = "#FF7F50"
BluCol = "#5A8AC6"
GryCol = "#CAE2EE" 

LinAlp = 0.7
NodAlp = 0.8

BigSiz = 2.5
SmlSiz = 1

WdtCIv = 0.010
WdtCIh = 0.010

WdtIQRv = 0.040
WdtIQRh = 0.065

# 3. Plot the data using ggplot2 and ggrepel
p <- ggplot(StsCoxTau, aes(x = Med.Tau, y = Med.Cox, label = Com)) +
  
  # For the error bars of "Other" group
  geom_errorbar(data = subset(StsCoxTau, Col == "Other"), 
                aes(ymin = LowCnf.Cox, ymax = UppCnf.Cox, linetype="CI"), width = WdtCIv, size = SmlSiz, color = GryCol, alpha = LinAlp) +
  geom_errorbar(data = subset(StsCoxTau, Col == "Other"), 
                aes(ymin = Q1.Cox, ymax = Q3.Cox, linetype="IQR"), width = WdtIQRv, size = BigSiz, color = GryCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "Other"), 
                 aes(xmin = LowCnf.Tau, xmax = UppCnf.Tau, linetype="CI"), height = WdtCIh, size = SmlSiz, color = GryCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "Other"), 
                 aes(xmin = Q1.Tau, xmax = Q3.Tau, linetype="IQR"), height = WdtIQRh, size=BigSiz, color = GryCol, alpha = LinAlp) +
  
  
  
  # For the error bars of "Other" group
  geom_errorbar(data = subset(StsCoxTau, Col == "GenAge"), 
                aes(ymin = LowCnf.Cox, ymax = UppCnf.Cox, linetype="CI"), width = WdtCIv, color = RedCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbar(data = subset(StsCoxTau, Col == "GenAge"), 
                aes(ymin = Q1.Cox, ymax = Q3.Cox, linetype="IQR"), width = WdtIQRv, size = BigSiz, color = RedCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "GenAge"), 
                 aes(xmin = LowCnf.Tau, xmax = UppCnf.Tau, linetype="CI"), height = WdtCIh, color = RedCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "GenAge"), 
                 aes(xmin = Q1.Tau, xmax = Q3.Tau, linetype="IQR"), height = WdtIQRh, size=BigSiz, color = RedCol, alpha = LinAlp) +
  
  
  
  # For the error bars of "Other" group
  geom_errorbar(data = subset(StsCoxTau, Com == "GenAge.Mod"), 
                aes(ymin = LowCnf.Cox, ymax = UppCnf.Cox, linetype="CI"), width = WdtCIv, color = ModCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbar(data = subset(StsCoxTau, Com == "GenAge.Mod"), 
                aes(ymin = Q1.Cox, ymax = Q3.Cox, linetype="IQR"), width = WdtIQRv, size = BigSiz, color = ModCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Com == "GenAge.Mod"), 
                 aes(xmin = LowCnf.Tau, xmax = UppCnf.Tau, linetype="CI"), height = WdtCIh, color = ModCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Com == "GenAge.Mod"), 
                 aes(xmin = Q1.Tau, xmax = Q3.Tau, linetype="IQR"), height = WdtIQRh, size=BigSiz, color = ModCol, alpha = LinAlp) +
  
  
  
  # For the error bars of "Other" group
  geom_errorbar(data = subset(StsCoxTau, Col == "Immunological"), 
                aes(ymin = LowCnf.Cox, ymax = UppCnf.Cox, linetype="CI"), width = WdtCIv, color = BluCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbar(data = subset(StsCoxTau, Col == "Immunological"), 
                aes(ymin = Q1.Cox, ymax = Q3.Cox, linetype="IQR"), width = WdtIQRv, size = BigSiz, color = BluCol, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "Immunological"), 
                 aes(xmin = LowCnf.Tau, xmax = UppCnf.Tau, linetype="CI"), height = WdtCIh, color = BluCol, size = SmlSiz, alpha = LinAlp) +
  geom_errorbarh(data = subset(StsCoxTau, Col == "Immunological"), 
                 aes(xmin = Q1.Tau, xmax = Q3.Tau, linetype="IQR"), height = WdtIQRh, size=BigSiz, color = BluCol, alpha = LinAlp) +
  
  geom_point(aes(color = Col), size = 4, alpha = NodAlp) +
  
 
  scale_color_manual(values = c("GenAge" = "red", "Immunological" = "blue", "ModAge" = "#FF4500", "Other" = "#67ADCF"), 
                     labels = c("GenAge.Hum\n(High iARC-Interactors\nin PPI and KEGG)\n", "Immunological\nDisorders\n(High ARC-Pleiotropy)\n",
                                "GenAge.Mod\n", "Other ARCs")) + #3E6C48, #398A50"
  
  
  scale_linetype_manual(values = c("IQR" = "solid", "CI" = "solid"), 
                        guide = guide_legend(override.aes = list(size = c(1.5, 1)))) +
  
  theme_minimal() +
  labs(title = "Specificity vs Self-Coexpression in GenAge and ARCs",
       x = "Specificity", y = "Self-Coexpression") +
  theme(panel.grid.major = element_line(linetype = 'dashed', color = "grey"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  labs(color = "Group Category") + #+
  
  guides(linetype = guide_legend(title = "Error Bar Type", 
                                 override.aes = list(color = "black", size = c(1.5, 1))))


print(p)


################################################################################
# FUNCTIONS
################################################################################

calculateStats <- function(array, com) {
  if (!is.vector(array) || is.null(array)) {
    stop("Input must be a non-null vector")
  }
  
  # Basic statistics
  stats <- data.frame(
    Com = com,
    Mean = mean(array, na.rm = TRUE),
    Median = median(array, na.rm = TRUE),
    Q1 = quantile(array, 0.25, na.rm = TRUE),
    Q3 = quantile(array, 0.75, na.rm = TRUE)
  )
  
  # Confidence values
  stats$Lower_Confidence_Value <- quantile(array, 0.05, na.rm = TRUE)
  stats$Upper_Confidence_Value <- quantile(array, 0.95, na.rm = TRUE)
  
  return(stats)
}

