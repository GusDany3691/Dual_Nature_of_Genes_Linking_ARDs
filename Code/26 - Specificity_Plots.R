library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(circlize)

# Srt - Sorted

################################################################################
# RETRIEVE VAIBALES
################################################################################

WidComComLst = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/WidComComLst.rds") 
ExtComMenFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/ExtComMenFrm.rds") 

### GENE SETS ##################################################################

ComArr = c("cardiovascular", "endocrine/diabetes", "gastrointestinal/abdominal", "renal/urology", 
           "haematology/dermatology", "immunological/systemic disorders", "musculoskeletal/trauma", "neurology/eye/psychiatry")

PltArr = c("Plt1", "Plt2", "Plt3", "Plt4", "Plt5", "Plt6")

AgeArr = c("GenAge", "ModAge")

ComPltAgeArr = c(ComArr, PltArr, AgeArr)
ComAgeArr = c(ComArr, AgeArr)
PltAgeArr = c(PltArr, AgeArr)

################################################################################
### BOXPLOT ####################################################################
################################################################################

ComComScrFrm = readRDS("C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/ComComScrFrm.rds") %>%
  filter(ExtCom == InnCom) %>% 
  mutate(Com = ExtCom) %>% 
  select(ExtCom, Scr) %>% 
  rename(Com=ExtCom)

TckAng = -30

# SORT COLUMNS BY MEDIAN
FrmComArr = ComComScrFrm$Com %>% unique()
FrmComLen = length(FrmComArr)
ComScrArr = c()
for(i in 1:FrmComLen){
  FrmComQry = FrmComArr[i]
  ComScrArr[i] = ComComScrFrm %>% filter(Com == FrmComQry) %>% pull(Scr) %>% median(na.rm=TRUE)
}
names(ComScrArr) = FrmComArr
SrtNam = ComScrArr %>% sort() %>% names() %>% rev()


### BOXPLOT ####################################################################

BoxComScrFrm = ComComScrFrm %>% filter(Com %in% ComAgeArr)

BoxComScrFrm$Com = ifelse(BoxComScrFrm$Com == "GenAge", "GenAge.Hum", BoxComScrFrm$Com)
BoxComScrFrm$Com = ifelse(BoxComScrFrm$Com == "ModAge", "GenAge.Mod", BoxComScrFrm$Com)

# Define custom color palette
custom_palette <- c("GenAge" = '#F8776D', "immunological/systemic disorders" = "#5A8AC6", "Other_ARCs" = "#CAE2EE")

# Add a new column for custom coloring
BoxComScrFrm$color_group <- ifelse(BoxComScrFrm$Com %in% c("GenAge.Hum", "GenAge.Mod"), "GenAge",
                                   ifelse(BoxComScrFrm$Com == "immunological/systemic disorders", "immunological/systemic disorders", "Other_ARCs"))


# ACTUAL PLOT
BoxComScrFrm$Lbl = "Tissue specificity in genes associated with ARCs and GenAge"
# Compute the counts for each box
counts <- as.data.frame(table(BoxComScrFrm$Com))

SrtNam = ifelse(SrtNam == "GenAge", "GenAge.Hum", SrtNam)
SrtNam = ifelse(SrtNam == "ModAge", "GenAge.Mod", SrtNam)


# Create the plot
p = ggboxplot(BoxComScrFrm, x = "Com", y = "Scr",
              outlier.size = 0.4,
              outlier.shape = NA,
              combine = TRUE,
              fill = "color_group",
              facet.by = "Lbl",
              order = SrtNam) + scale_fill_manual(values = custom_palette)

counts$FreqComma = format_with_commas(counts$Freq)

p = p+geom_dotplot(binaxis='y', stackdir='center',
                   position=position_dodge(1), dotsize=0.03, color="#202020")

#p = p + geom_text(data = counts, aes(x = Var1, y = -0.2, label = FreqComma), size=2.5)
p = p + geom_text(data = counts, aes(x = Var1, y = -0.0, label = FreqComma), size=3)

p = p + theme(axis.text.x = element_text(angle = TckAng, vjust = 0.6, hjust=0),
              strip.text = element_text(face="bold")) 

p = p + grids(axis = c("xy"), linetype = "dashed",  color = "grey", size = NULL)

p = ggpar(p, xlab ="ARC or GenAge Groups", ylab = "Tau", legend.title = "") +
  font("title", size = 10, color = "black", face = "bold") +
  font("xlab", size = 10, color = "black", face = "bold") +
  font("ylab", size = 10, color = "black", face = "bold") +
  font("xy.text", size = 9, color = "black") +
  theme(legend.position = "none")

# Ensure legend is visible
p = p + theme(legend.position = "top")  

p = p + theme(plot.margin = unit(c(0, 1, 0, 0), "cm"))  # Adjust plot margins

# Modify the x-axis labels
p = p + scale_x_discrete(labels = c(
  "immunological/systemic disorders" = "immunological/systemic disorders\n(High ARC-Pleiotropy genes)", 
  "haematology/dermatology" = "haematology/dermatology", 
  "endocrine/diabetes" = "endocrine/diabetes",
  "renal/urology" = "renal/urology",
  "gastrointestinal/abdominal" = "gastrointestinal/abdominal",
  "cancer" = "cancer",
  "neurology/eye/psychiatry" = "neurology/eye/psychiatry",
  "cardiovascular" = "cardiovascular",
  "musculoskeletal/trauma" = "musculoskeletal/trauma",
  "GenAge.Hum" = "GenAge.Hum",
  "GenAge.Mod" = "GenAge.Mod"
))

p

write.table(BoxComScrFrm,'C:/Users/Usuario/Desktop/Nature/Data/Generated/Specificity/BoxComScrFrm.csv', row.names = FALSE, sep=",")

################################################################################
### INTRA COMMUNITY DIFFERENCE (HEATMAP) #######################################
################################################################################

### INITIAL FRAMES AND PARAMETERS ##############################################

FntSiz = 8
TxtDeg = 60

PltFrm = WidComComLst$ExtInn$Dif
TxtFrm = WidComComLst$ExtInn$TstBnfTxt

RowNam = PltFrm$InnCom

row.names(PltFrm) = PltFrm$InnCom
row.names(TxtFrm) = TxtFrm$InnCom

PltFrm$InnCom = NULL
TxtFrm$InnCom = NULL


# SORTING ROWS AND COLUMNS 
TxtFrm = TxtFrm[ComPltAgeArr,ComPltAgeArr]
PltFrm = PltFrm[ComPltAgeArr,ComPltAgeArr]
ComMenFrm = data.frame(Men=ExtComMenFrm[row.names(PltFrm),])
row.names(ComMenFrm) = row.names(PltFrm)

# GENE GROUPS LABELLING
RowNam = row.names(PltFrm)
SptVec = rep("A",length(RowNam))
SptVec[RowNam %in% c("GenAge","ModAge")] = "H"
SptVec[RowNam %in% PltArr] = "S"
ExtComMenFrm = ExtComMenFrm[row.names(PltFrm),] %>% data.frame()
colnames(ExtComMenFrm) = "Men"
row.names(ExtComMenFrm) = row.names(PltFrm)

# TOP ANNOTATION OF HEATMAP ####################################################

TopAnn = HeatmapAnnotation(
  
  Coexpression_Mean = anno_simple(ComMenFrm$Men, 
                                  col = colorRamp2(c(min(ExtComMenFrm$Men, na.rm=TRUE), 
                                                     ( min(ExtComMenFrm$Men, na.rm=TRUE) + max(ExtComMenFrm$Men, na.rm=TRUE) )/2, 
                                                     max(ExtComMenFrm$Men, na.rm=TRUE)), 
                                                   c("#5A8AC6", "white", "#F8696B")),
                                  pch = as.character(ComMenFrm$Men),
                                  gp= gpar(col = "black"),
                                  height = unit(8, "mm"),,
                                  pt_size = unit(1, "pt")*6
  ),
  annotation_name_gp= gpar(fontsize = 8),
  
  text = anno_text(replace_dirplt(colnames(PltFrm)), rot = TxtDeg, gp = gpar(fontsize = FntSiz))
)


# HEATMAP ######################################################################

FntSiz = 8
TxtDeg = 60
Htm = Heatmap(as.matrix(PltFrm), name = "Differential Self-Coexpression\nColumn.Group - Row.Group", 
              column_title = "Column Group (ARC, Pleiotropy or GenAge Group)", 
              row_title = "Row Group (ARC, Pleiotropy or GenAge Group)",
              column_title_gp = gpar(fontsize = 12), 
              row_title_gp = gpar(fontsize = 12),     
              show_row_names = FALSE, 
              show_column_names = FALSE,
              border = "black",
              row_split =  SptVec,
              column_split =  SptVec,
              rect_gp = grid::gpar(col = "grey", lwd = 0.5), 
              
              bottom_annotation = TopAnn,
              
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(TxtFrm[i, j], x, y, gp = gpar(fontsize = FntSiz))
              }
) +
  Heatmap(ComMenFrm$Men,
          width = unit(9, "mm"),
          border = "black",
          rect_gp = grid::gpar(col = "black", lwd = 0.5), 
          name = "Self-Coexpression (Mean)",
          
          right_annotation = rowAnnotation(
            text = anno_text(replace_dirplt(row.names(PltFrm)), gp = gpar(fontsize = FntSiz))
          ),
          col = colorRamp2(c(min(ExtComMenFrm$Men, na.rm=TRUE), 
                             ( min(ExtComMenFrm$Men, na.rm=TRUE) + max(ExtComMenFrm$Men, na.rm=TRUE) )/2, 
                             max(ExtComMenFrm$Men, na.rm=TRUE)), 
                           c("#5A8AC6", "white", "#F8696B")),
          show_row_names = FALSE, show_column_names = FALSE,
          
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text( ComMenFrm[i, j], x, y, gp = gpar(fontsize = 7))
          },
          
          
  ) 


draw(Htm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

draw(Htm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     show_heatmap_legend = FALSE, show_annotation_legend = FALSE)


################################################################################
### FUNCTION
################################################################################

replace_dirplt <- function(arr) {
  # Use sub() function to replace all occurrences of "DirPlt" with "dPleiotropy"
  arr[arr=="GenAge"] = "GenAge.Hum" 
  arr[arr=="ModAge"] = "GenAge.Mod" 
  return(sub("Plt", "ARC-Pleiotropy_", arr))
}

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
