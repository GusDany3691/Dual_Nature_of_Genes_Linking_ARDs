#Nam.Keg = "GenNod.PhnNod_BRF_cv3o5i_a78g71sn67sp75_BthAnl2.csv"
#Nam.Pin = "GenNod.PhnNod_BRF_cv3o5i_a78g71sn67sp75_BthAnl2.csv"
library(scatterpie)
library(ggpubr)
library(cowplot)
library(ggpubr)

#Keg = paste("F:/UK_Biobank/Data/Generated/Sorted/MPA/Connectivity_analysis/DLHCP_dAcA_ThrAll/EndKgu/Pre/UadPinNgb/Inv/GenNod.PhnNod_BRF_cv3o5i_a78g71sn67sp75_BthAnl2.csv",sep="") %>% read.csv()
#Pin = read.csv("F:/UK_Biobank/Data/Generated/Sorted/MPA/Connectivity_analysis/DLHCP_dAcA_ThrAll/EndPin/Pre/UadPinNgb/Inv/GenNod.PhnNod_BRF_cv3o5i_a80g74sn75sp73_BthAnl2.csv")

#Keg = paste("D:/PhD/EndNtw/EndKgo/Pre/UadPinNgb/Inv/uGenNod.PhnNod_BRF_cv10o5i_a82g75sn75sp75_NoCancer_AllPlt.csv") %>% read.table(sep=";")
#Pin = paste("D:/PhD/EndNtw/EndPin/Pre/UadPinNgb/Inv/uGenNod.PhnNod_BRF_cv10o5i_a82g75sn76sp75_NoCancer_AllPlt.csv") %>% read.table(sep=";")

KegDirIndFrm = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Generated/EndKgo/DirIndIntFrm.rds') 
PinDirIndFrm = readRDS('C:/Users/Usuario/Desktop/Nature/Data/Generated/EndPin/DirIndIntFrm.rds') 

PinPre = read.table("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndPin/Pre/UadPinNgb/GenPhnCls_PrxFrm_BRF_cv10o5i_a81g75sn76sp74.csv", header = TRUE, sep=",")
PinPre$X = NULL
PinPre = PinPre %>% rename(Gen=Label)
Pin = merge(PinDirIndFrm,PinPre) %>% 
  select(Gen, Prob, Class, Test, DirIntCom, IndIntCom, DirIntPhn, IndIntPhn)

Pin %>% arrange(desc(Prob)) %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PredPPI_All.csv")
Pin %>% arrange(desc(Prob)) %>% filter(Class=="Age") %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PredPPI_Age.csv")
Pin %>% arrange(desc(Prob)) %>% filter(Class=="NotAge") %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PredPPI_NotAge.csv")


KegPre = read.table("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndKgo/Pre/UadPinNgb/GenPhnCls_PrxFrm_BRF_cv10o5i_a83g75sn75sp75.csv", header = TRUE, sep=",")
KegPre$X = NULL
KegPre = KegPre %>% rename(Gen=Label)
Keg = merge(KegDirIndFrm,KegPre) %>% 
  select(Gen, Prob, Class, Test, DirIntCom, IndIntCom, DirIntPhn, IndIntPhn)


Keg %>% arrange(desc(Prob)) %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PredKEGG_All.csv")
Keg %>% arrange(desc(Prob)) %>% filter(Class=="Age") %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PredKEGG_Age.csv")
Keg %>% arrange(desc(Prob)) %>% filter(Class=="NotAge") %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PredKEGG_NotAge.csv")



#Keg = paste("C:/Users/Usuario/Desktop/Nature/Data/Generated/EndKgo/Pre/UadPinNgb/GenPhnCls_PrxFrm_BRF_cv10o5i_a83g75sn75sp75.csv") %>% read.table(sep=",")
#Pin = paste("D:/PhD/EndNtw/EndPin/Pre/UadPinNgb/Inv/uGenNod.PhnNod_BRF_cv10o5i_a82g75sn76sp75_NoCancer_AllPlt.csv") %>% read.table(sep=";")


#Keg = Keg %>% select(Gen, Prob, Class, Test, UnqComInGenPt1, UnqComInGenPt2, PhnInGenPt1, PhnInGenPt2, GenAgeDst)
#Pin = Pin %>% select(Gen, Prob, UnqComInGenPt2, PhnInGenPt2, GenAgeDst)

#Keg = Keg %>% rename(Prb.Keg = Prob, UnqComInGenPt2.Keg=UnqComInGenPt2, PhnInGenPt2.Keg = PhnInGenPt2, GenAgeDst.Keg = GenAgeDst)
#Pin = Pin %>% rename(Prb.Pin = Prob, UnqComInGenPt2.Pin=UnqComInGenPt2, PhnInGenPt2.Pin = PhnInGenPt2, GenAgeDst.Pin = GenAgeDst)

#DirIntCom - UnqComInGenPt1
#IndIntCom - UnqComInGenPt2
#DirIntPhn - PhnInGenPt1
#IndIntPhn - PhnInGenPt2

Keg = Keg %>% select(Gen, Prob, Class, Test, DirIntCom, IndIntCom, DirIntPhn, IndIntPhn) %>%
  rename(Prb.Keg = Prob, IndIntCom.Keg=IndIntCom, IndIntPhn.Keg = IndIntPhn)

Pin = Pin %>% select(Gen, Prob, IndIntCom, IndIntPhn) %>%
  rename(Prb.Pin = Prob, IndIntCom.Pin=IndIntCom, IndIntPhn.Pin = IndIntPhn)

#Keg = Keg %>% rename(Prb.Keg = Prob, IndIntCom.Keg=IndIntCom, IndIntPhn.Keg = IndIntPhn)
#Pin = Pin %>% rename(Prb.Pin = Prob, IndIntCom.Pin=IndIntCom, IndIntPhn.Pin = IndIntPhn)


ShrGen = intersect(Keg$Gen, Pin$Gen)

SmlKeg = Keg %>% filter(Gen %in% ShrGen)
SmlKeg$Rnk.Keg = 1:nrow(SmlKeg)

SmlPin = Pin %>% filter(Gen %in% ShrGen)
SmlPin$Rnk.Pin = 1:nrow(SmlPin)

KegPin = merge(SmlKeg, SmlPin, by="Gen")
KegPin = KegPin %>% mutate(Rnk.Men = (Rnk.Keg+Rnk.Pin)/2)
KegPin = KegPin %>% mutate(Rnk.Gmn = sqrt(Rnk.Keg*Rnk.Pin))

KegPin = KegPin %>% mutate(Prb.Men = (Prb.Keg+Prb.Pin)/2)
KegPin = KegPin %>% mutate(Prb.Gmn = sqrt(Prb.Keg*Prb.Pin))
KegPin = KegPin %>% mutate(IndIntCom.Men = (IndIntCom.Keg+IndIntCom.Pin)/2 )
KegPin = KegPin %>% mutate(IndIntCom.Gmn = sqrt(IndIntCom.Keg+IndIntCom.Pin) )
KegPin = KegPin %>% mutate(IndIntPhn.Men = (IndIntPhn.Keg+IndIntPhn.Pin)/2 )
KegPin = KegPin %>% mutate(IndIntPhn.Gmn = sqrt(IndIntPhn.Keg+IndIntPhn.Pin) )
#KegPin = KegPin %>% mutate(GenAgeDst.Men = (GenAgeDst.Keg+GenAgeDst.Pin)/2 )
#KegPin = KegPin %>% mutate(GenAgeDst.Gmn = sqrt(GenAgeDst.Keg*GenAgeDst.Pin) )

KegPin = KegPin %>% arrange(desc(Prb.Men))

colnames(KegPin)

KegPin=
KegPin %>% select(Gen, Class, Prb.Men, Prb.Gmn, Prb.Keg, Prb.Pin, DirIntCom, DirIntPhn, 
                  IndIntCom.Keg, IndIntCom.Pin, IndIntCom.Men, IndIntCom.Gmn, 
                  IndIntPhn.Keg, IndIntPhn.Pin, IndIntPhn.Men, IndIntPhn.Gmn,
                  #GenAgeDst.Keg, GenAgeDst.Pin, GenAgeDst.Men, GenAgeDst.Gmn,
                  Rnk.Keg, Rnk.Pin, Rnk.Men, Rnk.Gmn, Test)



KegPin %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PairwisePred_All.csv")
KegPin %>% filter(Class=="Age") %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PairwisePred_Age.csv")
KegPin %>% filter(Class=="NotAge") %>% write.csv("C:/Users/Usuario/Desktop/Nature/Data/Generated/Pairwise_Predictions/PairwisePred_NotAge.csv")

#Keg$Rnk = 

cor(KegPin$IndIntCom.Keg, KegPin$IndIntCom.Pin)

cor(KegPin$Pred.x, KegPin$Pred.y)

