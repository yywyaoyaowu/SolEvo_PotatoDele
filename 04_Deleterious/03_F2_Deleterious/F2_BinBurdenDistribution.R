  # Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 2022/4/25
#input
setwd("/Users/wuyaoyao/work/03-evolutionary-dele/06_Sol_evo/03_Results/Paper/Output/NewDeleCutOff/F2")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette12 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                  "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
final_col=c("#C1CAC0","#DD8466","#FAE5BA","#96C9F7","#C9ACEC")
cbbPalette2=c("#56B4E9", "#C9ACEC")
cbbPalette3=cbbPalette[c(3,2,4)]
library(ggplot2)
library(tidyr)
library(dplyr)

F2_bin_burden_file="E463_A626_F2_2603BinBurden_GERP_AllDiffCutOffs.final_NewCutOff.txt"
cbbPalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

F2_bin_burden_all=read.table(F2_bin_burden_file,header=T)

AllCutOffs=c(2,2.75,3.5)
for (ConservationCutoff in AllCutOffs){
  print (ConservationCutoff)
F2_bin_burden_GERP2=subset(F2_bin_burden_all,Constraint_Cutoff==ConservationCutoff)
F2_bin_burden_GERP2_types=gather(F2_bin_burden_GERP2,DeleAllele,Burden,names(F2_bin_burden_GERP2[,9:10]))
p=ggplot(F2_bin_burden_GERP2_types, aes(x=Burden, color=DeleAllele))+ geom_histogram(fill="white", alpha=0.5, position="identity")+scale_colour_manual(values=cbbPalette[c(3,7)])+ theme_bw()
p.axis=p + theme(axis.text.x = element_text(size = 15))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("Bin deleterious burden") +ylab("Frequency") + theme(axis.title.x = element_text(size = 15))+ theme(axis.title.y = element_text(size = 15))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
ggsave(p.lenged,filename =paste0("F2_2603Bin.DeleAllele.hist_GERP",ConservationCutoff,".pdf"),height =4,width =7)
print(paste0("F2_2603Bin.DeleAllele.hist_GERP",ConservationCutoff,".pdf"))
}
