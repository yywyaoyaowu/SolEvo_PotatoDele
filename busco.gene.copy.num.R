# Title     : TODO
# Objective : TODO
# Created by: wuyaoyao
# Created on: 1/28/21
library(tidyr)
library(ggplot2)
library(ggthemes)
all.buscp=NULL
species_list=read.table("01_genome.assembly/01_hifiasm_assembly/busco/species_list")
species=as.character(species_list[-c(2,5),])
for (n in 1:length(species)){
  one=species[n]
  print(one)
  busco=read.table(paste("01_genome.assembly/01_hifiasm_assembly/busco/",one,"_hifiasm.p_ctg_busco_full_table_at.least.one.cp.tsv",sep=""),sep="\t")
  busco.gene_num=table(table(busco[,1]))
  cp.num=names(busco.gene_num)
  busco_infor<-as.data.frame(cbind(busco.gene_num,cp.num))
  busco_infor$species=one
all.buscp=rbind(all.buscp,busco_infor)
}
all.buscp$busco.gene_num=as.numeric(all.buscp$busco.gene_num)
p=ggplot(all.buscp, aes(x=cp.num, y=busco.gene_num, fill=species)) + geom_bar(stat='identity',position="dodge")+ theme_set(theme_bw())
#ggplot(mask, aes(x=Region, y=mask.percentage, fill=method)) + geom_bar(stat='identity',position="dodge")+scale_fill_manual(values =colors)
p.axis=p + theme(axis.text.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.text.y = element_text(size = 15))
##change xlab and ylab
p.lab=p.axis + xlab("Cp num of Busco gene") + ylab("Num of Buscso gene")+ theme(axis.title.x = element_text(size = 15, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 15, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 15))+theme(legend.title = element_text(size = 15))
p.lenged
#+ theme(legend.position = "none")
#+ geom_hline(yintercept=2104350183, lty='dotted')
ggsave(p.lenged,filename ="01_genome.assembly/01_hifiasm_assembly/busco/Busco.gene.cp.num.pdf",height = 4,width = 7)
