#!/usr/bin/env Rscript

HeterRegion=read.table("../RH_RH1015_heter_stat_50k_0.02_250k_merge.bed")
AllCutOffs=c(2.75,3.5,2)

all.RH_phased=NULL
all.C151_phased=NULL
all.E8669_phased=NULL
for (ConservationCutoff in AllCutOffs){
    allchrs=c("01","02","03","04","05","06","07","08","09","10","11","12")
    for (chr in allchrs){
        Chrom=paste0("chr",chr)
        RH_phased_one=read.table(paste0("../01_DeleInfor/JustPhasedDele/PhasedHap_JustPhasedDele_RH_",Chrom,"_GERP",ConservationCutoff,".txt"),header=T)
        RH_phased_one$ConservationCutoff=ConservationCutoff
        all.RH_phased=rbind(all.RH_phased,RH_phased_one)

        C151_phased_one=read.table(paste0("../01_DeleInfor/JustPhasedDele/PhasedHap_JustPhasedDele_C151_",Chrom,"_GERP",ConservationCutoff,".txt"),header=T)
        C151_phased_one$ConservationCutoff=ConservationCutoff
        all.C151_phased=rbind(all.RH_phased,RH_phased_one)

        E8669_phased_one=read.table(paste0("../01_DeleInfor/JustPhasedDele/PhasedHap_JustPhasedDele_E8669_",Chrom,"_GERP",ConservationCutoff,".txt"),header=T)
        E8669_phased_one$ConservationCutoff=ConservationCutoff
        all.E8669_phased=rbind(all.E8669_phased,E8669_phased_one)
    }
}

for (ConservationCutoff in AllCutOffs){

for (n in 1:nrow(HeterRegion)){
    Chr=HeterRegion[n,1]
    Start=HeterRegion[n,2]
    End=HeterRegion[n,3]
    locus=paste(Chr,Start,End,sep='_')
    figure_name=paste0("E8669_C151_RH_onlyPhased_",locus,"_GERP",ConservationCutoff,".pdf")

 h1 <- all.RH_phased[all.RH_phased$position>=Start & all.RH_phased$position<=End & all.RH_phased$Hap1==2,]
  h1$pos= h1$position-Start
    h2 <- all.RH_phased[all.RH_phased$position>=Start & all.RH_phased$position<=End & all.RH_phased$Hap2==2,]
  h2$pos= h2$position-Start
    h3 <- all.C151_phased[all.C151_phased$position>=Start & all.C151_phased$position<=End & all.C151_phased$Hap1==2,]
    h4 <- all.C151_phased[all.C151_phased$position>=Start & all.C151_phased$position<=End & all.C151_phased$Hap2==2,]
  h3$pos= h3$position-Start
  h4$pos= h4$position-Start
    h5 <- all.E8669_phased[all.E8669_phased$position>=Start & all.E8669_phased$position<=End & all.E8669_phased$Hap1==2,]
    h6 <- all.E8669_phased[all.E8669_phased$position>=Start & all.E8669_phased$position<=End & all.E8669_phased$Hap2==2,]
  h5$pos= h5$position-Start
  h6$pos= h6$position-Start

  pdf(figure_name,10,4)
par(mgp=c(2,0.75,0),oma=c(1,1,1,1),mar=c(3,1,1,1),las=1,cex.axis=1.15,cex.lab=1.3,cex.main=1.7)
chr_len <- End-Start+1
plot(1,1, type="n",xlab="Chromosome (Mb)",ylab="",axes=F, xlim=c(0, chr_len), ylim=c(0.25,5))
axis(1, at=seq(0, chr_len, 100000), label=seq(0, chr_len, 100000), tcl=-0.4,lwd=1.5,lwd.ticks=1.5)
for (i in h1$pos){
    segments(i, 1, i, 1.3, lwd=1.5, col=rgb(129,196,240,max=255))
}

for (i in h2$pos){
    segments(i, 0.5, i, 0.8, lwd=1.5, col=rgb(235,70,144,max=255))
}

rect(0,0.5,chr_len,0.8,border='black',lwd=1.5,density=0)
rect(0,1,chr_len,1.3,border='black',lwd=1.5,density=0)

for (i in h3$pos){
    segments(i, 2.5, i, 2.8, lwd=1.5, col=rgb(129,196,240,max=255))
}

for (i in h4$pos){
    segments(i, 2, i, 2.3, lwd=1.5, col=rgb(235,70,144,max=255))
}

rect(0,2,chr_len,2.3,border='black',lwd=1.5,density=0)
rect(0,2.5,chr_len,2.8,border='black',lwd=1.5,density=0)


for (i in h5$pos){
    segments(i, 4, i, 4.3, lwd=1.5, col=rgb(129,196,240,max=255))
}

for (i in h6$pos){
    segments(i, 3.5, i, 3.8, lwd=1.5, col=rgb(235,70,144,max=255))
}

rect(0,3.5,chr_len,3.8,border='black',lwd=1.5,density=0)
rect(0,4,chr_len,4.3,border='black',lwd=1.5,density=0)
  dev.off()
}
}
