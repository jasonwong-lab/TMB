library(GenomicRanges)
library(VariantAnnotation)

# load data
load("R/sysdata.rda")
# all coding region
exome.bed<-exome.final.bed
gr.exome<-GRanges(seqnames = Rle(exome.bed$V1),ranges = IRanges(exome.bed$V2,exome.bed$V3))

# main function
#' TMBpredict
#'
#' Correct panel TMB evaluation by modeling adjustment.
#' @param ttype Tumor type
#' @param mut Detected mutations in VCF format. Also accept multiple VCF compressed in .tar.gz.
#' @param panel.bed Panel region file in BED format
#' @param convert.type Convert Panel to WES or reversely
#' @return The adjusted TMB value and correlation figure
#' @examples 
#' TMBpredict("COAD","sample.vcf","panel.bed","panel2wes")
#' @export


TMBpredict<-function(ttype,mut,panel.bed,convert.type){
  
  #==============  
  # WES to Panel
  #=============
  if(convert.type=="panel2wes"){
  #panel bed
  panel.bed<-read.table(panel.bed)
  gr.panel<-GRanges(seqnames = Rle(panel.bed$V1),ranges = IRanges(panel.bed$V2,panel.bed$V3))
  tmp<-findOverlapPairs(gr.panel,gr.exome)
  gr.panel<-pintersect(tmp)
  
  # calculate tcga muations in given panel (all mutations)
  mut.exome.bed<-all.mut.list[[ttype]]
  gr.exome.mut<-GRanges(seqnames = Rle(mut.exome.bed$V1),ranges = IRanges(mut.exome.bed$V2,mut.exome.bed$V2),sample=mut.exome.bed$V4)
  tmp<-subsetByOverlaps(gr.exome.mut,gr.panel)
  panel.mut.cout<-as.data.frame(table(mcols(tmp)$sample))
  rownames(panel.mut.cout)<-panel.mut.cout$Var1
  
  # caluculate tcga muations in all coding region (non-syn mutations)
  exome.mut.cout<-wes.tmb.list[[ttype]]
  
  # prepare correlation data
  sname<-intersect(rownames(panel.mut.cout),rownames(exome.mut.cout))
  final.data<-data.frame(panel.mut.cout[sname,]$Freq,as.numeric(exome.mut.cout[sname,]))
  colnames(final.data)<-c("panel","wes")
  rownames(final.data)<-sname
  
  # model estimate
  region.exome=36.747178
  tmp<-data.frame(x=gr.panel)
  region.panel=sum(tmp[,3]-tmp[,2])/1000000
  x<-final.data$panel/region.panel
  y<-final.data$wes/region.exome
  fit<-lm(y~x)
  max<-max(x,y)
  
  # plot correlation of panel TMB and WES TMB for TCGA data
  #dev.new()
  pdf("TMB_correlation.pdf",4,4)
  par(mar=c(4,4,1,1),mgp=c(2,.5,0))
  plot(x,y,xlim=c(0,max),ylim=c(0,max),xlab="Panel (mut/Mb)",ylab="WES (mut/Mb)")
  fit<-lm(y~x)
  abline(fit,lty=2,col="#fc8d62")
  p<-signif(summary(fit)$coefficients[,4][2],5)
  tmp<-cor.test(x,y)
  r<-signif(tmp$estimate,3)
  mtext(paste0(ttype," (R=",r,", P=",p,")"),3,cex=1,line=-1)
  
  
  if(length(grep(".tar",mut,fixed=T))!=1){
    # upload mutations with vcf format
    vcf <- readVcf(mut, "hg19")
    gr.panel.mut<-rowRanges(vcf)
    seqlevelsStyle(gr.panel.mut) <- "UCSC"
  
    obs.panel<-sum(countOverlaps(gr.panel.mut, gr.panel))/region.panel
    z<-predict(fit,newdata=data.frame(x=obs.panel))
    if(z<0){z<-"Too few mut to estimate"}
  
    # output results
    write.out<-data.frame(PANEL=obs.panel,Predicted_WES=z)
    colnames(write.out)<-c("Observed mutations (mut/Mb)","Predicted WES TMB (mut/Mb)")
    sap.id<-gsub(".vcf","",basename(mut))
    rownames(write.out)<-sap.id
    write.table(write.out,"TMB_convert.txt",sep = "\t",quote = F)
    points(write.out[,1],write.out[,2],col="#e41a1c",pch=16)
    dev.off()
    packageStartupMessage("The conversion is successfully finished and the outputs are stored in working directory (TMB_correlation.pdf and TMB_convert.txt).")
    return(write.out)
  } else {
    # upload mutations with tar.gz format
    sap.list<-untar(mut,list=TRUE)
    untar(mut)
    n.sap<-length(sap.list)
    
    write.out<-data.frame(PANEL=rep(NA,n.sap),Predicted_WES=rep(NA,n.sap))
    colnames(write.out)<-c("Observed mutations (mut/Mb)","Predicted WES TMB (mut/Mb)")
    rownames(write.out)<-gsub(".vcf","",basename(sap.list))
    
    for (j in 1:n.sap){
      vcf <- readVcf(sap.list[j], "hg19")
      gr.panel.mut<-rowRanges(vcf)
      seqlevelsStyle(gr.panel.mut) <- "UCSC"
    
      # model estimate
      obs.panel<-sum(countOverlaps(gr.panel.mut, gr.panel))/region.panel
      z<-predict(fit,newdata=data.frame(x=obs.panel))
      if(z<0){z<-"Too few mut to estimate"}
 
      # output results
      write.out[j,1]<-obs.panel
      write.out[j,2]<-z
      }
    write.table(write.out,"TMB_convert.txt",sep = "\t",quote = F)
    points(write.out[,1],write.out[,2],col="#e41a1c",pch=16)
    dev.off()
    packageStartupMessage("The convertion is successfully finished and the outputs are stored in working directory (TMB_correlation.pdf and TMB_convert.txt).")
    return(write.out)
  }
}

#==============  
# WES to Panel
#=============
  if(convert.type=="wes2panel"){
    #panel bed
    panel.bed<-read.table(panel.bed)
    gr.panel<-GRanges(seqnames = Rle(panel.bed$V1),ranges = IRanges(panel.bed$V2,panel.bed$V3))
    tmp<-findOverlapPairs(gr.panel,gr.exome)
    gr.panel<-pintersect(tmp)
    
    # calculate tcga muations in given panel or exome(all mutations)
    mut.exome.bed<-all.mut.list[[ttype]]
    gr.exome.mut<-GRanges(seqnames = Rle(mut.exome.bed$V1),ranges = IRanges(mut.exome.bed$V2,mut.exome.bed$V2),sample=mut.exome.bed$V4)
    tmp<-subsetByOverlaps(gr.exome.mut,gr.exome)
    panel.mut.cout<-as.data.frame(table(mcols(tmp)$sample))
    rownames(panel.mut.cout)<-panel.mut.cout$Var1
    
    # caluculate tcga muations in all coding region (non-syn mutations)
    mut.exome.bed1<-ns.mut.list[[ttype]]
    gr.exome.mut1<-GRanges(seqnames = Rle(mut.exome.bed1$V1),ranges = IRanges(mut.exome.bed1$V2,mut.exome.bed1$V2),sample=mut.exome.bed1$V4)
    tmp1<-subsetByOverlaps(gr.exome.mut1,gr.f1cdx)
    exome.mut.cout<-as.data.frame(table(mcols(tmp1)$sample))
    rownames(exome.mut.cout)<-exome.mut.cout$Var1
    
    # prepare correlation data
    sname<-intersect(rownames(panel.mut.cout),rownames(exome.mut.cout))
    final.data<-data.frame(panel.mut.cout[sname,]$Freq,exome.mut.cout[sname,]$Freq)
    colnames(final.data)<-c("panel","wes")
    rownames(final.data)<-sname
    
    # model estimate
    region.exome=36.747178
    tmp<-data.frame(x=gr.panel)
    region.panel=sum(tmp[,3]-tmp[,2])/1000000
    x<-final.data$panel/region.exome
    y<-final.data$wes/region.panel
    fit<-lm(y~x)
    max<-max(x,y)
    
    # plot correlation of panel TMB and WES TMB for TCGA data
    #dev.new()
    pdf("TMB_correlation.pdf",4,4)
    par(mar=c(4,4,1,1),mgp=c(2,.5,0))
    plot(x,y,xlim=c(0,max),ylim=c(0,max),xlab="WES (mut/Mb)",ylab="Panel (mut/Mb)")
    fit<-lm(y~x)
    abline(fit,lty=2,col="#fc8d62")
    p<-signif(summary(fit)$coefficients[,4][2],5)
    tmp<-cor.test(x,y)
    r<-signif(tmp$estimate,3)
    mtext(paste0(ttype," (R=",r,", P=",p,")"),3,cex=1,line=-1)
    
    
    if(length(grep(".tar",mut,fixed=T))!=1){
      # upload mutations with vcf format
      vcf <- readVcf(mut, "hg19")
      gr.panel.mut<-rowRanges(vcf)
      seqlevelsStyle(gr.panel.mut) <- "UCSC"
      
      obs.panel<-sum(countOverlaps(gr.panel.mut, gr.panel))/region.exome
      z<-predict(fit,newdata=data.frame(x=obs.panel))
      if(z<0){z<-"Too few mut to estimate"}
      
      # output results
      write.out<-data.frame(PANEL=obs.panel,Predicted_WES=z)
      colnames(write.out)<-c("Observed mutations (mut/Mb)","Predicted Panel TMB (mut/Mb)")
      sap.id<-gsub(".vcf","",basename(mut))
      rownames(write.out)<-sap.id
      write.table(write.out,"TMB_convert.txt",sep = "\t",quote = F)
      points(write.out[,1],write.out[,2],col="#e41a1c",pch=16)
      dev.off()
      packageStartupMessage("The conversion is successfully finished and the outputs are stored in working directory (TMB_correlation.pdf and TMB_convert.txt).")
      return(write.out)
    } else {
      # upload mutations with tar.gz format
      sap.list<-untar(mut,list=TRUE)
      untar(mut)
      n.sap<-length(sap.list)
      
      write.out<-data.frame(PANEL=rep(NA,n.sap),Predicted_WES=rep(NA,n.sap))
      colnames(write.out)<-c("Observed mutations (mut/Mb)","Predicted Panel TMB (mut/Mb)")
      rownames(write.out)<-gsub(".vcf","",basename(sap.list))
      
      for (j in 1:n.sap){
        vcf <- readVcf(sap.list[j], "hg19")
        gr.panel.mut<-rowRanges(vcf)
        seqlevelsStyle(gr.panel.mut) <- "UCSC"
        
        # model estimate
        obs.panel<-sum(countOverlaps(gr.panel.mut, gr.panel))/region.exome
        z<-predict(fit,newdata=data.frame(x=obs.panel))
        if(z<0){z<-"Too few mut to estimate"}
        
        # output results
        write.out[j,1]<-obs.panel
        write.out[j,2]<-z
      }
      write.table(write.out,"TMB_convert.txt",sep = "\t",quote = F)
      points(write.out[,1],write.out[,2],col="#e41a1c",pch=16)
      dev.off()
      packageStartupMessage("The conversion is successfully finished and the outputs are stored in working directory (TMB_correlation.pdf and TMB_convert.txt).")
      return(write.out)
    }
  }
# wes2panel end
}
