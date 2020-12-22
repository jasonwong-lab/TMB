
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
#' @param mut Detected mutations in VCF format
#' @param panel.bed Panel region file in BED format
#' @return The adjusted TMB value and correlation figure
#' @examples 
#' tmb_pre("COAD","sample.vcf","panel.bed")
#' @export
tmb_pre<-function(ttype,mut,panel.bed){
  #panel bed
  panel.bed<-read.table(panel.bed)
  gr.panel<-GRanges(seqnames = Rle(panel.bed$V1),ranges = IRanges(panel.bed$V2,panel.bed$V3))
  tmp<-findOverlapPairs(gr.panel,gr.exome)
  gr.panel<-pintersect(tmp)
  
  # calculate tcga muations in given panel (all mutations)
  mut.exome.bed<-all.mut.list[[ttype]]
  gr.exome.mut<-GRanges(seqnames = Rle(mut.exome.bed$V1),ranges = IRanges(mut.exome.bed$V2,mut.exome.bed$V2),sample=mut.exome.bed$V6)
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
  
  # upload mutations with vcf format
  vcf <- readVcf(mut, "hg19")
  gr.panel.mut<-rowRanges(vcf)
  seqlevelsStyle(gr.panel.mut) <- "UCSC"
  
  # model estimate
  region.exome=36.747178
  tmp<-data.frame(x=gr.panel)
  region.panel=sum(tmp[,3]-tmp[,2])/1000000
  x<-final.data$panel/region.panel
  y<-final.data$wes/region.exome
  fit<-lm(y~x)
  max<-max(x,y)
  obs.panel<-sum(countOverlaps(gr.panel.mut, gr.panel))/region.panel
  z<-predict(fit,newdata=data.frame(x=obs.panel))
  if(z<0){z<-"Too few mut to estimate"}
  
  # plot correlation of panel TMB and WES TMB for TCGA data
  dev.new()
  pdf("TMB_correlation.pdf",4,4)
  par(mar=c(4,4,1,1),mgp=c(2,.5,0))
  plot(x,y,xlim=c(0,max),ylim=c(0,max),xlab="Panel (mut/Mb)",ylab="WES (mut/Mb)")
  fit<-lm(y~x)
  abline(fit,lty=2,col="#fc8d62")
  p<-signif(summary(fit)$coefficients[,4][2],5)
  tmp<-cor.test(x,y)
  r<-signif(tmp$estimate,3)
  mtext(paste0(ttype," (R=",r,", P=",p,")"),3,cex=1,line=-1)
  dev.off()
  
  # output results
  write.out<-data.frame(PANEL=obs.panel,Predicted_WES=z)
  colnames(write.out)<-c("Observed mutations","Predicted TMB")
  write.out[1,]<-paste(round(write.out[1,],3),"mut/Mb")
  write.table(write.out,"TMB_predicted_WES.txt",sep = "\t",quote = F,row.names = F)
}

