
# all coding region
exome.bed<-read.table("data/region_bed/exome.final.bed")
gr.exome<-GRanges(seqnames = Rle(exome.bed$V1),ranges = IRanges(exome.bed$V2,exome.bed$V3))

# upload panel region
f1cdx.bed<-read.table("www/msk_coding.bed")
gr.f1cdx<-GRanges(seqnames = Rle(f1cdx.bed$V1),ranges = IRanges(f1cdx.bed$V2,f1cdx.bed$V3))

# restrict panel within coding regions
tmp<-findOverlapPairs(gr.f1cdx,gr.exome)
gr.f1cdx<-pintersect(tmp)



vcf <- readVcf("TCGA-AX-A0J1.vcf", "hg19")
gr.panel.mut<-rowRanges(vcf)
seqlevelsStyle(gr.panel.mut) <- "UCSC"

# upload mutations and calculate mutations in the given panel
mut.panel.bed<-read.table("www/BRCA_TCGA-AR-A256_msk.bed")
gr.panel.mut<-GRanges(seqnames = Rle(mut.panel.bed$V1),ranges = IRanges(mut.panel.bed$V2,mut.panel.bed$V2),sample=mut.panel.bed$V6)
obs.panel<-sum(countOverlaps(gr.panel.mut, gr.f1cdx))



# calculate tcga muations in given panel (all mutations)
mut.exome.bed<-read.table("data/exome_all_bed/BRCA.bed")
gr.exome.mut<-GRanges(seqnames = Rle(mut.exome.bed$V1),ranges = IRanges(mut.exome.bed$V2,mut.exome.bed$V2),sample=mut.exome.bed$V6)
tmp<-subsetByOverlaps(gr.exome.mut,gr.f1cdx)
panel.mut.cout<-as.data.frame(table(mcols(tmp)$sample))
rownames(panel.mut.cout)<-panel.mut.cout$Var1

# caluculate tcga muations in all coding region (non-syn mutations)
exome.mut.cout<-as.data.frame(table(mut.exome.bed$V6))
rownames(exome.mut.cout)<-exome.mut.cout$Var1

# prepare correlation data
sname<-intersect(rownames(panel.mut.cout),rownames(exome.mut.cout))
final.data<-data.frame(panel.mut.cout[sname,]$Freq,exome.mut.cout[sname,]$Freq)
colnames(final.data)<-c("panel","wes")
rownames(final.data)<-sname

#---model estimate----
region.exome=36.747178
tmp<-data.frame(x=gr.f1cdx)
region.panel=sum((tmp[,3])-(tmp[,2]))/1000000
#obs.panel<-nrow(mut.panel.bed)/region.panel

