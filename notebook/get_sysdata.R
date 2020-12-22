ttype<-scan("data/ttype_33",what="s")
n<-length(ttype)

# all mut bed file
all.mut.list<-list()
for(i in 1:n){
  tmp<-read.table(paste0("data/exome_all_bed/",ttype[i],".bed"))
  all.mut.list[[ttype[i]]]=tmp
}

# WES TMB file
wes.tmb.list<-list()
for(i in 1:n){
  tmp<-read.table(paste0("data/exome/",ttype[i],".exome.tmb.txt"),row.names = 1)
  wes.tmb.list[[ttype[i]]]=tmp
}

# WES bed file
exome.final.bed <- read.table("data/region_bed/exome.final.bed")

usethis::use_data(all.mut.list,wes.tmb.list, exome.final.bed,internal = TRUE,overwrite =T)