region.exome=36.747178
ttype<-scan("~/Projects/tmb/shiny/package/TMBpredict/data/ttype_33",what="s")
exome.final.bed<-read.table("~/Projects/tmb/shiny/package/TMBpredict/data/region_bed/exome.final.bed")
n<-length(ttype)

#=========== get exome tmb==================
wes.tmb.list<-list()
pwd<-"/home/fanghu/Projects/tmb/shiny/package/TMBpredict/data/exome/"
for(i in 1:n){
  tmp<-paste0(pwd,ttype[i],".exome.tmb.txt")
  tmp1<-read.table(tmp,row.names = 1)
  wes.tmb.list[[ttype[i]]]=tmp1[,1]/region.exome
}

#=========== get exome all mut bed==================
all.mut.list<-list()
pwd<-"/home/fanghu/Projects/tmb/shiny/package/TMBpredict/data/exome_all_bed/"
for(i in 1:n){
  tmp<-paste0(pwd,ttype[i],".bed")
  tmp1<-read.table(tmp,header = F)
  all.mut.list[[ttype[i]]]=tmp1
}

#=========== get exome ns mut bed==================
all.ns.list<-list()
pwd<-"/home/fanghu/Projects/tmb/shiny/package/TMBpredict/data/exome_ns_bed/"
for(i in 1:n){
  tmp<-paste0(pwd,ttype[i],".bed")
  tmp1<-read.table(tmp,header = F)
  all.ns.list[[ttype[i]]]=tmp1
}

save(exome.final.bed,all.mut.list,wes.tmb.list,all.ns.list,file = "sysdata1.rda")

