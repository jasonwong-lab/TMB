region.exome=36.747178
region.msk=1.338520
region.f1cdx=0.830852

ttype<-scan("data/ttype_33",what="s")
n<-length(ttype)
#============get tumor type choices========
tty_list<-list()
for(i in 1:n){
  tty_list[ttype[i]]=ttype[i]
}

#=========== get exome tmb==================
exome.tmb.list<-list()
pwd<-"data/exome/"
for(i in 1:n){
  tmp<-paste0(pwd,ttype[i],".exome.tmb.txt")
  tmp1<-read.table(tmp,row.names = 1)
  exome.tmb.list[[ttype[i]]]=tmp1[,1]/region.exome
}

#=========== get f1cdx tmb==================
f1cdx.tmb.list<-list()
pwd<-"data/f1cdx/"
for(i in 1:n){
  tmp<-paste0(pwd,ttype[i],".f1cdx.tmb.txt")
  tmp1<-read.table(tmp,row.names = 1)
  f1cdx.tmb.list[[ttype[i]]]=tmp1[,1]/region.f1cdx
}

#=========== get msk tmb==================
msk.tmb.list<-list()
pwd<-"data/msk/"
for(i in 1:n){
  tmp<-paste0(pwd,ttype[i],".msk.tmb.txt")
  tmp1<-read.table(tmp,row.names = 1)
  msk.tmb.list[[ttype[i]]]=tmp1[,1]/region.msk
}
