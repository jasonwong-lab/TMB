#-------------------------------------------------------------------------
#                      model development
#-------------------------------------------------------------------------
# mutation type specific Poisson model
train_poisson_smb <- function(counts.exome, exposures.exome, rates.panel) {
  mut.types <- names(rates.panel);
  names(mut.types) <- mut.types;
  
  fits <- lapply(
    mut.types,
    function(mut.type) {
      d <- data.frame(
        count = counts.exome[[mut.type]],
        log_exposure = log(exposures.exome[[mut.type]]),
        log_x = log(rates.panel[[mut.type]])
      );
      glm(count ~ offset(log_exposure) + log_x, family="quasipoisson", data=d)
    }
  );
  
  params <- list(
    exposure = exposures.exome
  );
  
  structure(list(models=fits, params=params), class="tmb_calib_poisson_smb")
}

# mutation type specific Poisson model
train_poisson_smb_2l <- function(counts.exome, exposures.exome, rates.panel) {
  fits <- train_poisson_smb(counts.exome, exposures.exome, rates.panel)$models;
  
  fitteds <- lapply(fits, fitted);
  tmc.exome <- rowSums(counts.exome);

  ensemble <- glm(
    tmc.exome ~ c_a + c_g + c_t + t_a + t_c + t_g + indel,
    data = fitteds
    #		family="gaussian"
  );
  
  params <- list(
    exposure = exposures.exome
  );
  
  structure(list(models.l1=fits, model=ensemble, params=params),
            class="tmb_calib_poisson_smb_2l")
}

# @param smb.panel  specific mutation burden measured by a targeted panel
predict_poisson_smb_2l <- function(model, smb.panel) {
  mut.types <- names(smb.panel);
  names(mut.types) <- mut.types;
  
  preds <- lapply(
    mut.types,
    function(mut.type) {
      d <- data.frame(
        log_exposure = log(model$params$exposure[[mut.type]]),
        log_x = log(smb.panel[[mut.type]])
      );
      # NB response is the count
      predict(model$models.l1[[mut.type]], d, type="response")
    }
  );
  exposure <- model$params$exposure$indel;
  pmax(0, predict(model$model, newdata=preds, type="response") / exposure)
}

#----------------------------------------------------------------------------
#                             PREPARE DATA 
#----------------------------------------------------------------------------
#STEP 1
counts.exome<- read.table("~/Projects/tmb/JNCC/David/tmb-calib-main_modify/data/mut_count/Exome_nonsilent_mut.txt", sep="\t", row.names=1, header=FALSE);
counts.panel<- read.table("~/Projects/tmb/JNCC/data/region_mut/MSK-IMPACT410/MSK-IMPACT410_all_mut.txt", sep="\t", row.names=1, header=FALSE);
sap.name<-intersect(rownames(counts.exome),rownames(counts.panel))
counts.exome<-counts.exome[sap.name,]
colnames(counts.exome) <- c("c_a", "c_g", "c_t", "t_a", "t_c", "t_g", "indel");
counts.panel<-counts.panel[sap.name,]
colnames(counts.panel) <- c("c_a", "c_g", "c_t", "t_a", "t_c", "t_g", "indel");

#STEP 2
opps.exome <- read.table("~/Projects/tmb/JNCC/David/tmb-calib-main_modify/data/nuc_content/Exome.final.bed.nuc", sep="\t", header=FALSE);
colnames(opps.exome) <- c("a", "c", "g", "t")
opps.exome <- as.matrix(opps.exome)[1,];
exposures.exome <- list(
  c_a = unname(opps.exome["c"] + opps.exome["g"]),
  c_g = unname(opps.exome["c"] + opps.exome["g"]),
  c_t = unname(opps.exome["c"] + opps.exome["g"]),
  t_a = unname(opps.exome["t"] + opps.exome["a"]),
  t_c = unname(opps.exome["t"] + opps.exome["a"]),
  t_g = unname(opps.exome["t"] + opps.exome["a"]),
  indel = sum(opps.exome)
);

#STEP 3
opps.panel <- read.table("~/Projects/tmb/JNCC/data/region_bed/MSK-IMPACT410.bed.nuc", sep="\t", header=FALSE);
colnames(opps.panel) <- c("a", "c", "g", "t")
opps.panel <- as.matrix(opps.panel)[1,];
exposures.panel <- list(
  c_a = unname(opps.panel["c"] + opps.panel["g"]),
  c_g = unname(opps.panel["c"] + opps.panel["g"]),
  c_t = unname(opps.panel["c"] + opps.panel["g"]),
  t_a = unname(opps.panel["t"] + opps.panel["a"]),
  t_c = unname(opps.panel["t"] + opps.panel["a"]),
  t_g = unname(opps.panel["t"] + opps.panel["a"]),
  indel = sum(opps.panel)
);

#STEP 4
rates.panel<-mapply(
  function(x, n) {
    (x + 0.5) / (n + 1)
  },
  counts.panel,
  exposures.panel,
  SIMPLIFY=FALSE
)

#-----------------------------------------------------------------------
#                                     train model
#-----------------------------------------------------------------------
possion.model<-train_poisson_smb_2l(counts.exome, exposures.exome, rates.panel)


#-----------------------------------------------------------------------
#                                     PREDICTION
#-----------------------------------------------------------------------
#   INPUT DATA
library(GenomicRanges)
vcf <- readVcf("COAD_test_sap.vcf", "hg19")
gr.panel.mut<-rowRanges(vcf)
seqlevelsStyle(gr.panel.mut) <- "UCSC"
test<-data.frame(ref=gr.panel.mut$REF,alt=unlist(gr.panel.mut$ALT))
mut.count<-as.data.frame(matrix(rep(0,7),nrow=1))
colnames(mut.count) <- c("c_a", "c_g", "c_t", "t_a", "t_c", "t_g", "indel");
for(i in 1:nrow(test)){
  if ((test[i,1]=="C" & test[i,2]=="A") |(test[i,1]=="G" & test[i,2]=="T")){mut.count[1,1]=mut.count[1,1]+1}
  if ((test[i,1]=="C" & test[i,2]=="G") |(test[i,1]=="G" & test[i,2]=="C")){mut.count[1,2]=mut.count[1,2]+1}
  if ((test[i,1]=="C" & test[i,2]=="T") |(test[i,1]=="G" & test[i,2]=="A")){mut.count[1,3]=mut.count[1,3]+1}
  if ((test[i,1]=="T" & test[i,2]=="A") |(test[i,1]=="A" & test[i,2]=="T")){mut.count[1,4]=mut.count[1,4]+1}
  if ((test[i,1]=="T" & test[i,2]=="C") |(test[i,1]=="A" & test[i,2]=="G")){mut.count[1,5]=mut.count[1,5]+1}
  if ((test[i,1]=="T" & test[i,2]=="G") |(test[i,1]=="A" & test[i,2]=="C")){mut.count[1,6]=mut.count[1,6]+1}
  if (nchar(test[i,1])!=nchar(test[i,2])){mut.count[1,7]=mut.count[1,7]+1}
}

rates.panel<-mapply(
  function(x, n) {
    (x + 0.5) / (n + 1)
  },
  mut.count,
  exposures.panel,
  SIMPLIFY=FALSE
)

predict.value<-predict_poisson_smb_2l(possion.model, rates.panel)
pre.tmb<-data.frame(tmb=predict.value*1000000)
rownames(pre.tmb)<-rownames(mut.count)
write.table(pre.tmb,"Predict_TMB.txt",sep = "\t",quote = F)