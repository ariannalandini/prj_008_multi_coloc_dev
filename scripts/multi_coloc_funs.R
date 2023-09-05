###library(qgraph)
###library(susieR)
#suppressMessages(library(corrplot))
#suppressMessages(library(coloc))
#suppressMessages(library(data.table))
#suppressMessages(library(bigsnpr))
#suppressMessages(library(ggplot2))
#suppressMessages(library(easyGgplot2))
#suppressMessages(library(igraph))
#suppressMessages(library(RColorBrewer))
#suppressMessages(library(ggnet))
#suppressMessages(library(patchwork))
#suppressMessages(library(stringi))
#suppressMessages(library(plyr))
#suppressMessages(library(dplyr))


### dataset.munge ###
dataset.munge=function(sumstats.file
                        ,map=mappa.loc
                        ,snp.lab="SNP"
                        ,chr.lab="CHR"
                        ,pos.lab="BP"
                        ,a1.lab="A1"
                        ,a0.lab="A2"
                        ,beta.lab="BETA"
                        ,se.lab="SE"
                        ,pval.lab="P"
                        ,freq.lab="FRQ"
                        ,n.lab="N"
                        ,type=NULL
                        ,sdY=NULL
                        ,s=NULL
){

# Load sumstat
  if(is.character(sumstats.file)){
    dataset=fread(sumstats.file)
  }else{
    dataset=sumstats.file
  }

  if(!is.null(a1.lab) & a1.lab%in%names(dataset) & !is.null(a0.lab) & a0.lab%in%names(dataset) ){
    names(dataset)[match(c(a1.lab,a0.lab),names(dataset))]=c("A1","A2")
  }else{
    stop("a0.lab or a1.lab have not been defined or the column is missing")
  }
  if(!is.null(beta.lab)& beta.lab%in%names(dataset)){
    names(dataset)[names(dataset)==beta.lab]="BETA"
  }else{
    stop("beta.lab has not been defined or the column is missing")
  }
  
  if(!is.null(snp.lab) & snp.lab%in%names(dataset)){
    names(dataset)[names(dataset)==snp.lab]="SNP"
  }else{
    stop("snp.lab has not been defined or the column is missing")
  }

### Add also parallel possibility to already have type as column of the GWAS sum stats (instead of specified in the input table)   
  if(exists("type")){
    if(type=="cc" | type=="quant"){
      dataset <- dataset %>% mutate(type=type)
    } else {
    stop("Type has not been defined or has been defined incorrectly - it has to be either 'cc' or 'quant'")
    }
  }
 
#### Map/plink files have unnamed SNP as CHROM:GENPOS_A1_A0, while the GWAS summary statistics as CHROM:GENPOS only. Add also alleles info to avoid losing too many SNPs in merging

##### TEMPORARY FIX FOR SARA - NEED TO DOUBLE CHECK THIS STEP
  dataset$SNP <- gsub("(.*)_\\w+_\\w+$", "\\1", dataset$SNP)
#####
  
  if(!is.null(chr.lab) & chr.lab%in%names(dataset)){
    names(dataset)[names(dataset)==chr.lab]="CHR"
  }else{
    dataset$CHR=map$CHR[match(dataset$SNP,map$SNP)]
  }
  
  if(!is.null(pos.lab) & pos.lab%in%names(dataset)){
    names(dataset)[names(dataset)==pos.lab]="BP"
  }else{
    dataset$BP=map$BP[match(dataset$SNP,map$SNP)]
  }
  
  if(!is.null(freq.lab) & freq.lab%in%names(dataset)){
    names(dataset)[names(dataset)==freq.lab]="FRQ"
  }else{
    dataset$FRQ=map$MAF[match(dataset$SNP,map$SNP)]
  }
  
  if("FRQ" %in% colnames(dataset)){
    dataset$MAF=dataset$FRQ
    dataset$MAF[dataset$MAF>0.5]=1-dataset$MAF[dataset$MAF>0.5]
  }
  
  if(!is.null(n.lab) & n.lab%in%names(dataset)){
    names(dataset)[names(dataset)==n.lab]="N"
  }else{
    N_hat<-median(1/((2*dataset$MAF*(1-dataset$MAF))*dataset$SE^2),na.rm = T)
    dataset$N=ceiling(N_hat)
  }
  
  if(!is.null(pval.lab) & pval.lab%in%names(dataset)){
    names(dataset)[names(dataset)==pval.lab]="P"
### Check if p-value column provided is log10 transformed. If yes, compute original p-value
    if (!all(dataset$P >= 0 & dataset$P <= 1)) {
      dataset <- dataset %>% mutate(P=10^(-P))
    }
  }else{
    dataset$P=pchisq((dataset$BETA/dataset$SE)^2,df=1,lower=F)
  }
# Add variance of beta  
  dataset$varbeta=dataset$SE^2
  
# Add sdY/s
  if(type=="cc" && !(is.null(s))){
    dataset$s=s
  } else if(type=="cc" && is.null(s)){
    #### Is this correct?? Is "s" strictly necessary for cc traits??
    stop("Please provide s, the proportion of samples who are cases")
  }

  if(type=="quant" && !(is.null(sdY))){
    dataset$sdY <- sdY
  } else if(type=="quant" && is.null(sdY)){
    dataset$sdY <- sdY.est(dataset$varbeta, dataset$MAF, dataset$N)
  }
  

# Match with locus reference map
  dataset <- dataset[which(dataset$SNP %in% map$SNP),]
  dataset <- dataset[match(map$SNP,dataset$SNP),]
  
  flip=dataset[,c("SNP","CHR","BP","A2","A1","BETA")]
  names(flip)=c("rsid","chr","pos","a0","a1","beta")
  names(map)=c("rsid","chr","pos","maf","a1","a0")

  flip.t=snp_match(sumstats=flip,
                   info_snp=map,
                   join_by_pos=FALSE,
                   strand_flip=FALSE,
                   match.min.prop=0)
    
  dataset=dataset[match(flip.t$rsid,dataset$SNP),]
  dataset$A1=flip.t$a1
  dataset$A2=flip.t$a0
  dataset$BETA=flip.t$beta
  
  dataset <- dataset %>%
    select("SNP","CHR","BP","A1","A2","BETA","varbeta","P","MAF","N", "type", any_of(c("s", "sdY")))
  
  if(type=="cc"){
    names(dataset)=c("snp","chr","pos","a1","a0","beta","varbeta","pvalues","MAF","N","type","s")
  } else if(type=="quant"){
    names(dataset)=c("snp","chr","pos","a1","a0","beta","varbeta","pvalues","MAF","N","type","sdY")
  }
  dataset
}


### cojo.ht ###
### Performs --cojo-slct first to identify all indipendent SNPs and --cojo-cond then to condition upon identified SNPs
cojo.ht=function(D=datasets[[1]]
                 ,plink.bin="/project/alfredo/software/plink/1.90_20210606/plink"
                 ,gcta.bin="/project/alfredo/software/GCTA/1.94.0beta/gcta64"
                 ,bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/p01_output/ukbb_all_30000_random_unrelated_white_british"
                 ,p.tresh=1e-4
                 ,maf.thresh=0.0001){
  
  random.number=stri_rand_strings(n=1, length=20, pattern = "[A-Za-z0-9]")
  
  write(D$snp,ncol=1,file=paste0(random.number,".snp.list"))
  system(paste0(plink.bin," --bfile ",bfile," --extract ",random.number,".snp.list --maf ", maf.thresh, " --make-bed --freqx --out ",random.number))
  
  freqs=fread(paste0(random.number,".frqx"))
  freqs$FreqA1=(freqs$'C(HOM A1)'*2+freqs$'C(HET)')/(2*(rowSums(freqs[,c("C(HOM A1)", "C(HET)", "C(HOM A2)")])))
  D$FREQ=freqs$FreqA1[match(D$snp,freqs$SNP)]
  idx=which(D$a1!=freqs$A1)
  D$FREQ[idx]=1-D$FREQ[idx]
  D$se=sqrt(D$varbeta)

  D <- D %>%
    select("snp","a1","a0","FREQ","beta","se","pvalues","N","type", any_of(c("sdY","s"))) %>%
    rename("SNP"="snp","A1"="a1","A2"="a0","freq"="FREQ","b"="beta","p"="pvalues")
  
  write.table(D,file=paste0(random.number,"_sum.txt"),row.names=F,quote=F,sep="\t")

# step1 determine independent snps
  system(paste0(gcta.bin," --bfile ", random.number, " --cojo-p ", p.tresh, " --maf ", maf.thresh, " --extract ", random.number, ".snp.list --cojo-file ", random.number, "_sum.txt --cojo-slct --out ", random.number, "_step1"))
  
  ind.snp=fread(paste0(random.number,"_step1.jma.cojo"))
  dataset.list=list()
  dataset.list$ind.snps=ind.snp
  dataset.list$results=list()
  
  
  if(nrow(ind.snp)>1){
    for(i in 1:nrow(ind.snp)){
      
      write(ind.snp$SNP[-i],ncol=1,file=paste0(random.number,"_independent.snp"))
      print(ind.snp$SNP[-i])
      
      system(paste0(gcta.bin," --bfile ",random.number, " --maf ", maf.thresh, " --extract ",random.number,".snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
      
#### KEEP WORKING FROM HERE ###### try catch COJO errors - commented out for not breaking the code until this piece in finished
      
#      test_res <- tryCatch(
#        system(paste0(gcta.bin," --bfile ",random.number,"  --extract ",random.number,".snp.list  --cojo-file ",random.number,"_sum.txt  --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"), intern=TRUE),
#        warning = function(war) {
          # Handle the warning here
          # You can print an error message or perform any other desired action
#          print(paste("Hey hey! Error:", conditionMessage(war)))
          # Return a default value or NULL to indicate failure
#          NULL
#        }
#      )
      
############

# Re-add type and sdY/s info
      type <- unique(D$type)
      if(type=="quant"){
        step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
          mutate(type=type, sdY=unique(D$sdY))
      } else if(type=="cc"){
        step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
          mutate(type=type, s=unique(D$s))
      }
      dataset.list$results[[i]]=step2.res
      names(dataset.list$results)[i]=ind.snp$SNP[i]
    }
    
  } else {

### NB: COJO here is performed ONLY for formatting sakes - No need to condition if only one signal is found!!        
    write(ind.snp$SNP[1],ncol=1,file=paste0(random.number,"_independent.snp"))
    system(paste0(gcta.bin," --bfile ",random.number," --cojo-p ",p.tresh, " --maf ", maf.thresh, " --extract ",random.number,".snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
    step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE)

#### Add back top SNP, removed from the data frame with the conditioning step
    step2.res <- rbind.fill(ind.snp,step2.res) %>%
        select("Chr","SNP","bp","refA","freq","b","se","p","n","freq_geno")
    step2.res <- cbind(
      step2.res,
      D %>% select(type, any_of(c("sdY", "s"))) %>% slice_sample(n=nrow(step2.res))
    )
    
    dataset.list$results[[1]]=step2.res
    names(dataset.list$results)[1]=ind.snp$SNP[1]
  }
  system(paste0("rm *",random.number,"*"))
  dataset.list
}


### plot.cojo.ht ###
plot.cojo.ht=function(cojo.ht.obj){
#  library(ggplot2)
#  library(patchwork)
  
  if(nrow(cojo.ht.obj$ind.snps)>1){
    
    whole.dataset=c()
    for(i in 1:nrow(cojo.ht.obj$ind.snps)){
      
      tmp=cojo.ht.obj$results[[i]]
      tmp$signal=cojo.ht.obj$ind.snps$SNP[i]
      whole.dataset=rbind(whole.dataset,tmp)
    }
    
    p1 <- ggplot(cojo.ht.obj$results[[i]], aes(x=bp,y=-log10(p))) +
      geom_point(alpha=0.6,size=3)+
      theme_minimal()+
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=SNP),size=6,shape=23)
    
    p2 <- ggplot(whole.dataset,aes(x=bp,y=-log10(pC),color=signal)) +
      facet_grid(signal~.) +
      geom_point(alpha=0.8,size=3) +
      theme_minimal() +
      ggtitle("Conditioned results")
    
    p3 <- p1/p2 + plot_layout(heights = c(1, nrow(cojo.ht.obj$ind.snps)+0.2))
    
  }else{
    
    p3 <- ggplot(cojo.ht.obj$results[[1]], aes(x=bp,y=-log10(p))) +
      geom_point(alpha=0.6,size=3)+
      theme_minimal()+
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=SNP),size=6,shape=23)
  }
  (p3)
}


### coloc.prep.table
coloc.prep.table=function(pairwise.list, conditional.datasets, loci.table.tmp){
  ### Select only 1) genome-wide significant or 2) with conditioned p-value < 1e-6 SNPs
  final.locus.table.tmp <- as.data.frame(rbindlist(lapply(names(conditional.datasets), function(x){
    tmp <- conditional.datasets[[x]]$ind.snps
    if(nrow(tmp)>1){tmp <- tmp %>% filter(p<5e-8 | pJ<1e-6)}
    tmp <- tmp %>%
      mutate(
        start=unique(loci.table.tmp$start),
        end=unique(loci.table.tmp$end),
        trait=x,
        othA=NA,
        pan.locus=locus,
        sub_locus=NA
      )
    ### Add other allele from map    
    for(n in 1:nrow(tmp)){
      alleles=unlist(mappa.loc[mappa.loc$SNP==tmp$SNP[n],c("A1","A2")])
      tmp$othA[n]=alleles[!(alleles%in%tmp$refA[n])]
    }
    tmp
  })))
  return(final.locus.table.tmp)
}


### colo.cojo.ht ###
colo.cojo.ht=function(conditional.dataset1=conditional.datasets[[pairwise.list[i,1]]]
                      ,conditional.dataset2=conditional.datasets[[pairwise.list[i,2]]]
                      ,p.threshold.cond=1e-6
                      ,p.threshold.orig=5e-8){
  
  if(length(grep("pJ",names(conditional.dataset1$ind.snps)))>0){
    
    hits.t1=conditional.dataset1$ind.snps$SNP[conditional.dataset1$ind.snps$pJ<p.threshold.cond | conditional.dataset1$ind.snps$p <p.threshold.orig]
  }else{
    hits.t1=conditional.dataset1$ind.snps$SNP
  }
  
  if(length(grep("pJ",names(conditional.dataset2$ind.snps)))>0){
    
    hits.t2=conditional.dataset2$ind.snps$SNP[conditional.dataset2$ind.snps$pJ<p.threshold.cond | conditional.dataset2$ind.snps$p<p.threshold.orig]
  }else{
    hits.t2=conditional.dataset2$ind.snps$SNP
  }
  
  if(length(hits.t2)>0 & length(hits.t1)>0){
    
    coloc.final=list()
    for(i in hits.t1){
      for(j in hits.t2){
        
        D1=conditional.dataset1$results[[i]]     
        D2=conditional.dataset2$results[[j]] 
        
        if(length(grep("bC", names(D1)))>0){
          D1 <- D1 %>%
            select("SNP","Chr","bp","bC","bC_se","n","pC","freq","type",any_of(c("sdY","s"))) %>%
            rename("snp"="SNP","chr"="Chr","position"="bp","beta"="bC","varbeta"="bC_se","N"="n","pvalues"="pC","MAF"="freq")
        }else{
          D1 <- D1 %>%
            select("SNP","Chr","bp","b","se","n","p","freq","type",any_of(c("sdY","s"))) %>%
            rename("snp"="SNP","chr"="Chr","position"="bp","beta"="b","varbeta"="se","N"="n","pvalues"="p","MAF"="freq")
        }
        D1$varbeta=D1$varbeta^2
        D1=na.omit(D1)
        
        if(length(grep("bC", names(D2)))>0){
          D2 <- D2 %>%
            select("SNP","Chr","bp","bC","bC_se","n","pC","freq","type",any_of(c("sdY","s"))) %>%
            rename("snp"="SNP","chr"="Chr","position"="bp","beta"="bC","varbeta"="bC_se","N"="n","pvalues"="pC","MAF"="freq")
        }else{
          D2 <- D2 %>%
            select("SNP","Chr","bp","b","se","n","p","freq","type",any_of(c("sdY","s"))) %>%
            rename("snp"="SNP","chr"="Chr","position"="bp","beta"="b","varbeta"="se","N"="n","pvalues"="p","MAF"="freq")
        }
        D2$varbeta=D2$varbeta^2
        D2=na.omit(D2)
        
        colo.res <- suppressWarnings(coloc.abf(D1,D2)) ### BETTER NOT TO SUPPRESS WARNINGS?
## Save coloc summary        
        colo.sum=data.frame(t(colo.res$summary))
        colo.sum$hit1=i
        colo.sum$hit2=j
#        coloc.summary=rbind(coloc.summary,colo.sum)

## Save coloc result by SNP
        colo.full_res <- colo.res$results %>% 
          select(snp,position,lABF.df1,lABF.df2,SNP.PP.H4) %>% 
          mutate(hit1=i, hit2=j)
        colo.all <- list(summary=colo.sum, results=colo.full_res)
## Organise all in a list of lists (each list is composed of summary + results)
        coloc.final <- c(coloc.final, list(colo.all))
      }
    }
  }else{
    coloc.final=NULL
  }
  return(coloc.final)
}


### coloc.subgrouping 
coloc.subgrouping <- function(final.colocs.H4, final.locus.table.tmp){
  # Add sublocus to locus table
  for(k in 1:nrow(final.colocs.H4)){
    
    final.locus.table.tmp$sub_locus[
      final.locus.table.tmp$trait==final.colocs.H4$t1[k] &
        final.locus.table.tmp$SNP==final.colocs.H4$hit1[k]]=final.colocs.H4$g1[k]
    
    final.locus.table.tmp$sub_locus[
      final.locus.table.tmp$trait==final.colocs.H4$t2[k] &
        final.locus.table.tmp$SNP==final.colocs.H4$hit2[k]]=final.colocs.H4$g1[k]
  }
  
  # Add sublocus also for non-colocalising loci
  if(any(is.na(final.locus.table.tmp$sub_locus))){
    idx=which(is.na(final.locus.table.tmp$sub_locus))
    pri=max(final.locus.table.tmp$sub_locus,na.rm=T)+1
    final.locus.table.tmp$sub_locus[idx]=pri:(pri+length(idx)-1)
  }
  
  final.locus.table.tmp=as.data.frame(final.locus.table.tmp)[,col.order]
  return(final.locus.table.tmp)
} 


### coloc.plot ###
coloc.plot <- function(x, outpath=NULL){
  
  if(is.null(x)){
    cat("\nError: table summarising all colocalizing results (PP.H4 >= 0.9) is missing\n")
  } else {
    
    per.plot.data=c()
    
    for(i in 1:nrow(x)){
      tmp=conditional.datasets[[x$t1[i]]]$results[[x$hit1[i]]]
      tmp$label=paste(x$t1[i],x$hit1[i],sep="-")
      tmp$group=x$g1[i]
      
      if(length(grep("pC",names(tmp)))>0){
        
        tmp=tmp[,c("bp","pC","label","group")]
        names(tmp)=c("bp","p","label","group")
        
      } else {
        
        tmp=tmp[,c("bp","p","label","group")]
      }  
      
      per.plot.data=rbind(per.plot.data,tmp)
      tmp=conditional.datasets[[x$t2[i]]]$results[[x$hit2[i]]]
      tmp$label=paste(x$t2[i],x$hit2[i],sep="-")
      tmp$group=x$g1[i]
      
      if(length(grep("pC",names(tmp)))>0){
        
        tmp=tmp[,c("bp","pC","label","group")]
        names(tmp)=c("bp","p","label","group")
      } else {
        tmp=tmp[,c("bp","p","label","group")]
      }
      
      per.plot.data=rbind(per.plot.data,tmp)
      per.plot.data=per.plot.data[order(per.plot.data$group),]
      per.plot.data$label=factor(per.plot.data$label,levels=unique(per.plot.data$label))
      
      x$locus=locus
    }
    
### Plot
    pdf(paste0(outpath, "locus_", locus, "_colocalization_plot.pdf"),
        width=14, height=4*length(unique(per.plot.data$label)))  
    p1 <- ggplot(per.plot.data, aes(x=bp,y=-log10(p), fill=as.character(group))) +
      geom_point(shape=21, alpha=0.9, size=4) +
      geom_hline(yintercept = 0, linewidth=0.5) +
      facet_grid(label~.,scales="free") +
      theme_minimal() +
#      scale_fill_manual(values=viridis(length(unique(per.plot.data$group)))) +
#      scale_color_manual(values=viridis(length(unique(per.plot.data$group)))) +
      theme(strip.text.y.right = element_text(angle = 270, size=20),
            axis.line.y = element_line(),
            legend.position = "none") +
      ggtitle(paste("Locus ",locus," conditional regional plot"))
    print(p1)
    dev.off()
  }    
}


### p.minim ###
p.minim=function(x)min(pchisq(as.numeric(unlist((dataset[[i]][x,6:ncol(dataset[[i]])]^2))),df=1,lower=F))
minimiser=function(x){
  p1=unlist(susie.list[[unlist(x["d1"])]]$sets$min.p[paste0("L",unlist(x["idx1"]))])
  p2=unlist(susie.list[[unlist(x["d2"])]]$sets$min.p[paste0("L",unlist(x["idx2"]))])
  min(c(p1,p2))
  
}


### locus.joyplot ###
locus.joyplot=function(x, window.size=10000, susie.res=susie.list){
  start=min(x$pos)
  final.end=max(x$pos)
  log10p=x[,6:ncol(x)]
  res.all=c()
  
  while(start<final.end){
    end=start+window.size
    idx=which(x$pos>=start & x$pos<=end)
    if(nrow(log10p[idx,])>0){
      means=apply(log10p[idx,],2,function(x)x[which.max(abs(x))])
      res.tmp=c(mean(c(start,end)),means)
      res.all=rbind(res.all,res.tmp)
    }
    start=start+window.size/2
  }
  
#  library(reshape2)
  res.all=as.data.frame(res.all)
  ordine=hclust(as.dist(1-cor(as.matrix(res.all[,-1]),use ="pairwise.complete.obs")),method = "ward.D2")$order
  ordine=(names(res.all)[-1]) [ordine]
  names(res.all)[1]="pos"
  
  res.all=melt(res.all,id.vars = "pos")
  
  names(res.all)=c("pos","trait","log10p")
  res.all$trait=factor(res.all$trait,levels=ordine)
  
  ggplot(res.all,aes(x=pos,y=log10p))+geom_area(fill="blue",alpha=0.4)+
    geom_hline(yintercept = -log10(5e-8),linetype="dashed",color="red")+
    geom_hline(yintercept = log10(5e-8),linetype="dashed",color="red")+
    facet_grid(scales="free",rows = "trait")+
    theme_minimal()+theme(panel.grid = element_blank())
}







### bin2lin2 ###
bin2lin2=function (D, dotplot = FALSE){
  if (D$type != "cc") 
    stop("type != cc")
  D1 <- D
  beta = (D1$beta/sqrt(D1$varbeta)) / sqrt (D1$N * 2*(D1$MAF)*(1-D1$MAF)) 
  D1$varbeta <- (beta/(D1$beta/sqrt(D1$varbeta)))^2
  D1$beta <-beta
  
  z1 <- D1$beta/sqrt(D1$varbeta)
  z <- D$beta/sqrt(D$varbeta)
  if(dotplot==TRUE){
    plot(z1,z)
  }
  D1$quality <- 1
  D1
}








pleio.table=function(conditional.datasets=conditional.datasets,loc.table=NA,plot=FALSE,plot.file=NULL, index.trait=NULL){
  
  
  ## Remove duplicates
  duplicati=table(loc.table$trait)
  if(any(duplicati>1)){
    doppi=names(duplicati)[duplicati>1]
    for(i in doppi){
      tmp=conditional.datasets[[i]]$ind.snps
      tmp=tmp[which(tmp$SNP%in%loc.table$SNP[loc.table$trait==i]),]
      min.snp=tmp$SNP[which.min(tmp$pJ)]
      loc.table=loc.table[which(loc.table$trait!=i | (loc.table$trait==i & loc.table$SNP==min.snp)),]
    }
    
  }
  
  
  if(nrow(loc.table)>1){
    for(i in 1:nrow(loc.table)){
      
      if(i==1){
        merged.datasets=conditional.datasets[[loc.table$trait[i]]]$results[[loc.table$SNP[i]]]
        if("bC"%in% names(merged.datasets)){
          
          merged.datasets=merged.datasets[,c("SNP","bC","bC_se","pC")]
          names(merged.datasets)=c("SNP",paste(c("beta","se","pval"),loc.table$trait[i],sep=".."))
        }else{
          
          merged.datasets=merged.datasets[,c("SNP","b","se","p")]
          names(merged.datasets)=c("SNP",paste(c("beta","se","pval"),loc.table$trait[i],sep=".."))
          
          
        }
        
        
      }else{
        
        merged.datasets2=conditional.datasets[[loc.table$trait[i]]]$results[[loc.table$SNP[i]]]
        if("bC"%in% names(merged.datasets2)){
          
          merged.datasets2=merged.datasets2[,c("SNP","bC","bC_se","pC")]
          names(merged.datasets2)=c("SNP",paste(c("beta","se","pval"),loc.table$trait[i],sep=".."))
        }else{
          
          merged.datasets2=merged.datasets2[,c("SNP","b","se","p")]
          names(merged.datasets2)=c("SNP",paste(c("beta","se","pval"),loc.table$trait[i],sep=".."))
          
          
        }
        merged.datasets=merge(merged.datasets,merged.datasets2,by="SNP",all = FALSE,suffixes =c(""))
        
      }
      
    }
    
    merged.datasets=na.omit(merged.datasets)
    if(is.null(index.trait)){
      
      top.snp=as.data.frame(merged.datasets[which.min(apply(merged.datasets[,grep("pval..",names(merged.datasets))],1,min)),])
    }else if (length(grep(index.trait,names(merged.datasets)))==0){
      top.snp=as.data.frame(merged.datasets[which.min(apply(merged.datasets[,grep("pval..",names(merged.datasets))],1,min)),])
    }else{
      top.snp=as.data.frame(merged.datasets[which.min(merged.datasets[,paste0("pval..",index.trait)]),])
    }
    
    top.snp=reshape2::melt(top.snp)
    matric=matrix(unlist(strsplit(as.character(top.snp$variable),split="\\.\\.")),ncol=2,byrow = T)
    
    top.snp$variable=matric[,1]
    top.snp$trait=matric[,2]
    if(plot==TRUE){
     
      mappa=conditional.datasets[[1]]$results[[1]][,c("SNP","bp")]
      mappa=na.omit(mappa)
      mappa=merge(mappa,merged.datasets[,c(1,grep("pval..",names(merged.datasets)))],by="SNP")
      per.plot=reshape2::melt(mappa[,-1],id.vars = "bp")
      per.plot$variable=gsub("pval..","",per.plot$variable)
      pdf(plot.file,width=14,height=7*nrow(loc.table))
      p1=ggplot(per.plot,aes(x=bp,y=-log10(value),fill=variable))+
        geom_point(shape=21,color="black",size=4,alpha=0.9)+
        scale_fill_manual(values=RColorBrewer::brewer.pal(n=nrow(loc.table),name = "Paired"))+
        facet_grid(variable~., scales='free')+
        theme_minimal()+
        theme(axis.line = element_line(colour="black"),panel.spacing = unit(3, "lines"))
      print(p1)
      dev.off()
    }
    
  }else{
    
    merged.datasets=conditional.datasets[[loc.table$trait[1]]]$results[[loc.table$SNP[1]]]
    if("bC"%in% names(merged.datasets)){
      
      merged.datasets=merged.datasets[,c("SNP","bC","bC_se","pC")]
      names(merged.datasets)=c("SNP","beta","se","pval")
    }else{
      
      merged.datasets=merged.datasets[,c("SNP","b","se","p")]
      names(merged.datasets)=c("SNP","beta","se","pval")
    }
    top.snp=merged.datasets[which.min(merged.datasets$p),]
    names(top.snp)=c("SNP",paste(c("beta","se","pval"),loc.table$trait[1],sep=".."))
    top.snp=reshape2::melt(top.snp)
    matric=matrix(unlist(strsplit(as.character(top.snp$variable),split="\\.\\.")),ncol=2,byrow = T)
    top.snp$variable=matric[,1]
    top.snp$trait=matric[,2]
  }
  top.snp
}


## Taken from coloc package https://github.com/chr1swallace/coloc/blob/HEAD/R/claudia.R
##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##' 
##' @title Estimate trait variance, internal function
##' @param vbeta vector of variance of coefficients
##' @param maf vector of MAF (same length as vbeta)
##' @param n sample size
##' @return estimated standard deviation of Y
##' 
##' @author Chris Wallace
sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}


#### package.loader - Load packages if available, install them first if not
# Check if the package is already installed
package.loader <- function(package_name){
  if(!require(package_name, character.only = TRUE)) {
    # If not installed, install the package
    install.packages(package_name)
    # Load the package
    library(package_name, character.only = TRUE)
  } else {
    # If already installed, just load the package
    library(package_name, character.only = TRUE)
  }
}


#### final.plot - Final summary plot function
final.plot <- function(locus,
                       final.locus.table.tmp,
                       conditional.datasets,
                       by_snp_PPH3=NULL, ### not necessarily present!
                       inter=NULL, ### not necessarily present!
                       output=opt$output
){
  
  ### Among NOT colocalising traits, extract SNP having the highest lABF (scaled)
    if(length(by_snp_PPH3)>0){
      by_snp_PPH3_final <- rbindlist(lapply(by_snp_PPH3, function(x){
      x %>%
        dplyr::select(-position, -SNP.PP.H4) %>%
        gather("trait", "lABF", -snp,-hit1,-hit2,-t1,-t2,-pan.locus) %>%
        mutate(trait=ifelse(trait=="lABF.df1", unique(t1), unique(t2))) %>%
        mutate(cojo_snp=ifelse(trait==t1, unique(hit1), unique(hit2))) %>%
        dplyr::select(-hit1,-hit2,-t1,-t2)
    })) %>% group_split(trait, cojo_snp)
    
    by_snp_PPH3_final <- as.data.frame(rbindlist(lapply(by_snp_PPH3_final, function(x){
      x %>% distinct(snp, .keep_all = T) %>%
        mutate(bf=exp(lABF)) %>% 
        arrange(desc(bf)) %>% 
        mutate(joint.pp.cv = bf/sum(bf)) %>%
        dplyr::filter(joint.pp.cv==max(joint.pp.cv))
    }))) %>% 
      dplyr::select(pan.locus, trait, cojo_snp, snp, joint.pp.cv) %>%
      dplyr::rename(causal_snp=snp)
    }
  
  
  #### Final locus table - flag colocalising and not colocalising groups
  loci_table <- final.locus.table.tmp %>%
    dplyr::select(pan.locus, sub_locus, trait, SNP) %>%
    group_by(sub_locus) %>%
    mutate(coloc_out=ifelse(n()>1, "H4", "H3")) %>%
    dplyr::rename(cojo_snp=SNP)
  
  ### H4 traits - add SNPs with highest joint.pp.cv among intersection cs
  #loci_table %>% 
  #  filter(coloc_out=="H4") %>%
  #  left_join(inter, by=c("pan.locus", "sub_locus"="g1"))
  
  ### H4 - extract cs intersection SNP having highest joint lABF
  #inter <- inter %>% 
  #  group_by(g1) %>%
  #  filter(joint.pp.cv==max(joint.pp.cv)) %>%
  #  rename(causal_snp=snp, pp.cv=joint.pp.cv) %>%
  #  left_join(inter_info %>% select(-bf,-pp.cv,-cred.set,-joint.pp.cv,-lABF), by=c("causal_snp"="snp","g1","pan.locus"), multiple = "all")
  
  
  ### Extract beta info from conditional datasets  
  data_sub <- unlist(conditional.datasets, recursive=F)
  
  ## Retrieve only conditioned results
  data_sub <- data_sub[grep("results", names(data_sub))]
  
  # Another round of unlisting
  data_sub <- unlist(data_sub, recursive=F)
  
  # Add trait and cojo hit as dataframe columns
  for(i in 1:length(data_sub)){
    data_sub[[i]] <- data_sub[[i]] %>% 
      mutate(trait=gsub("(\\w+).results.(.*$)", "\\1", names(data_sub)[[i]]),
             cojo_snp=gsub("(\\w+).results.(.*$)", "\\2", names(data_sub)[[i]])
      )
  }
  data_sub <- as.data.frame(rbindlist(data_sub, fill=TRUE))
  
  
  #### ASSEMBLE FINAL TABLE FOR PLOTTING
  
  if(any(loci_table$coloc_out=="H4")){
    # H4  
    temp_H4 <- loci_table %>%
      dplyr::filter(coloc_out=="H4") %>%
      left_join(inter %>% dplyr::rename(causal_snp=snp),
                by=c("pan.locus", "sub_locus"="g1"), multiple="all")
  }
  
  if(any(loci_table$coloc_out=="H3")){
    # H3  
    temp_H3 <- loci_table %>% 
      dplyr::filter(coloc_out=="H3") %>%
      left_join(by_snp_PPH3_final, by=c("pan.locus","trait","cojo_snp"))
  }
    
  if(exists("temp_H3") & exists("temp_H4")){ final <- rbind(temp_H3,temp_H4) }  
  if(exists("temp_H3") & !exists("temp_H4")){ final <- temp_H3 }  
  if(!exists("temp_H3") & exists("temp_H4")){ final <- temp_H4 }  
  
  final <- as.data.frame(
    final %>% left_join(data_sub, by=c("causal_snp"="SNP", "trait", "cojo_snp")) %>%
      dplyr::select(pan.locus,sub_locus,trait,causal_snp,Chr,bp,freq,any_of(c("b","bC")),joint.pp.cv))

## Adjustment for plotting
    if("bC" %in% names(final)){
      final <- final %>% mutate(beta=ifelse(is.na(bC), b, bC))
    } else {
      final <- final %>% mutate(beta=b)
    }
        
  final <- final %>%
        mutate(dir=ifelse(beta<0, "-", "+")) %>%
        mutate(group=paste0(sub_locus, " ", causal_snp)) %>%
        arrange(bp)

  #### To delete - just for script developing sake
  #  fwrite(final,
  #         paste0(opt$output, "/results/locus_", locus, "_table_for_final_plot.tsv"),
  #         sep="\t", quote=F, na=NA)
  

#### PLOT #### 
  
# Set plot boundiaries  
  bp_max = max(final$bp) + 250000
  bp_min = min(final$bp) - 250000
  chr = unique(final$Chr)

### Plot 1 - causal SNPs beta
  
# lock in factor level order
  final$group <- factor(final$group, levels = unique(final$group))
# Set x axis labels  
  x_axis <- gsub("\\d+ (.*)", "\\1", unique(final$group))
# Set spacing between each data point
  spacing <- (bp_max-bp_min)/(length(x_axis)+1)
# Calculate position of each data point based on spacing just calculated - so it's easier to match in the next plot!!
  final <- as.data.frame(
    final %>% 
    group_by(group) %>%
    mutate(group_index = cur_group_id()) %>%
    mutate(breaks=bp_min+(spacing*group_index))
  )

  p1 <- ggplot(final) + 
    geom_point(aes(x=breaks, y=trait, size=beta, color=dir)) +
    scale_x_continuous(
      limits=c(bp_min,bp_max),
      expand=c(0,0),
      breaks=unique(final$breaks),
      labels=x_axis) +
    ylab("Magnitude and\ndirection of effect\n") +
    scale_color_manual(values = c("-"="#10b090", "+"="#dc143c")) +
    guides(size = "none") +
    guides(color=guide_legend(title="Direction of effect   ", override.aes=list(size=4)
    )) +
    theme_bw() +
    theme(
      legend.background=element_rect(linewidth=0.5, linetype="solid", colour="black"),
      legend.text = element_text(size = 12, vjust=0.5),
      legend.title = element_text(size = 10, vjust=0.5),
      legend.position = "bottom",
      legend.direction="horizontal",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=9),#, angle=45, hjust=1),
      axis.text.y = element_text(size=9),
      panel.border = element_rect(color = "black", fill=NA, linewidth=1),
      panel.grid.minor = element_blank()
    )
  
  

### Plot 2 - Likely causal SNPs labels and lines connecting the equally spaced betas to actual position on chromosome
  
  p2 <- ggplot(final %>% distinct(group, .keep_all=T), aes(x = bp)) +
    geom_segment(aes(y=0.5, yend=1, x=breaks, xend=breaks), color="black") +
    geom_segment(aes(y=0, yend=0.5, x=bp, xend=breaks), color="black") +
#    geom_text(aes(y=1.02, x=breaks, label=x_axis), vjust = 0) +
    scale_x_continuous(
      limits=c(bp_min,bp_max),
      expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    )
  

### Plot 3 - gene annotation
  
  # Retrieve info on the genes
  ref_genes <- genes(EnsDb.Hsapiens.v75)
  ref_genes <- ref_genes[ref_genes@elementMetadata$gene_biotype=='protein_coding',]
  ind <- findOverlaps(GRanges(seqnames=chr,IRanges(start=bp_min,end=bp_max)),ref_genes,type="any")
  a <- ref_genes[ind@to,]
  
  granges2df <- function(x) {
    df <- as(x, "data.frame")
    df <- df[,c("seqnames","start","end","strand","group_name",'exon_id')]
    colnames(df)[1] <- "chromosome"
    colnames(df)[5] <- "transcript"
    df
  }
  
  txdf <- ensembldb::select(EnsDb.Hsapiens.v75,
                            keys=keys(EnsDb.Hsapiens.v75, "GENEID"),
                            columns=c("GENEID","TXID", 'SYMBOL'),
                            keytype="GENEID")
  ebt <- exonsBy(EnsDb.Hsapiens.v75, by="tx")
  
  ## Arrange info about all transcripts
  d <- list()
  
  for(i in 1:length(a@elementMetadata$gene_id)){
    idx <- txdf$GENEID ==  a@elementMetadata$gene_id[i]
    txs <- txdf$TXID[idx]
    #all the xons for these transcripts
    ebt2 <- ebt[txs]
    df <- granges2df(ebt2)
    df$gene <- a@elementMetadata$gene_id[i]
    df$symbol <-  txdf[match(df$gene, txdf$GENEID), ]$SYMBOL 
    d[[i]] <- df
  }
  
  d <- do.call(rbind, d) %>%
    mutate(y=as.numeric(ifelse(strand=="+", 1, -1)))
  
  gene.starts <- by(d$start,d$symbo,FUN = min)
  gene.end <- by(d$end,d$symbo,FUN = max)
  gene.strand <- by(d$y,d$symbo,FUN = unique)
  
  g.table=as.data.frame(cbind(gene.starts,gene.end,gene.strand))
  g.table=g.table[order(g.table$gene.starts),]
#  g.table$gene.starts <- round(g.table$gene.starts/1e+6,2)
#  g.table$gene.end <- round(g.table$gene.end/1e+6,2)
#  g.table$gene.strand <- ifelse(g.table$gene.strand=="1", "+", "-")
  
  g.table <- g.table %>%
    mutate(gene.starts=round(gene.starts/1e+6,2), gene.end=round(gene.end/1e+6,2)) %>%
    mutate(gene.starts=ifelse(gene.starts<round(bp_min/1e+6,2),
      round(bp_min/1e+6,2),gene.starts)) %>%
    mutate(gene.end=ifelse(gene.end>round(bp_max/1e+6,2),
                              round(bp_max/1e+6,2),gene.end)) 
  
  g.table$gene.strand <- as.character(g.table$gene.strand)
  lab.table=rowMeans(g.table[,1:2])
  
  if(length(names(lab.table))>1){
  names(lab.table)[c(TRUE, FALSE)] <- paste0("\n\n\n", names(lab.table)[c(TRUE, FALSE)])
  names(lab.table)[c(FALSE,TRUE)] <- paste0("\n\n\n\n\n", names(lab.table)[c(FALSE,TRUE)])
  } else {
    names(lab.table) <- paste0("\n\n\n", names(lab.table))
  }
  
  
  ### Gene position - doesn't make a lot of sense with the two strands
  p3 <- ggplot(g.table) +
    geom_hline(yintercept=0.5, color="black", linewidth=1.5) +
    geom_rect(aes(xmin=gene.starts, xmax=gene.end, ymin=0.499, ymax=0.501, fill=gene.strand), color="black") +
    geom_text(aes(y=0.5, x=lab.table, label=names(lab.table)),size=3.5) +
    xlab(paste0("\nGenomic position on chromosome ", chr, " (Mb)")) +
    scale_x_continuous(
      limits=c(round(bp_min/1e+6,2), round(bp_max/1e+6,2)),
      breaks=round(seq(bp_min, bp_max, by=100000)/1e+6,2),
      expand=c(0.01,0.01)
    ) +
    scale_y_continuous(limits=c(0.495,0.501), expand = c(0,0)) +
    scale_fill_manual(
      values = c("-1"="#ff8000", "1"="#007fff"),
      labels = c("-","+")
      ) +
    theme_bw() +
    guides(fill=guide_legend(title="DNA strand   ")) +
    theme(
      legend.background=element_rect(linewidth=0.5, linetype="solid", colour="black"),
      legend.text = element_text(size = 12, vjust=0.5),
      legend.title = element_text(size = 10, vjust=0.5),
      legend.position = "bottom",
      legend.direction="horizontal",
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.5, "cm"),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color="black", linewidth=0.5),
      axis.line.y = element_blank()
    )
  
  ### Align all plots and legends
  p1_legend <- get_legend(p1)
  p3_legend <- get_legend(p3)
  
  arrange_p <- plot_grid(
    p1 + theme(legend.position = "none", plot.margin = unit(c(0.2, 1, 0, 0.2), "cm")),
    p2 + theme(plot.margin = unit(c(0.2, 1, 0, 2), "cm")),
    p3 + theme(legend.position = "none", plot.margin = unit(c(0, 1, 0.5, 0.2), "cm")),
    plot_grid(p1_legend,p3_legend, ncol=2),
    align="v", ncol=1, rel_heights=c(1,0.5,0.4,0.1))
  
  ggsave(paste0(opt$output, "/plots/locus_", locus, "_results_summary_plot_2.png"),
         arrange_p, width=45, height=22, units="cm", bg="white")
} 

