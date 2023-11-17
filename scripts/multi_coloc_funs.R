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


### load_and_check_input ###
load_and_check_input <- function(opt, locus){
## Throw error message - munged GWAS summary statistics file MUST be provided!
  if(is.null(opt$input)){
    stop("Please specify the file name and path of your pan loci table in --input option\n", call.=FALSE)
    #  print_help(opt_parser)
  } 
## Load-in pan locus table
  loci.table <- fread(opt$input)

## Check that the strictly required info are present in the loci table    
  minimal_info <- c("chr","start","end","path","trait","type")
  if(any(minimal_info %in% names(loci.table))==FALSE){
    stop(c("The following mandatory columns are missing from the loci table specified:\n", setdiff(minimal_info, names(loci.table))), call.=FALSE)
  }

# If pan locus table is produced through loci identification script, it will report the start and the end of the pan locus ONLY in the pan locus name - fix this
  if("pan_locus_name" %in% names(loci.table)){
    # Assign pan locus start/end to trait-specific locus start/end
    loci.table <- loci.table %>% mutate(
      start=as.numeric(gsub("(\\d+)_(\\d+)_(\\d+)", "\\2", pan_locus_name)),
      end=as.numeric(gsub("(\\d+)_(\\d+)_(\\d+)", "\\3", pan_locus_name))
    )
  }
# If pan locus table is NOT produced through loci identification script, it will not report the pan_locus index - add it
  if(!("pan_locus" %in% names(loci.table))){
    # Takes forever!  
    #  loci.table <- loci.table %>%
    #    arrange(chr,start,end) %>%
    #    group_by(chr,start,end) %>%
    #    mutate(pan_locus2=group_indices())
    loci.table <- loci.table %>% arrange(chr,start,end) %>% group_split(chr,start,end)
    loci.table <- lapply(loci.table, function(x) as.data.frame(x))
    for(i in 1:length(loci.table)){
      loci.table[[i]] <- as.data.frame(loci.table[[i]]) %>% mutate(pan_locus=i)}
    loci.table <- as.data.frame(rbindlist(loci.table))
  }
  
# Define genomic region
  loci.table.tmp <- loci.table %>% filter(pan_locus==locus)

  ## Check that strictly required info are not NA
  if(any(is.na(loci.table.tmp %>% select(all_of(minimal_info))))){
    stop(c("NA are not allowed in the following columns:\n", minimal_info), call.=FALSE)
  }
  
## Make sure genomic build (grch) info is provided
  if(is.null(opt$grch) & is.null(loci.table.tmp$grch)){
    stop("Please specify the genomic build either using the --grch option or adding a grch in your input loci table\n", call.=FALSE)
  }
  if(is.null(opt$grch) & !is.null(loci.table.tmp$grch) & any(is.na(loci.table.tmp$grch))){
    stop("Please specify the genomic build either using the --grch option or adding a grch in your input loci table\n", call.=FALSE)
  }
  if(!is.null(opt$grch) & is.null(loci.table.tmp$grch)){
    loci.table.tmp <- loci.table.tmp %>% mutate(grch=opt$grch)
  }
  if(!is.null(opt$grch) & !is.null(loci.table.tmp$grch) & any(is.na(loci.table.tmp$grch))){
    loci.table.tmp <- loci.table.tmp %>% mutate(grch=ifelse(is.na(grch), opt$grch,grch))
  }
  if(!(unique(loci.table.tmp$grch) %in% c(37,38))){
    stop("Sorry but we support only GRCh 37 and 38 at the moment!\n", call.=FALSE)
  }
## NB: all traits in the same locus MUST be in the same build! Otherwise, how could you define the start and end of the locus??
  if(length(unique(loci.table.tmp$grch))>1){
    stop("All traits having an association in the same locus MUST have the same build!\n", call.=FALSE)
  }
## If custom LD reference bfiles are not provided (either in input loci table or as argument), assign UKBB one
  if(is.null(opt$bfile) & is.null(loci.table.tmp$bfile)){
    print("Warning: since no custom LD reference was provided, defualt UKBB one will be used\n")
    loci.table.tmp <- loci.table.tmp %>% mutate(bfile=ifelse(grch==38, 
      "/ssu/bsssu/ghrc38_reference/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british",
      "/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"))
    }
  if(is.null(opt$bfile) & !is.null(loci.table.tmp$bfile) & any(is.na(loci.table.tmp$bfile))){
    print("Warning: since no custom LD reference was provided for some traits, default UKBB one will be used for those\n")
      loci.table.tmp <- loci.table.tmp %>% mutate(bfile=ifelse(is.na(bfile), 
        ifelse(grch==38, "/ssu/bsssu/ghrc38_reference/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british", "/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"), bfile))
  }
  if(!is.null(opt$bfile) & is.null(loci.table.tmp$bfile)){
    loci.table.tmp <- loci.table.tmp %>% mutate(bfile=opt$bfile)
  }
  if(!is.null(opt$bfile) & !is.null(loci.table.tmp$bfile) & any(is.na(loci.table.tmp$bfile))){
    loci.table.tmp <- loci.table.tmp %>% mutate(bfile=ifelse(is.na(bfile), opt$bfile, bfile))
  }

### NB: for larger pan loci, multiple loci from the same trait have been collapsed?! Doesn't make sense to munge and perform cojo more than once on the same combo of trait and locus
  loci.table.tmp <- loci.table.tmp %>% select(all_of(minimal_info), "pan_locus", "grch", "bfile", any_of(c("sdY","s"))) %>% distinct()
  
# Check tyoe column  
  ### Add also parallel possibility to already have type as column of the GWAS sum stats (instead of specified in the input table)   
  if(!(unique(loci.table.tmp$type) %in% c("cc", "quant"))){
    stop("Type has been defined incorrectly - it has to be either 'cc' or 'quant'", call.=FALSE)
  }
  return(loci.table.tmp)
  cat("\nInput check performed!\n")
}



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
    dataset=fread(sumstats.file, data.table=F)
  }else{
    dataset=as.data.frame(sumstats.file)
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

#### Map/plink files have unnamed SNP as CHROM:GENPOS_A1_A0, while the GWAS summary statistics as CHROM:GENPOS only

##### TEMPORARY FIX FOR SARA - NEED TO DOUBLE CHECK THIS STEP
# dataset$SNP <- gsub("(.*)_\\w+_\\w+$", "\\1", dataset$SNP)
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
#    dataset$FRQ=map$MAF[match(dataset$SNP,map$SNP)]
    stop("For the moment, frequency of effect allele MUST be provided in the GWAS summary statistics!\n", call.=FALSE)
  
### You should be able to calculate frequency from the bfiles provided (either default or custom). PROBLEM is that at this stage the SNP ids of GWAS and bfiles are still not matching! bfiles ones in fact should be the same of the map
### Find a way to fix this!    
  }
  
  if("FRQ" %in% colnames(dataset)){
    dataset$MAF=dataset$FRQ
    dataset <- dataset %>% mutate(MAF=ifelse(MAF<0.5, MAF, 1-MAF))
  }
  
  if(!is.null(n.lab) & n.lab%in%names(dataset)){
    names(dataset)[names(dataset)==n.lab]="N"
  }else{
    N_hat<-median(1/((2*dataset$MAF*(1-dataset$MAF))*dataset$SE^2),na.rm = T) ### But where does SE comes from?? No checks performed earlier
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
  
# Add type and sdY/s
  dataset$type <- type
  if(type=="cc" & !(is.null(s)) && !is.na(s)){ # && prevents to return "logical(0)" when s is null
    dataset$s=s
  } else if(type=="cc" & (is.null(s) || is.na(s))){
    #### Is this correct?? Is "s" strictly necessary for cc traits??
    stop("Please provide s, the proportion of samples who are cases")
  }

  if(type=="quant" & !(is.null(sdY)) && !is.na(sdY)){
    dataset$sdY <- sdY
  } else if(type=="quant" & (is.null(sdY) || is.na(sdY))){ #### Gives back "logical(0)" - FIX!! Append sdY to the dataset table, even if null
    dataset$sdY <- coloc:::sdY.est(dataset$varbeta, dataset$MAF, dataset$N)
  }
  
## Match with locus reference map
#  dataset <- dataset[which(dataset$SNP %in% map$SNP),]
#  dataset <- dataset[match(map$SNP,dataset$SNP),]
  dataset <- dataset %>%
    filter(
      CHR==unique(mappa.loc$CHR),
      BP >= min(mappa.loc$BP, na.rm=T) & BP <= max(mappa.loc$BP, na.rm=T)
    )

  flip=dataset[,c("SNP","CHR","BP","A2","A1","BETA")]
  names(flip)=c("rsid","chr","pos","a0","a1","beta")
#  names(map)=c("rsid","chr","pos","maf","a1","a0") ### Do we really need to provide MAF?!
  names(map)=c("rsid","chr","pos","a1","a0")
  
  flip.t=snp_match(sumstats=flip,
                   info_snp=map,
#                   join_by_pos=FALSE,
                   remove_dups = TRUE,
                   join_by_pos=TRUE,
                   strand_flip=FALSE,
                   match.min.prop=0)
    
  #dataset=dataset[match(flip.t$rsid,dataset$SNP),]
  dataset <- dataset[flip.t$`_NUM_ID_.ss`,]

# Keep both original and map SNP id
  dataset$snp_map <- flip.t$rsid
  dataset$A1=flip.t$a1
  dataset$A2=flip.t$a0
  dataset$b=flip.t$beta

  if(type=="cc"){
    dataset <- dataset %>%
      select("snp_map","SNP","CHR","BP","A1","A2","b","varbeta","SE","P","MAF","N","type","s") %>%
      rename(se=SE, p=P)
  } else if(type=="quant"){
    dataset <- dataset %>%
      select("snp_map","SNP","CHR","BP","A1","A2","b","varbeta","SE","P","MAF","N","type","sdY") %>%
      rename(se=SE, p=P)
  }
  dataset
}



### cojo.ht ###
### Performs --cojo-slct first to identify all independent SNPs and --cojo-cond then to condition upon identified SNPs
cojo.ht=function(D=datasets[[1]]
                 ,plink.bin="/ssu/gassu/software/plink/2.00_20211217/plink2"
                 ,gcta.bin="/ssu/gassu/software/GCTA/1.94.0beta/gcta64"
                 ,bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/p01_output/ukbb_all_30000_random_unrelated_white_british"
                 ,p.tresh=1e-4
                 ,maf.thresh=0.0001){
  
  random.number=stri_rand_strings(n=1, length=20, pattern = "[A-Za-z0-9]")
  
  write(D$SNP,ncol=1,file=paste0(random.number,".snp.list"))
  system(paste0(plink.bin," --bfile ",bfile," --extract ",random.number,".snp.list --maf ", maf.thresh, " --make-bed --geno-counts --out ",random.number))
  
  freqs=fread(paste0(random.number,".gcount"))
  freqs$FreqREF=(freqs$HOM_REF_CT*2+freqs$HET_REF_ALT_CTS)/(2*(rowSums(freqs[,c("HOM_REF_CT", "HET_REF_ALT_CTS", "TWO_ALT_GENO_CTS")])))  #### Why doing all this when plink can directly calculate it with --frq?
  
# Assign allele frequency from the LD reference  
  D <- D %>%
    left_join(freqs %>% select(ID,FreqREF,REF), by=c("SNP"="ID")) %>%
    mutate(FREQ=ifelse(REF==A1, FreqREF, (1-FreqREF))) %>%
    select(-FreqREF,-REF)
### Following code was causing a lot of mismatch in allele frequency between GWAS and reference  
#  D$FREQ=freqs$FreqA1[match(D$SNP,freqs$ID)]
#  idx=which(D$A1!=freqs$REF)
#  D$FREQ[idx]=1-D$FREQ[idx]
#  D$se=sqrt(D$varbeta) # why not keeping se from the munging? se is anyway required to calculate varbeta
  D <- D %>% select("SNP","A1","A2","FREQ","b","se","p","N","snp_map","type", any_of(c("sdY","s")))
  write.table(D,file=paste0(random.number,"_sum.txt"),row.names=F,quote=F,sep="\t")

# step1 determine independent snps
  system(paste0(gcta.bin," --bfile ", random.number, " --cojo-p ", p.tresh, " --maf ", maf.thresh, " --extract ", random.number, ".snp.list --cojo-file ", random.number, "_sum.txt --cojo-slct --out ", random.number, "_step1"))
  
  if(file.exists(paste0(random.number,"_step1.jma.cojo"))){
    ind.snp=fread(paste0(random.number,"_step1.jma.cojo")) %>%
      left_join(D %>% select(SNP,snp_map,type,any_of(c("sdY", "s"))), by="SNP")

    dataset.list=list()
    dataset.list$ind.snps <- data.frame(matrix(ncol = ncol(ind.snp), nrow = 0))
    colnames(dataset.list$ind.snps) <- colnames(ind.snp)
    dataset.list$results=list()
  
    if(nrow(ind.snp)>1){
      for(i in 1:nrow(ind.snp)){
      
        write(ind.snp$SNP[-i],ncol=1,file=paste0(random.number,"_independent.snp"))
        print(ind.snp$SNP[-i])
      
        system(paste0(gcta.bin," --bfile ",random.number, " --maf ", maf.thresh, " --extract ",random.number,".snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))
      
#### STOP ANALYSIS FOR THAT TOP SNP IN CASE OF COLLINEARITY
        if(!file.exists(paste0(random.number,"_step2.cma.cojo"))){
          cat(paste0("\n****WARNING: COJO has encountered a collinearty problem. Affected SNP will be removed from following analysis****\n\n"))
        } else {
############
  # Re-add type and sdY/s info, and map SNPs!
          step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
            left_join(D %>% select(SNP,snp_map,type,any_of(c("sdY", "s"))), by="SNP")
# Add SNPs to the ind.snps dataframe         
          dataset.list$ind.snps <- rbind(dataset.list$ind.snps, ind.snp[i,])
# Add conditioned gwas to the results list          
          dataset.list$results[[i]]=step2.res
          names(dataset.list$results)[i]=ind.snp$snp_map[i]
          system(paste0("rm ",random.number,"_step2.cma.cojo"))
        }
      }
    } else {

### NB: COJO here is performed ONLY for formatting sakes - No need to condition if only one signal is found!!        
    write(ind.snp$SNP[1],ncol=1,file=paste0(random.number,"_independent.snp"))
    system(paste0(gcta.bin," --bfile ",random.number," --cojo-p ",p.tresh, " --maf ", maf.thresh, " --extract ",random.number,".snp.list --cojo-file ",random.number,"_sum.txt --cojo-cond ",random.number,"_independent.snp --out ",random.number,"_step2"))

    step2.res <- fread(paste0(random.number, "_step2.cma.cojo"), data.table=FALSE) %>%
      left_join(D %>% select(SNP,snp_map,A1,type,any_of(c("sdY", "s"))), by=c("SNP", "refA"="A1"))

#### Add back top SNP, removed from the data frame with the conditioning step
    step2.res <- rbind.fill(step2.res, ind.snp %>% select(-bJ,-bJ_se,-pJ,-LD_r))
    step2.res$bC <- NA
    step2.res$bC_se <- NA
    step2.res$pC <- NA

    dataset.list$ind.snps <- rbind(dataset.list$ind.snps, ind.snp[1,])
    dataset.list$results[[1]]=step2.res
    names(dataset.list$results)[1]=ind.snp$snp_map[1]
    }
  }
# Remove results df possibly empty (in case of collinearity issue)
  dataset.list$results <- dataset.list$results %>% discard(is.null)
  system(paste0("rm *",random.number,"*"))
  if(exists("dataset.list")){dataset.list}
}



### plot.cojo.ht ###
plot.cojo.ht=function(cojo.ht.obj){
#  library(ggplot2)
#  library(patchwork)
  
  if(nrow(cojo.ht.obj$ind.snps)>1){
    
    whole.dataset=c()
    for(i in 1:nrow(cojo.ht.obj$ind.snps)){
      
      tmp=cojo.ht.obj$results[[i]]
      tmp$signal=cojo.ht.obj$ind.snps$snp_map[i]
      whole.dataset=rbind(whole.dataset,tmp)
    }
    
    p1 <- ggplot(cojo.ht.obj$results[[i]], aes(x=bp,y=-log10(p))) +
      geom_point(alpha=0.6,size=3)+
      theme_minimal()+
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=snp_map),size=6,shape=23) +
      guides(fill=guide_legend(title="SNP"))
    
    p2 <- ggplot(whole.dataset,aes(x=bp,y=-log10(pC),color=signal)) +
      facet_grid(signal~.) +
      geom_point(alpha=0.8,size=3) +
      theme_minimal() +
      ggtitle("Conditioned results")
    
    p3 <- p1/p2 + plot_layout(heights = c(1, nrow(cojo.ht.obj$ind.snps)+0.2))
    
  } else {
    
    p3 <- ggplot(cojo.ht.obj$results[[1]], aes(x=bp,y=-log10(p))) +
      geom_point(alpha=0.6,size=3)+
      theme_minimal()+
      geom_point(data=cojo.ht.obj$ind.snps,aes(x=bp,y=-log10(p),fill=snp_map),size=6,shape=23)
  }
  (p3)
}



### coloc.prep.table
coloc.prep.table=function(pairwise.list, conditional.datasets, loci.table.tmp,mappa.loc){
  ### Select only 1) genome-wide significant or 2) with conditioned p-value < 1e-6 SNPs
  final.locus.table.tmp <- as.data.frame(rbindlist(lapply(names(conditional.datasets), function(x){
    tmp <- conditional.datasets[[x]]$ind.snps
#    if(nrow(tmp)>1){tmp <- tmp %>% filter(p<5e-8 | pJ<1e-6)} # BUT WHY?! Anyway filtering later and here is fucking up everything
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
      alleles=unlist(mappa.loc[mappa.loc$SNP==tmp$snp_map[n],c("A1","A2")])
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
  
#  if(length(grep("pJ",names(conditional.dataset1$ind.snps)))>0){
  if(any(!is.na(conditional.dataset1$ind.snps$pJ))){
    hits.t1=conditional.dataset1$ind.snps$snp_map[conditional.dataset1$ind.snps$pJ<p.threshold.cond | conditional.dataset1$ind.snps$p <p.threshold.orig]
  }else{
    hits.t1=conditional.dataset1$ind.snps$snp_map
  }
  
#  if(length(grep("pJ",names(conditional.dataset2$ind.snps)))>0){
  if(any(!is.na(conditional.dataset2$ind.snps$pJ))){
    hits.t2=conditional.dataset2$ind.snps$snp_map[conditional.dataset2$ind.snps$pJ<p.threshold.cond | conditional.dataset2$ind.snps$p<p.threshold.orig]
  }else{
    hits.t2=conditional.dataset2$ind.snps$hits.t2
  }
  
  if(length(hits.t2)>0 & length(hits.t1)>0){
    
    coloc.final=list()
    for(i in hits.t1){
      for(j in hits.t2){
        
        D1=conditional.dataset1$results[[i]]     
        D2=conditional.dataset2$results[[j]] 
        
        if(any(!is.na(D1$bC))){
          D1 <- D1 %>%
            select("snp_map","Chr","bp","bC","bC_se","n","pC","freq","type",any_of(c("sdY","s")),"SNP") %>%
            rename("snp"="snp_map","chr"="Chr","position"="bp","beta"="bC","varbeta"="bC_se","N"="n","pvalues"="pC","MAF"="freq","snp_original"="SNP")
        } else {
          D1 <- D1 %>%
            select("snp_map","Chr","bp","b","se","n","p","freq","type",any_of(c("sdY","s")),"SNP") %>%
            rename("snp"="snp_map","chr"="Chr","position"="bp","beta"="b","varbeta"="se","N"="n","pvalues"="p","MAF"="freq","snp_original"="SNP")
        }
        D1$varbeta=D1$varbeta^2
        D1=na.omit(D1)
        
        if(any(!is.na(D2$bC))){
          D2 <- D2 %>%
            select("snp_map","Chr","bp","bC","bC_se","n","pC","freq","type",any_of(c("sdY","s")),"SNP") %>%
            rename("snp"="snp_map","chr"="Chr","position"="bp","beta"="bC","varbeta"="bC_se","N"="n","pvalues"="pC","MAF"="freq","snp_original"="SNP")
        }else{
          D2 <- D2 %>%
            select("snp_map","Chr","bp","b","se","n","p","freq","type",any_of(c("sdY","s")),"SNP") %>%
            rename("snp"="snp_map","chr"="Chr","position"="bp","beta"="b","varbeta"="se","N"="n","pvalues"="p","MAF"="freq","snp_original"="SNP")
        }
        D2$varbeta=D2$varbeta^2
        D2=na.omit(D2)
        
        colo.res <- coloc.abf.ht(D1,D2)
## Save coloc summary        
        colo.sum=data.frame(t(colo.res$summary))
        colo.sum$hit1=i
        colo.sum$hit2=j
#       coloc.summary=rbind(coloc.summary,colo.sum)

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
coloc.subgrouping <- function(final.colocs.H4, final.locus.table.tmp, col.order){
  # Add sublocus to locus table
  for(k in 1:nrow(final.colocs.H4)){
    
    final.locus.table.tmp$sub_locus[
      final.locus.table.tmp$trait==final.colocs.H4$t1[k] &
        final.locus.table.tmp$snp_map==final.colocs.H4$hit1[k]]=final.colocs.H4$g1[k]
    
    final.locus.table.tmp$sub_locus[
      final.locus.table.tmp$trait==final.colocs.H4$t2[k] &
        final.locus.table.tmp$snp_map==final.colocs.H4$hit2[k]]=final.colocs.H4$g1[k]
  }
  
  # Add sublocus also for non-colocalising loci
  if(any(is.na(final.locus.table.tmp$sub_locus))){
    idx=which(is.na(final.locus.table.tmp$sub_locus))
    pri=max(final.locus.table.tmp$sub_locus,na.rm=T)+1
    final.locus.table.tmp$sub_locus[idx]=pri:(pri+length(idx)-1)
  }
  
  final.locus.table.tmp <- final.locus.table.tmp %>% rename(snp_original=SNP, SNP=snp_map)
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
      
      if(any(!is.na(tmp$pC))){
        tmp=tmp[,c("bp","pC","label","group")]
        names(tmp)=c("bp","p","label","group")
      } else {
        tmp=tmp[,c("bp","p","label","group")]
      }  
      
      per.plot.data=rbind(per.plot.data,tmp)
      tmp=conditional.datasets[[x$t2[i]]]$results[[x$hit2[i]]]
      tmp$label=paste(x$t2[i],x$hit2[i],sep="-")
      tmp$group=x$g1[i]
      
      if(any(!is.na(tmp$pC))){
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



### pleio.table ###
pleio.table=function(conditional.datasets=conditional.datasets,loc.table=NA,plot=FALSE,plot.file=NULL, index.trait=NULL){
  
## Remove duplicates
  duplicati=table(loc.table$trait)
  if(any(duplicati>1)){
    doppi=names(duplicati)[duplicati>1]
    for(i in doppi){
      tmp=conditional.datasets[[i]]$ind.snps
      tmp=tmp[which(tmp$snp_map%in%loc.table$SNP[loc.table$trait==i]),]
      min.snp=tmp$snp_map[which.min(tmp$pJ)]
      loc.table=loc.table[which(loc.table$trait!=i | (loc.table$trait==i & loc.table$SNP==min.snp)),]
    }
  }
  
  if(nrow(loc.table)>1){
    for(i in 1:nrow(loc.table)){
      
      if(i==1){
        merged.datasets=conditional.datasets[[loc.table$trait[i]]]$results[[loc.table$SNP[i]]]
        
        if(any(!is.na(merged.datasets$bC))){
          merged.datasets=merged.datasets[,c("snp_map","bC","bC_se","pC","n")]
          names(merged.datasets)=c("SNP",paste(c("beta","se","pval","n"),loc.table$trait[i],sep=".."))
        } else {
          merged.datasets=merged.datasets[,c("snp_map","b","se","p","n")]
          names(merged.datasets)=c("SNP",paste(c("beta","se","pval","n"),loc.table$trait[i],sep=".."))
        }
        
      } else {
        
        merged.datasets2=conditional.datasets[[loc.table$trait[i]]]$results[[loc.table$SNP[i]]]
        
        if(any(!is.na(merged.datasets2$bC))){
          merged.datasets2=merged.datasets2[,c("snp_map","bC","bC_se","pC","n")]
          names(merged.datasets2)=c("SNP",paste(c("beta","se","pval","N"),loc.table$trait[i],sep=".."))
        } else {
          merged.datasets2=merged.datasets2[,c("snp_map","b","se","p","n")]
          names(merged.datasets2)=c("SNP",paste(c("beta","se","pval","n"),loc.table$trait[i],sep=".."))
        }
        merged.datasets=merge(merged.datasets,merged.datasets2,by="SNP",all = FALSE,suffixes =c(""))
      }
    }
    
    merged.datasets=na.omit(merged.datasets)
    if(is.null(index.trait)){
      top.snp=as.data.frame(merged.datasets[which.min(apply(merged.datasets[,grep("pval..",names(merged.datasets))],1,min)),])
    } else if (length(grep(index.trait,names(merged.datasets)))==0){
      top.snp=as.data.frame(merged.datasets[which.min(apply(merged.datasets[,grep("pval..",names(merged.datasets))],1,min)),])
    } else {
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
    
  } else {
    
    merged.datasets=conditional.datasets[[loc.table$trait[1]]]$results[[loc.table$SNP[1]]]
    
    if(any(!is.na(merged.datasets$bC))){
      merged.datasets=merged.datasets[,c("snp_map","bC","bC_se","pC","n")]
      names(merged.datasets)=c("SNP","beta","se","pval","n")
    } else {
      merged.datasets=merged.datasets[,c("snp_map","b","se","p","n")]
      names(merged.datasets)=c("SNP","beta","se","pval","n")
    }
    top.snp=merged.datasets[which.min(merged.datasets$p),]
    names(top.snp)=c("SNP",paste(c("beta","se","pval","n"),loc.table$trait[1],sep=".."))
    top.snp=reshape2::melt(top.snp)
    matric=matrix(unlist(strsplit(as.character(top.snp$variable),split="\\.\\.")),ncol=2,byrow = T)
    top.snp$variable=matric[,1]
    top.snp$trait=matric[,2]
  }
  top.snp
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
                       data_sub,
                       finemap,
                       output=opt$output
){

# Integrate finemapping (lABF) and GWAS (beta...and p-value?) info
  full_df <- as.data.frame(rbindlist(lapply(names(data_sub), function(x){
    left_join(finemap[[x]],
              data_sub[[x]] %>% select(SNP,Chr,bp,any_of(c("b","bC")),any_of(c("p","pC")),trait,cojo_snp),
              by=c("snp"="SNP","trait","cojo_snp"))
  }), fill=TRUE))
  
# Add sublocus info to conditional dataset - keep only SNPs actually submitted to coloc
  loci_table <- full_df %>%
    right_join(final.locus.table.tmp %>% filter(flag=="keep") %>% select(sub_locus, trait, SNP),
      by=c("trait", "cojo_snp"="SNP"))

# Find representative SNP to plot
  final <- as.data.frame(rbindlist(lapply(loci_table %>% group_split(sub_locus), function(x){
    x %>%
      group_by(snp) %>%
      mutate(joint.pp=sum(lABF)) %>%
      ungroup() %>%
      filter(joint.pp==max(joint.pp))
  })))

## Adjustment for plotting
  if("bC" %in% names(final)){
    final <- final %>% mutate(beta=ifelse(is.na(bC), b, bC))
  } else {
    final <- final %>% mutate(beta=b)
  } ### All this should be automatically fixed by adjustment I made to cojo output in the dev branch!!
  final <- final %>%
    mutate(dir=ifelse(beta<0, "-", "+")) %>%
    mutate(group=paste0(sub_locus, " ", snp)) %>%
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
  
### If SNP labels are too many, plot them 45 degrees angles
  if(length(unique(final$group))>10){
    p1 <- p1 + theme(axis.text.x = element_text(size=9, angle=45, hjust=1))
  }
  
  
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
  
# Plot only if at least one gene is found  
  if(length(a$gene_id)>0){
  
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
    
    d <- do.call(rbind, d) %>% mutate(y=as.numeric(ifelse(strand=="+", 1, -1)))

####### In case of the same gene having transcripts on both strands, take the strand with most transcript
      d <- d %>% semi_join(
        d %>% count(symbol, y) %>% group_by(symbol) %>% slice(which.max(n)),
        by = c("symbol", "y"))

    gene.starts <- by(d$start,d$symbol,FUN = min)
    gene.end <- by(d$end,d$symbol,FUN = max)
    gene.strand <- by(d$y,d$symbol,FUN = unique)
    
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

# If no genes has been found in the region  
  } else {
    g.table <- data.frame(chr=chr, position=c(bp_min,bp_max))

    p3 <- ggplot(g.table, aes(position,1)) +
      geom_blank() +
      geom_hline(yintercept=0.5, color="black", linewidth=1.5) +
      xlab(paste0("\nGenomic position on chromosome ", chr, " (Mb)")) +
      scale_x_continuous(
        limits=c(round(bp_min/1e+6,2), round(bp_max/1e+6,2)),
        breaks=round(seq(bp_min, bp_max, by=100000)/1e+6,2),
        expand=c(0.01,0.01)
      ) +
      scale_y_continuous(limits=c(0.495,0.501), expand = c(0,0)) +
      theme_bw() +
      theme(
        legend.background=element_rect(linewidth=0.5, linetype="solid", colour="black"),
        legend.text = element_text(size = 12, vjust=0.5),
        legend.title = element_text(size = 10, vjust=0.5),
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

    arrange_p <- plot_grid(
      p1 + theme(legend.position = "none", plot.margin = unit(c(0.2, 1, 0, 0.2), "cm")),
      p2 + theme(plot.margin = unit(c(0.2, 1, 0, 2), "cm")),
      p3 + theme(legend.position = "none", plot.margin = unit(c(0, 1, 0.5, 0.2), "cm")),
      plot_grid(p1_legend, ncol=1),
      align="v", ncol=1, rel_heights=c(1,0.5,0.4,0.1))
  }
  
  ggsave(paste0(opt$output, "/plots/locus_", locus, "_results_summary_plot.png"),
         arrange_p, width=45, height=22, units="cm", bg="white")
} 


#################### Coloc functions modified to accept dataframes (and not only lists)

#### coloc.abf.ht
coloc.abf.ht <- function(dataset1, dataset2, MAF=NULL, 
                         p1=1e-4, p2=1e-4, p12=1e-5) {
  
  if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF
  check_dataset.ht(d=dataset1,1)
  check_dataset.ht(d=dataset2,2)
  
  df1 <- process.dataset.ht(d=dataset1, suffix="df1")
  df2 <- process.dataset.ht(d=dataset2, suffix="df2")
  merged.df <- merge(df1,df2)

  if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")
  
  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  my.denom.log.abf <- coloc:::logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)
  
  pp.abf <- coloc:::combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)
  
  output<-list(summary=results,
               results=merged.df,
               priors=c(p1=p1,p2=p2,p12=p12))
  class(output) <- c("coloc_abf",class(output))
  return(output)
}


#### check_dataset.ht
check_dataset.ht <- function (d, suffix = "", req = c("snp"), warn.minp = 1e-06) 
{
  if (!is.list(d)) 
    stop("dataset ", suffix, ": is not a list")
  nd <- names(d)
  n <- 0
  for (v in nd) {
    if (v %in% req && !(v %in% nd)) 
      stop("dataset ", suffix, ": missing required element ", 
           v)
    if (any(is.na(d[[v]]))) 
      stop("dataset ", suffix, ": ", v, " contains missing values")
  }
  if ("snp" %in% nd && any(duplicated(d$snp))) 
    stop("dataset ", suffix, ": duplicated snps found")
  if ("snp" %in% nd && is.factor(d$snp)) 
    stop("dataset ", suffix, ": snp should be a character vector but is a factor")
  if ("MAF" %in% nd && (!is.numeric(d$MAF) || any(is.na(d$MAF)) || 
                        any(d$MAF <= 0) || any(d$MAF >= 1))) 
    stop("dataset ", suffix, ": MAF should be a numeric, strictly >0 & <1")
  l <- -1
  shouldmatch <- c("pvalues", "MAF", "beta", "varbeta", "snp", 
                   "position")
  for (v in shouldmatch) if (v %in% nd) 
    if (l < 0) {
      l <- length(d[[v]])
    }
  else {
    if (length(d[[v]]) != l) {
      stop("dataset ", suffix, ": lengths of inputs don't match: ")
      print(intersect(nd, shouldmatch))
    }
  }
  if (("N" %in% req) && (!("N" %in% nd) || is.null(d$N) || 
                         is.na(d$N))) 
    stop("dataset ", suffix, ": sample size N not set")
  if (!("type" %in% nd)) 
    stop("dataset ", suffix, ": variable type not set")
  if (!(d$type[[1]] %in% c("quant", "cc"))) 
    stop("dataset ", suffix, ": ", "type must be quant or cc")
  if (("s" %in% nd) && (!is.numeric(d$s) || d$s <= 0 || d$s >= 
                        1)) 
    stop("dataset ", suffix, ": ", "s must be between 0 and 1")
  if (!("beta" %in% nd) || !("varbeta" %in% nd)) {
    if (!("pvalues" %in% nd) || !("MAF" %in% nd)) 
      stop("dataset ", suffix, ": ", "require p values and MAF if beta, varbeta are unavailable")
    if (d$type[[1]] == "cc" && !("s" %in% nd)) 
      stop("dataset ", suffix, ": ", "require, s, proportion of samples who are cases, if beta, varbeta are unavailable")
    p = d$pvalues
  }
  else {
    p = pnorm(-abs(d$beta/sqrt(d$varbeta))) * 2
  }
  if (min(p) > warn.minp) 
    warning("minimum p value is: ", format.pval(min(p)), 
            "\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check the 02_data vignette.")
  if (d$type[[1]] == "quant" && !("sdY" %in% nd)) 
    if (!("MAF" %in% nd && "N" %in% nd)) 
      stop("dataset ", suffix, ": ", "must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
  if ("LD" %in% nd) {
    if (nrow(d$LD) != ncol(d$LD)) 
      stop("LD not square")
    if (!identical(colnames(d$LD), rownames(d$LD))) 
      stop("LD rownames != colnames")
    if (length(setdiff(d$snp, colnames(d$LD)))) 
      stop("colnames in LD do not contain all SNPs")
  }
  NULL
}


#### process.dataset.ht
process.dataset.ht <- function(d, suffix) {
  nd <- names(d)
  if (!"type" %in% nd) 
    stop("dataset ", suffix, ": ", "The variable type must be set, otherwise the Bayes factors cannot be computed")
  if (!(d$type[[1]] %in% c("quant", "cc"))) 
    stop("dataset ", suffix, ": ", "type must be quant or cc")
  if (d$type[[1]] == "cc" & "pvalues" %in% nd) {
    if (!("s" %in% nd)) 
      stop("dataset ", suffix, ": ", "please give s, proportion of samples who are cases, if using p values")
    if (!("MAF" %in% nd)) 
      stop("dataset ", suffix, ": ", "please give MAF if using p values")
    if (d$s[[1]] <= 0 || d$s[[1]] >= 1) 
      stop("dataset ", suffix, ": ", "s must be between 0 and 1")
  }
  if (d$type[[1]] == "quant") {
    if (!("sdY" %in% nd || ("MAF" %in% nd && "N" %in% nd))) 
      stop("dataset ", suffix, ": ", "must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
  }
  if ("beta" %in% nd && "varbeta" %in% nd) {
    if (length(d$beta) != length(d$varbeta)) 
      stop("dataset ", suffix, ": ", "Length of the beta vectors and variance vectors must match")
    if (!("snp" %in% nd)) 
      d$snp <- sprintf("SNP.%s", 1:length(d$beta))
    if (length(d$snp) != length(d$beta)) 
      stop("dataset ", suffix, ": ", "Length of snp names and beta vectors must match")
    if (d$type[[1]] == "quant" && !("sdY" %in% nd)) 
      d$sdY <- coloc:::sdY.est(d$varbeta, d$MAF, d$N)
    df <- coloc:::approx.bf.estimates(z = d$beta/sqrt(d$varbeta), 
                              V = d$varbeta, type = d$type[[1]], suffix = suffix, sdY = d$sdY[[1]])
    df$snp <- as.character(d$snp)
    if ("position" %in% nd) 
      df <- cbind(df, position = d$position)
    return(df)
  }
  if ("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) {
    if (length(d$pvalues) != length(d$MAF)) 
      stop("Length of the P-value vectors and MAF vector must match")
    if (!("snp" %in% nd)) 
      d$snp <- sprintf("SNP.%s", 1:length(d$pvalues))
    df <- data.frame(pvalues = d$pvalues, MAF = d$MAF, N = d$N, 
                     snp = as.character(d$snp))
    snp.index <- which(colnames(df) == "snp")
    colnames(df)[-snp.index] <- paste(colnames(df)[-snp.index], 
                                      suffix, sep = ".")
    keep <- which(df$MAF > 0 & df$pvalues > 0)
    df <- df[keep, ]
    abf <- coloc:::approx.bf.p(p = df$pvalues, f = df$MAF, type = d$type[[1]], 
                       N = df$N, s = d$s, suffix = suffix)
    df <- cbind(df, abf)
    if ("position" %in% nd) 
      df <- cbind(df, position = d$position[keep])
    return(df)
  }
  stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}

################# coloc::finemap.abf() from coloc_5.1.0 was causing in some cases to have NA SNP.PP. Bug seem fixed in latest version (coloc_5.2.2)

adjust_prior=function(p,nsnps,suffix="") {
  if(nsnps * p >= 1) { ## for very large regions
    warning(paste0("p",suffix," * nsnps >= 1, setting p",suffix,"=1/(nsnps + 1)"))
    1/(nsnps + 1)
  } else {
    p
  }
}

finemap.abf.new <- function(dataset, p1=1e-4) {
  
  coloc:::check_dataset(dataset,"")
  
  df <- coloc:::process.dataset(d=dataset, suffix="")
  nsnps <- nrow(df)
  p1=adjust_prior(p1,nsnps,"1")
  
  dfnull <- df[1,]
  for(nm in colnames(df))
    dfnull[,nm] <- NA
  dfnull[,"snp"] <- "null"
  dfnull[,"lABF."] <- 0
  df <- rbind(df,dfnull)
  ## data.frame("V."=NA,
  ##            z.=NA,
  ##            r.=NA,
  ##            lABF.=1,
  ##            snp="null"))
  df$prior <- c(rep(p1,nsnps),1-nsnps*p1)
  
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  ## BUGFIX 16/5/19
  ## my.denom.log.abf <- coloc:::logsum(df$lABF + df$prior)
  ## df$SNP.PP <- exp(df$lABF - my.denom.log.abf)
  my.denom.log.abf <- coloc:::logsum(df$lABF + log(df$prior))
  df$SNP.PP <- exp(df$lABF + log(df$prior) - my.denom.log.abf)
  
  return(df)
}

