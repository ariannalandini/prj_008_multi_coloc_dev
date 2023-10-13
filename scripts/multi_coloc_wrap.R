### Multicoloc
### detach all loaded packages
#lapply(paste("package:", names(sessionInfo()$otherPkgs), sep=""), 
#       detach, 
#       character.only = TRUE, 
#       unload = TRUE)


source("prj_008_multi_coloc_dev/scripts/multi_coloc_funs.R")

### Load necessary packages, if not available install them first
package_list <- c("optparse","data.table","tidyr","corrplot","coloc","bigsnpr","ggplot2","easyGgplot2","cowplot","igraph","RColorBrewer","ggnet","patchwork","stringi","reshape2","plyr","Gviz","EnsDb.Hsapiens.v75","purrr","dplyr")
for(package in package_list){
  package.loader(package)
}


option_list <- list(
  make_option("--input", type="character", default=NULL, 
              help="Pan loci table created by the locuslister R function", metavar="character"),
  make_option("--chr", type="character", default="CHROM",
              help="Name of chromosome column", metavar="integer"),
  make_option("--pos", type="character", default="GENPOS", 
              help="Name of genetic position column", metavar="character"),  
  make_option("--rsid", type="character", default="ID", 
              help="Name of rsid column", metavar="character"),  
  make_option("--a1", type="character", default="ALLELE1", 
              help="Name of effect allele column", metavar="character"),  
  make_option("--a0", type="character", default="ALLELE0", 
              help="Name of NON effect allele column", metavar="character"),  
  make_option("--freq", type="character", default="A1FREQ", 
              help="Name of effect allele frequency column", metavar="character"),  
  make_option("--n", type="character", default="N", 
              help="Name of sample size column", metavar="character"),  
  make_option("--effect", type="character", default="BETA", 
              help="Name of effect size column", metavar="character"),  
  make_option("--se", type="character", default="SE", 
              help="Name of standard error of effect column", metavar="character"), 
  make_option("--pvalue", type="character", default="P", 
              help="Name of p-value of effect column", metavar="character"),
  make_option("--output", type="character", default="./multi_coloc", 
              help="Path and name of output directory", metavar="character"),
  make_option("--grch", type="integer", default=38, 
              help="Genomic build of GWAS summary statistics", metavar="character"),
  make_option("--maf", type="numeric", default=0.0001, 
              help="MAF filter", metavar="character"),
  make_option("--save_inter_files", type="numeric", default=FALSE, 
              help="Whether to save intermediate datasets as R objects", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

########### to delete
#opt$input="/group/pirastu/prj_004_variant2function/gwas_topmed_rap/mh_and_loci/ukbb_topmed_all_loci.tsv"
#opt$output="/group/pirastu/prj_004_variant2function/coloc/multi_coloc"
###########

## Throw error message - munged GWAS summary statistics file MUST be provided!
if(is.null(opt$input)){
  stop("Please specify the file name and path of your pan loci table in --input option\n",
    call.=FALSE)
#  print_help(opt_parser)
}
  
## Create output folder
system(paste0("mkdir -p ", opt$output, "/temporary"))
system(paste0("mkdir -p ", opt$output, "/plots"))
system(paste0("mkdir -p ", opt$output, "/results"))

## Load-in pan locus table and sort by pan locus number
loci.table <- fread(opt$input)

# If pan locus tabel is produced through loci identification script, it will report the start and the end of the pan locus ONLY in the pan locus name - fix this
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

### NB: for larger pan loci, multiple loci from the same trait have been collapsed?! Doesn't make sense to munge and perform cojo more than once on the same combo of trait and locus
loci.table <- loci.table %>% select(any_of(c("chr","start","end","trait","path","pan_locus","type","sdY","s"))) %>% distinct()

## Locus defined by array job
locus <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

## Set HLA coordinates. See:
# GRCh38 - https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13
# GRCh37 - https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37

# Set start and end of HLA locus 
if(opt$grch==38){
  hla_start=28510120
  hla_end=33480577
}
if(opt$grch==37){
  hla_start=28477797
  hla_end=33448354
}

## Identify loci falling in HLA region
hla_locus <- unique((
  loci.table %>%
    filter(chr==6) %>%
    mutate(flag=data.table::between(hla_start, start, end) | data.table::between(hla_end, start, end)) %>%
    filter(flag==TRUE))$pan_locus)

## Don't run for HLA loci as cojo will take forever
if(locus %in% hla_locus){
  cat("\nLocus falling in the HLA region, colocalization not performered - script stops here\n")
} else {

## Define reference files for munging and cojo based on genomic build
  if(opt$grch==38){
    ## Reference map for munging
    mappa <- fread("/ssu/bsssu/ghrc38_reference/ukbb_grch38_map.tsv") ## temporary location?
    ## LD reference panel (30k random unrelated british UKBB)
    bfile="/ssu/bsssu/ghrc38_reference/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british"
  }
  
  if(opt$grch==37){
    ## Reference map for munging
    mappa <- fread("/ssu/bsssu/ghrc37_reference/UKBB_30k_map.tsv") ## temporary location?
    ## LD reference panel (30k random unrelated british UKBB)
    bfile="/processing_data/shared_datasets/ukbiobank/genotypes/LD_reference/ld_reference_bfiles/ukbb_all_30000_random_unrelated_white_british"
  }
  
  
## Set up starting input  
  final.locus.table=c()  
  col.order=c("trait","Chr","start","end","SNP","bp","refA","othA","freq","b","se","p","bJ","bJ_se","pJ","LD_r","n","pan.locus","sub_locus")
  
# Define genomic region
  loci.table.tmp=loci.table[loci.table$pan_locus==locus,]
  locus.info=loci.table.tmp
  
  start=min(loci.table.tmp$start)-100000
  end=max(loci.table.tmp$end)+100000
  chr=loci.table.tmp$chr[1]
  mappa.loc=mappa[which(mappa$CHR==chr & mappa$BP>=start & mappa$BP<=end),]
#  n.table=c()  ### Who uses this?
  
  cat("\nAll set and ready to start!\n")

## Munge files
  datasets <- lapply(loci.table.tmp %>% group_split(trait), function(x){
    dataset.munge(
      x$path
      ,map = mappa.loc
      ,snp.lab = opt$rsid
      ,chr.lab = opt$chr
      ,pos.lab = opt$pos
      ,a1.lab = opt$a1
      ,a0.lab = opt$a0
      ,beta.lab = opt$effect
      ,se.lab = opt$se
      ,freq.lab = opt$freq
      ,pval.lab = opt$pvalue
      ,n.lab = opt$n
      ,type = x$type ### better to add this directly to the GWAS sum stats?
      ,sdY = x$sdY ### better to add this directly to the GWAS sum stats?
      ,s = x$s
    )
  })
  names(datasets) <- unique(loci.table.tmp$trait)
  cat("\nGWAS summary statistics succesfully munged\n")
  
  if(opt$save_inter_files==TRUE){
    saveRDS(datasets, file=paste0(opt$output, "/temporary/locus_", locus, "_datasets.rds"))
#    datasets <- readRDS(file=paste0(opt$output, "/temporary/locus_", locus, "_datasets.rds"))
  }  

# Perform cojo
  conditional.datasets=list()
  max.loci=1
  
  for(i in 1:length(datasets)){
    tmp=cojo.ht(D=datasets[[i]], p.tresh=1e-4, maf.thresh=opt$maf, bfile=bfile)
    if(!is.null(tmp)){
      conditional.datasets[[i]]=tmp
      names(conditional.datasets)[i]=names(datasets)[i]
      max.loci=max(max.loci,nrow(tmp$ind.snps))
    }
  }
## Remove eventually empty dataframes (maf filter removed all independent SNPs)  
  conditional.datasets <- conditional.datasets %>% discard(is.null)
  
  if(length(conditional.datasets)==0){
    cat(paste0("No independent signal identified for any of the traits at locus ", locus, " for the p-value (1e-4) and maf (", opt$maf, ") thresholds specified" ))
  } else {
  
    cat("\nSecondary associations signals identified with COJO\n")
     
    if(opt$save_inter_files==TRUE){
      saveRDS(conditional.datasets, file=paste0(opt$output, "/temporary/locus_", locus, "_conditional.datasets.rds"))
#      conditional.datasets <- readRDS(file=paste0(opt$output, "/temporary/locus_", locus, "_conditional.datasets.rds"))
    }
    
  # Plot of all independent associations for each trait
    pdf(paste0(opt$output, "/plots/locus_", locus, "_conditioned_loci.pdf"), height=3.5*max.loci, width=10)
    for(i in 1:length(conditional.datasets)){
      p4 <- plot.cojo.ht(conditional.datasets[[i]]) +
        plot_annotation(paste("Locus", locus, names(conditional.datasets)[i]))
      print(p4)
    }
    dev.off()
  
  
  # Only if there are multiple traits at the same locus
    if(length(conditional.datasets)>1){
      
  # Identify all pairwise combination of traits to test
      pairwise.list=t(combn(names(conditional.datasets),2))
  # Prepare final locus.table for coloc
      final.locus.table.tmp <- coloc.prep.table(pairwise.list, conditional.datasets, loci.table.tmp)
  
  ### Run colocalisation and store results
      final.colocs.summary=c()
      final.colocs.results=c()
  # Loop over each pair-wise combination of traits (stored in pairwise.list)   
      for(i in 1:nrow(pairwise.list)){
  # Perform colocalisation for each combination of independent SNPs
        coloc.res=colo.cojo.ht(
          conditional.dataset1=conditional.datasets[[pairwise.list[i,1]]],
          conditional.dataset2 = conditional.datasets[[pairwise.list[i,2]]],
          p.threshold.cond = 1e-6,
          p.threshold.orig = 5e-8)
  # Check if coloc was actually performed        
        if(!is.null(coloc.res)){
  # Store the summary output in a data frame, adding tested traits column         
          only_summary_df <- as.data.frame(rbindlist(lapply(coloc.res, function(x) {
            x$summary <- x$summary %>% 
              mutate(t1=pairwise.list[i,1], t2=pairwise.list[i,2])
          }))) %>% mutate(pan.locus=locus)
  # Store the results output in a list, adding tested traits column     
          only_results_list <- lapply(coloc.res, function(x) {
            x$results <- x$results %>% 
              mutate(t1=pairwise.list[i,1], t2=pairwise.list[i,2], pan.locus=locus)
          })
  # Append loop result
            final.colocs.summary=rbind(final.colocs.summary,only_summary_df)
            final.colocs.results=c(final.colocs.results,only_results_list)
        }
      }
  
  # Get the index of columns where PP.H4 >= 0.9
      index <- which(round(as.numeric(final.colocs.summary$PP.H4.abf),2) >= 0.90)
  
  # Keep only traits colocalising (PP.H4 >= 0.9) for both summary and results coloc output
      final.colocs.H4 <- final.colocs.summary[index,]
      by_snp_PPH4 <- final.colocs.results[index]
      by_snp_PPH3 <- final.colocs.results[setdiff(seq(1,length(final.colocs.results)), index)]
      
   
  ### Define colocalisation groups
      colocalization.table.all=c()
  #    colocalization.table.H4=c() (needed?)
  #   k=1 (needed? Left from Nicola's code)
      cs_threshold=0.99 #### Set as user provided parameter?
          
  # If at least a couple of traits successfully colocalised 
      if(!is.null(final.colocs.H4)){
  # Create a graph from the "hit1" and "hit2" columns of the final.colocs.H4 data frame
        a.graph=graph_from_data_frame(final.colocs.H4[,c("hit1","hit2")],directed=F)
  # Identify connected components in the graph      
        groups=components(a.graph)
        groups=data.frame(snp=names(groups$membership),group=groups$membership)
  # Assign group to coloc results   
        final.colocs.H4$g1=groups$group[match(final.colocs.H4$hit1,groups$snp)]
  # Assign group to indipendent SNPs  
        final.locus.table.tmp <- coloc.subgrouping(final.colocs.H4, final.locus.table.tmp)
  
  ### Flag SNPs not passing p-value filtering condition applied to colocalisation groups
  ######## Filtering criteria to double-check with Nicola - SNPs with extremely significant pJ (but not p) are removed   
        final.locus.table.tmp <- as.data.frame(final.locus.table.tmp %>%
          group_by(sub_locus) %>%
          mutate(flag=ifelse(any((p < 5e-8 & pJ < 1e-6) | pJ < 5e-8), "keep", "remove")) %>%
          mutate(flag=ifelse((p < 5e-8 & pJ < 1e-4) | pJ < 5e-8, flag, "remove"))
        )
  
  ### Create pleiotropy table
        pleio.all=c()
        sub.loci=sort(unique(final.locus.table.tmp$sub_locus))
        
        for(i in sub.loci){
            tmp=pleio.table(conditional.datasets = conditional.datasets,
                            loc.table = final.locus.table.tmp[final.locus.table.tmp$sub_locus==i,])
            tmp$sublocus=i
            pleio.all=rbind(pleio.all,tmp)
        }
        pleio.all=unique(pleio.all)
          
        pleio.all=reshape2::dcast(pleio.all,formula = SNP+sublocus+trait~variable,fill=NA)
        pleio.all=pleio.all[order(pleio.all$sublocus), ]
        pleio.all$Z=pleio.all$b/pleio.all$se
          
        a=reshape2::dcast(pleio.all[,c("SNP","trait","Z")],SNP~trait,fill = 0)
        row.names(a)=a$SNP
                
        pdf(paste0(opt$output, "/plots/locus_",locus,"_pleiotropy_table.pdf"),
          width=ifelse((dim(a)[1]/dim(a)[2])*7>4, (dim(a)[1]/dim(a)[2])*7, 4)
        )
        corrplot(t(as.matrix(a[,-1])),
          is.corr = F,
          method = "color",
          addgrid.col = 'darkgrey',
          col=COL2('RdBu', 200),
          col.lim=c(max(abs(a[,-1]))*-1,max(abs(a[,-1]))))
        dev.off()
  
     
  ### Plot coloc
        coloc.plot(final.colocs.H4, outpath=paste0(opt$output, "/plots/"))  
      }else{
        idx=which(is.na(final.locus.table.tmp$sub_locus))
        pri=1
        final.locus.table.tmp$sub_locus[idx]=pri:(pri+length(idx)-1)
        final.locus.table.tmp=as.data.frame(final.locus.table.tmp)[,col.order]
      } 
  
  ### Save ALL colocalisation summary output    
      if(!is.null(final.colocs.summary)){
  #      final.colocs.summary$pan.locus=locus
        colocalization.table.all=rbind(colocalization.table.all,final.colocs.summary) ### Necessary??
        write.table(colocalization.table.all,
          file=paste0(opt$output, "/results/locus_", locus, "_colocalization.table.all.tsv"),
          row.names=F,quote=F,sep="\t")
      }
      
      if(!is.null(final.colocs.H4)){
  
  ### Join H4 coloc info with flagged SNPs info to remove SNPs failing above p-value filtering
  # Summary output of coloc      
        colocalization.table.H4 <- final.colocs.H4 %>%
          inner_join(final.locus.table.tmp %>% select(SNP, trait, sub_locus, flag),
          by=c("t1"="trait", "hit1"="SNP", "g1"="sub_locus")
          , multiple = "all"
        ) %>%
        inner_join(final.locus.table.tmp %>% select(SNP, trait, sub_locus, flag),
          by=c("t2"="trait", "hit2"="SNP", "g1"="sub_locus")
          ,multiple = "all"
        ) 
        
  # Get the index of rows where at least one "remove" flag is present
        index2 <- which(apply(colocalization.table.H4, 1, function(x) sum(x == "keep"))==2)
  
  # Remove all SNPs flagged and save from coloc summary output     
        colocalization.table.H4 <- colocalization.table.H4 %>%
          filter(flag.x=="keep" & flag.y=="keep") %>%
          select(-flag.x, -flag.y)
    
  # Result output of coloc
        by_snp_PPH4_final <- by_snp_PPH4[index2]
  
  ### Formatting of coloc results output for merge of same SNP across multiple traits         
        by_snp_PPH4_final <- lapply(by_snp_PPH4_final, function(x){
          x %>%
  # Add sub locus info in result output of coloc    
            left_join(colocalization.table.H4 %>% select(t1,hit1,t2,hit2,g1),
              by=c("t1","hit1","hit2","t2")) %>% 
            select(snp, matches("lABF"), t1,hit1,t2,hit2,g1,pan.locus) %>%
  # Move lABF values in a single column (SNPs duplicated by trait)          
            gather("trait", "lABF", -snp,-hit1,-hit2,-t1,-t2,-pan.locus,-g1) %>%
            mutate(trait=ifelse(trait=="lABF.df1", unique(t1), unique(t2))) %>%
            mutate(cojo_hit=ifelse(trait==t1, unique(hit1), unique(hit2))) %>%
            select(-hit1,-hit2,-t1,-t2)
        })
  
  # Merge in single data frame and then re-split by trait            
        by_snp_PPH4_final <- rbindlist(by_snp_PPH4_final) %>% group_split(trait,cojo_hit)
        
        by_snp_PPH4_final <- lapply(by_snp_PPH4_final, function(x){
          temp <- x %>% distinct(snp, .keep_all = T) %>%
            mutate(bf=exp(lABF)) %>% 
            arrange(desc(bf)) %>% 
            mutate(pp.cv = bf/sum(bf)) %>% 
            mutate(cred.set = cumsum(pp.cv)) %>%
            arrange(cred.set) ### just to be extra sure
          w <- which(temp$cred.set > cs_threshold)[1] ### do not hard code 99%
          
          ### CUMULATIVE SUM PLOT (to keep?)
          ggplot(temp %>% mutate(label=seq(1:nrow(temp))), aes(x=label, y=cred.set)) +
            geom_line() +
            geom_point(temp %>% mutate(label=seq(1:nrow(temp))) %>% slice(1:w), mapping=aes(x=label, y=cred.set, color="red"), show.legend = FALSE) +
            #          scale_x_continuous(breaks=seq(0, nrow(temp), 500)) + #### automatic breaks?
            ggtitle(paste0("99% credible set (in red) for ", unique(temp$trait)," (locus ", unique(temp$pan.locus), ", sublocus ",  unique(temp$g1), ")")) +
            xlab("Number of SNPs") + theme(axis.title.y=element_blank())
          ggsave(paste0(opt$output, "/plots/locus_", locus, "_sublocus_", unique(temp$g1), "_", unique(temp$trait), "_99_credible_set.png"), width=30, height = 12, units = "cm")
          
          temp <- temp %>% slice(1:w)
          temp
        })
        
  ### Credible set intersection
        inter_info <- rbindlist(by_snp_PPH4_final) %>%
          group_by(g1) %>%
          mutate(n_traits=length(unique(trait))) %>%
          group_by(g1,snp) %>%
          mutate(
            n_snps=length(unique(trait)),
            flag=n_traits==n_snps,
            joint.pp.cv=sum(pp.cv)
            ) %>%
          ungroup() %>%
          filter(flag==TRUE) %>%
          select(-n_traits,-n_snps, -flag) 
        
        inter <- inter_info %>%
          distinct(snp, .keep_all=T) %>%
          select(pan.locus,g1,snp,joint.pp.cv) %>%
          group_by(g1) %>%
          mutate(joint.pp.cv=joint.pp.cv/sum(joint.pp.cv)) %>%
          filter(joint.pp.cv==max(joint.pp.cv, na.rm=T)) %>%
          ungroup() #%>%
          #arrange(g1, desc(joint.pp.cv))
          
    ### Add cs to H4 coloc table
        colocalization.table.H4 <- colocalization.table.H4 %>%
          left_join(
            rbindlist(lapply(inter %>% group_split(g1), function(x){
              x %>% 
                mutate(cs=paste0(x$snp, collapse=",")) %>%
                distinct(pan.locus,g1,cs)
            })), by=c("pan.locus", "g1"))
          
        write.table(colocalization.table.H4,
          file=paste0(opt$output, "/results/locus_", locus, "_colocalization.table.H4.tsv"), row.names=F,quote=F,sep="\t")    
   
             
  ############################################################## Joint PP      
  # Merge in single data frame and then re-split by sub locus            
  #      by_snp_PPH4_final <- rbindlist(by_snp_PPH4_final) %>% group_split(g1)
  
  #### Fine-mapping likely causal variant ~~~~ MOVE TO FUNCTION(?) ~~~~~
  #      fine.mapping.table <- as.data.frame(rbindlist(lapply(by_snp_PPH4_final, function(x){
         
  #         merged_df <- x %>% 
  # Merge by SNP dataframes from the same sub locus
  #          distinct(trait, snp, lABF, .keep_all = T) %>% ## also lABF, just to check that it stays the same for the same SNP-trait pairs
  #          arrange(snp, trait) %>%
  #          mutate(lABF_exp=exp(lABF)) %>% ### exp of log to go back to ABF
  #          group_by(trait) %>% ### group by trait
  #          mutate(lABF_std=log(lABF_exp/sum(lABF_exp))) %>%
  #          group_by(snp) %>%
  #          mutate(lABF_sum=exp(sum(lABF_std))) %>%
  #          ungroup()
           
  #        temp <- merged_df %>%
  #          distinct(snp, g1, pan.locus, lABF_sum) %>%
  #          mutate(tot=sum(lABF_sum)) %>% 
  #          mutate(lABF_sum_scaled=lABF_sum/tot) %>%
  #          filter(lABF_sum_scaled==max(lABF_sum_scaled)) %>%
  #          select(snp, lABF_sum_scaled, pan.locus, g1)
          
  #        final <- temp %>% left_join(merged_df %>% select(snp, trait, pan.locus, g1),
  #          by=c("snp", "pan.locus", "g1"), multiple = "all")
          
  #        final
  #      })))
  ###############################################################
  
    }
    
  ### Final summary plot    
      if(!is.null(final.colocs.summary)){
        tryCatch({
          final.plot(locus,final.locus.table.tmp,conditional.datasets,by_snp_PPH3=by_snp_PPH3,inter=inter,output=opt$output)}, error = function(e) {
          cat("final.plot function failed for some reasons and the plot was not produced. Ask Arianna")
          print(e)
        })
      }
  
    } else {
      final.locus.table.tmp=conditional.datasets[[1]]$ind.snps
      final.locus.table.tmp$start=unique(locus.info$start)
      final.locus.table.tmp$end=unique(locus.info$end)
      final.locus.table.tmp$pan.locus=locus
      final.locus.table.tmp$sub_locus=1
      final.locus.table.tmp$freq_geno=NA
      final.locus.table.tmp$bJ=NA
      final.locus.table.tmp$bJ_se=NA
      final.locus.table.tmp$pJ=NA
      final.locus.table.tmp$LD_r=NA
#      alleles=unlist(mappa.loc[mappa.loc$SNP %in% final.locus.table.tmp$SNP,c("A1","A2")]) ### to remove if next two lines work
#      final.locus.table.tmp$othA=alleles[!(alleles%in%final.locus.table.tmp$refA)] ### to remove if next two lines work
      alleles=mappa.loc[mappa.loc$SNP %in% final.locus.table.tmp$SNP,c("A1","A2")]
      final.locus.table.tmp <- final.locus.table.tmp %>%
        mutate(othA=ifelse(refA==alleles$A1, alleles$A2, alleles$A1))
      final.locus.table.tmp$trait=names(datasets)[1]
      final.locus.table.tmp=as.data.frame(final.locus.table.tmp)
      final.locus.table.tmp=final.locus.table.tmp[,col.order]
    }
    
    final.locus.table <- as.data.frame(rbind(final.locus.table,final.locus.table.tmp)) 
    if("flag" %in% names(final.locus.table)){
      final.locus.table <- final.locus.table %>% select(-flag) }
    
    if(exists("colocalization.table.H4")){
      final.locus.table <- final.locus.table %>% left_join(rbind(
        colocalization.table.H4 %>% select(hit1,t1,pan.locus,g1,cs) %>% rename(SNP=hit1,trait=t1),
        colocalization.table.H4 %>% select(hit2,t2,pan.locus,g1,cs) %>% rename(SNP=hit2,trait=t2)
        ) %>% distinct(), by=c("SNP","trait","pan.locus","sub_locus"="g1"))
    }
    
    print(final.locus.table)
    
    write.table(final.locus.table, file=paste0(opt$output, "/results/locus_", locus, "_final_locus_table.tsv"),
      row.names=F,quote=F,sep="\t")
  }
}
cat("\n**** DONE!! ****\n")
