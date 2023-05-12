### Multicoloc
library(optparse)
#library(corrplot)
source("/group/pirastu/prj_004_variant2function/scripts/multi_coloc_funs.R")

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
              help="Path and name of output directory", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

########### to delete
#opt$input="/group/pirastu/prj_004_variant2function/gwas_topmed_rap/loci_and_mh/ukbb_topmed_all_loci.tsv"
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

## Load-in pan_locus table and sort by pan locus number
loci.table <- fread(opt$input) %>%
  arrange(pan_locus) %>%
## Assign pan locus start/end to trait-specific locus start/end
  mutate(
    start=as.numeric(gsub("(\\d+)_(\\d+)_(\\d+)", "\\2", pan_locus_name)),
    end=as.numeric(gsub("(\\d+)_(\\d+)_(\\d+)", "\\3", pan_locus_name))
  )

## Locus defined by array job
locus <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

## Identify loci falling in HLA region
## Coordinates taken from this paper https://www.sciencedirect.com/science/article/pii/S1357272520301990
hla_locus <- unique((
  loci.table %>%
    filter(chr==6) %>%
    mutate(flag=data.table::between(28510120, start, end) | data.table::between(33480577, start, end)) %>%
    filter(flag==TRUE))$pan_locus)

## Don't run for HLA loci as cojo will take forever - possibly add other regions?
if(locus %in% hla_locus){
  cat("\nLocus falling in the HLA region, colocalization not performered - script stops here\n")
} else {

## Reference map for munging
  mappa <- fread("/group/pirastu/prj_004_variant2function/coloc/ghrc38_reference/ukbb_grch38_map.tsv")
  
## GWAS info
  file.list <- unique(loci.table$path)
  labels <- unique(loci.table$trait)
  
  #proportions=readRDS("/project/aid_sharing/AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/04_sumstat_outputs/case_controls_fraction_nic.RDS")
  
  t.type=rep("quant", length(labels))
  prev.var=rep(1, length(labels))
  
  file.table <- data.frame(file=file.list, trait=labels, type=t.type, prev.var=prev.var)
  col.order=c("trait","Chr","start","end","SNP","bp","refA","othA","freq","b","se","p","bJ","bJ_se","pJ","LD_r","n","pan.locus","sub_locus")
  final.locus.table=c()
  
  loci.table.tmp=loci.table[loci.table$pan_locus==locus,]
  locus.info=loci.table.tmp
    
# Define genomic region
  start=min(loci.table.tmp$start)-100000
  end=max(loci.table.tmp$end)+100000
  chr=loci.table.tmp$chr[1]
  n.table=c()
  
  cat("\nAll set and ready to start!\n")

# Munge files
  mappa.loc=mappa[which(mappa$CHR==chr & mappa$BP>=start & mappa$BP<=end),]
  
  datasets=list()
  for(i in unique(loci.table.tmp$trait)){
    datasets[[i]]=dataset.munge(
    sumstats.file=file.table$file[file.table$trait==i]
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
      ,type = file.table$type[file.table$trait==i]
      ,sdY = file.table$prev.var[file.table$trait==i]
      ,s = file.table$prev.var[file.table$trait==i]
    )
  }
  cat("\nGWAS summary statistics succesfully munged\n")
  
############################################################## To delete (?)
  saveRDS(datasets, file=paste0(opt$output, "/temporary/locus_", locus, "_datasets.RData"))
  #datasets <- readRDS(file=paste0(opt$output, "/temporary/locus_", locus, "_datasets.RData"))
############################################################## 
  
## LD reference panel (30k random unrelated british UKBB)
  bfile="/group/pirastu/prj_004_variant2function/coloc/ghrc38_reference/ukbb_all_chrs_grch38_maf0.01_30000_random_unrelated_white_british"
  
  conditional.datasets=list()
  max.loci=1

# Perform cojo  
  for(i in 1:length(datasets)){
    tmp=cojo.ht(D=datasets[[i]],p.tresh = 1e-4, bfile=bfile)
    conditional.datasets[[i]]=tmp
    names(conditional.datasets)[i]=names(datasets)[i]
    max.loci=max(max.loci,nrow(tmp$ind.snps))
  }
  cat("\nSecondary associations signals identified with COJO\n")
  
############################################################## To delete (?)
  saveRDS(conditional.datasets, file=paste0(opt$output, "/temporary/locus_", locus, "_conditional.datasets.RData"))
  #conditional.datasets <- readRDS(file=paste0(opt$output, "/temporary/locus_", locus, "_conditional.datasets.RData"))
##############################################################
  
  
# Plot of all indipendent associations for each trait
  pdf(paste0(opt$output, "/plots/locus_", locus, "_conditioned_loci.pdf"), height=3.5*max.loci, width=10)
  for(i in 1:length(conditional.datasets)){
    p4 <- plot.cojo.ht(conditional.datasets[[i]]) +
      plot_annotation(paste("Locus", locus, names(conditional.datasets)[i]))
    print(p4)
  }
  dev.off()



# Only if there are multiple traits
  if(length(conditional.datasets)>1){
    
# Identify all pairwise combination of tested traits 
    pairwise.list=t(combn(names(conditional.datasets),2))
# Prepare final locus.table for coloc
    final.locus.table.tmp <- coloc.prep.table(pairwise.list, conditional.datasets, loci.table.tmp)

    
### Run colocalisation and store results
    final.colocs=c()
    for(i in 1:nrow(pairwise.list)){
      coloc.res=colo.cojo.ht(
        conditional.dataset1=conditional.datasets[[pairwise.list[i,1]]],
        conditional.dataset2 = conditional.datasets[[pairwise.list[i,2]]],
        p.threshold.cond = 1e-6,
        p.threshold.orig = 5e-8 
      )
        
      if(!is.null(coloc.res)){
          coloc.res$t1=pairwise.list[i,1]
          coloc.res$t2=pairwise.list[i,2]
          final.colocs=rbind(final.colocs,coloc.res)
        }
    }

    
# Identify traits colocalising (PP.H4 >= 0.9)
    final.colocs.H4=final.colocs[round(final.colocs$PP.H4.abf,digits = 2)>=0.90,]
    final.colocs.H4=as.data.frame(final.colocs.H4)

    colocalization.table.all=c()
    colocalization.table.H4=c()
#   k=1 (needed?)
        
# If any traits colocalisation is present 
    if(nrow(final.colocs.H4)>0){
      
# Define groups based on coloc results 
      a.graph=graph_from_data_frame(final.colocs.H4[,c("hit1","hit2")],directed=F)
      groups=components(a.graph)
      groups=data.frame(snp=names(groups$membership),group=groups$membership)
      final.colocs.H4$g1=groups$group[match(final.colocs.H4$hit1,groups$snp)]
      
      final.locus.table.tmp <- coloc.subgrouping(final.colocs.H4, final.locus.table.tmp)

      
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
    } else {
      idx=which(is.na(final.locus.table.tmp$sub_locus))
      pri=1
      final.locus.table.tmp$sub_locus[idx]=pri:(pri+length(idx)-1)
      final.locus.table.tmp=as.data.frame(final.locus.table.tmp)[,col.order]
    } 
    
    if(nrow(final.colocs)>0){
        final.colocs$locus=locus
        colocalization.table.all=rbind(colocalization.table.all,final.colocs)
    }
    if(nrow(final.colocs.H4)>0){
        colocalization.table.H4=rbind(colocalization.table.H4,final.colocs.H4)
    }

### Summary coloc tables
    write.table(colocalization.table.all,
      file=paste0(opt$output, "/results/locus_", locus, "_colocalization.table.all.tsv"),
      row.names=F,quote=F,sep="\t")
    
    write.table(colocalization.table.H4,
      file=paste0(opt$output, "/results/locus_", locus, "_colocalization.table.H4.tsv"),
      row.names=F,quote=F,sep="\t")
    
  }else{
    final.locus.table.tmp=conditional.datasets[[1]]$ind.snps
    final.locus.table.tmp$start=locus.info$start
    final.locus.table.tmp$end=locus.info$end
    final.locus.table.tmp$pan.locus=locus
    final.locus.table.tmp$sub_locus=1
    final.locus.table.tmp$freq_geno=NA
    final.locus.table.tmp$bJ=NA
    final.locus.table.tmp$bJ_se=NA
    final.locus.table.tmp$pJ=NA
    final.locus.table.tmp$LD_r=NA
          
    alleles=unlist(mappa.loc[mappa.loc$SNP==final.locus.table.tmp$SNP,c("A1","A2")])
    final.locus.table.tmp$othA=alleles[!(alleles%in%final.locus.table.tmp$refA)]
    final.locus.table.tmp$trait=names(datasets)[1]
    final.locus.table.tmp=as.data.frame(final.locus.table.tmp)
    final.locus.table.tmp=final.locus.table.tmp[,col.order]
  }
  
  final.locus.table=rbind(final.locus.table,final.locus.table.tmp)
  print(final.locus.table)
  
  write.table(final.locus.table, file=paste0(opt$output, "/results/locus_", locus, "_final_locus_table.tsv"),
    row.names=F,quote=F,sep="\t")
}

cat("\n**** DONE!! ****\n")
