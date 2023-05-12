#### Functions for loci identification


### locus.breaker
locus.breaker <- function(
  res,
  p.sig = 5e-08,
  p.limit = 1e-05,
  hole.size = 250000,
  p.label = "P",
  chr.label = "CHR",
  pos.label = "BP"){
  
  res = res[order(as.numeric(res[, chr.label]), as.numeric(res[,pos.label])), ]
  res = res[which(res[, p.label] < p.limit), ]
  trait.res = c()
  for (j in 1:22) {
    res.chr = res[which(res[, chr.label] == j), ]
    if (nrow(res.chr) > 1) {
      holes = res.chr[, pos.label][-1] - res.chr[, pos.label][-length(res.chr[,pos.label])]
      gaps = which(holes > hole.size)
      if (length(gaps) > 0) {
        for (k in 1:(length(gaps) + 1)) {
          if (k == 1) {
            res.loc = res.chr[1:(gaps[k]), ]
          }
          else if (k == (length(gaps) + 1)) {
            res.loc = res.chr[(gaps[k - 1] + 1):nrow(res.chr), 
            ]
          } else {
            res.loc = res.chr[(gaps[k - 1] + 1):(gaps[k]), 
            ]
          }
          if (min(res.loc[, p.label]) < p.sig) {
            start.pos = min(res.loc[, pos.label], na.rm = T)
            end.pos = max(res.loc[, pos.label], na.rm = T)
            chr = j
            best.snp = res.loc[which.min(res.loc[, p.label]), 
            ]
            line.res = c(chr, start.pos, end.pos, unlist(best.snp))
            trait.res = rbind(trait.res, line.res)
          }
        }
      } else {
        res.loc = res.chr
        if (min(res.loc[, p.label]) < p.sig) {
          start.pos = min(res.loc[, pos.label], na.rm = T)
          end.pos = max(res.loc[, pos.label], na.rm = T)
          chr = j
          best.snp = res.loc[which.min(res.loc[, p.label]), 
          ]
          line.res = c(chr, start.pos, end.pos, unlist(best.snp))
          trait.res = rbind(trait.res, line.res)
        }
      }
    }
    else if (nrow(res.chr) == 1) {
      res.loc = res.chr
      if (min(res.loc[, p.label]) < p.sig) {
        start.pos = min(res.loc[, pos.label], na.rm = T)
        end.pos = max(res.loc[, pos.label], na.rm = T)
        chr = j
        best.snp = res.loc[which.min(res.loc[, p.label]), 
        ]
        line.res = c(chr, start.pos, end.pos, unlist(best.snp))
        trait.res = rbind(trait.res, line.res)
      }
    }
  }
  trait.res = as.data.frame(trait.res, stringsAsFactors = FALSE)
  trait.res = trait.res[, -(which(names(trait.res) == chr.label))]
  names(trait.res)[1:3] = c("chr", "start", "end")
  rownames(trait.res) <- NULL
  return(trait.res)
}


### locus.lister
locus.lister <- function(
  my_path_gwas=NULL,
  traits_pref="",
  traits_suf=".regenie.gz",
  p.sig = 5e-08,
  p.limit = 1e-05,
  hole.size = 250000,
  my_path_loci=NULL,
  SNP="ID",
  CHR="CHROM",
  BP="GENPOS",
  A1="ALLELE1",
  A2="ALLELE0",
  BETA="BETA",
  SE="SE",
  P="P",
  out="./panlocus_table.tsv"){
  
  require(data.table)
  require(GenomicRanges)
  require(dplyr)
#  source("/group/pirastu/prj_004_variant2function/scripts/locusbreaker_AL.R")
  
  if(is.null(my_path_loci)){
    gwas <- list.files(path=my_path_gwas, pattern="*")
    cat(paste0("\n", length(gwas), " files found in ",my_path_gwas, " folder\n"))
    SNPs <- list()
    list_of_files <- list()
    # load the sumstats and put them into a list
    for(i in 1:length(gwas)){
      list_of_files[[i]] <- fread(paste0(my_path_gwas, "/", gwas)[[i]], data.table = F)
      cat(paste0("\n", gwas[[i]], " loaded\n"))
      SNPs[[i]] <- list_of_files[[i]][[SNP]]
    }
    cat("\nAll GWAS summary statistics loaded!\n")
    shared_SNPs <- Reduce(intersect, SNPs) # compute the SNPs that are present in all the summary stats
    
    if(is.null(gwas_traits)){gwas_traits <- gsub(paste0(traits_pref, "(.*)", traits_suf), "\\1", gwas)}
    names(list_of_files) <- unlist(gwas_traits)
    
    loci <- list()
    for(i in 1:length(gwas)){
      to_take <- gwas_traits[[i]]
      sstats <- list_of_files[[to_take]]
      colnames(sstats) <- toupper(colnames(sstats))   
      sstats <- sstats %>% select(one_of(c(SNP, CHR, BP, A2, A1, BETA, SE, P)))
      sstats <- sstats[which(sstats[[SNP]] %in% shared_SNPs),] # select only the SNPs that are shared among all the traits     
      loci[[i]] <- locus.breaker(sstats, p.label=P, chr.label=CHR, pos.label=BP, p.sig=p.sig, p.limit=p.limit, hole.size=hole.size)
      loci[[i]]$trait <- rep(to_take, nrow(loci[[i]]))
      loci[[i]]$path <- rep(paste0(my_path_gwas, "/", gwas[[i]]), nrow(loci[[i]]))
      cat(paste0("\nTrait-specific loci identified for ", to_take, "\n"))
    }
    cat("\nTrait-specific loci identified for all traits!\n")
  } else{
# If you have precomputed trait-specific loci, generated with locus.breaker function
    loci_list <- list.files(path=my_path_loci, pattern = "*_loci.tsv", full.names=T)
    loci <- lapply(loci_list, function(x) fread(x, data.table=F))
  }

  all_loci <- do.call(rbind, loci) # create the list of all the loci
  cat("\nIdentify overlapping loci across traits\n")
  
  pan_loci <- reduce(GRanges(seqnames = all_loci$chr, IRanges(as.numeric(all_loci$start), as.numeric(all_loci$end)))) 
  pan_loci_non_reduced <- GRanges(seqnames = all_loci$chr, IRanges(as.numeric(all_loci$start), as.numeric(all_loci$end)))
  
  overlapping <- findOverlaps(pan_loci_non_reduced, pan_loci) #find overlaps between all the loci and use it as an index of the unique non overlapping loci
  all_loci$pan_locus <- rep(0, nrow(all_loci)) #allocate the column 
  all_loci[overlapping@from,]$pan_locus <- overlapping@to  # assigning the number as index of which macro loci is overlapping 
  all_loci$pan_locus_name <- rep(0, nrow(all_loci))
  
  # assign a name refereed to the position of each pan_locus
  for(k in 1:length(unique(all_loci$pan_locus))){
    all_loci[which(all_loci$pan_locus==k), ]$pan_locus_name <- paste0(all_loci[which(all_loci$pan_locus==k), ]$chr,'_',  min(all_loci[which(all_loci$pan_locus==k), ]$start), '_',  max(all_loci[which(all_loci$pan_locus==k), ]$end))
  }
  rownames(all_loci) <- NULL
  cat(paste0("\nTotal of ", length(unique(all_loci$pan_locus)), " pan loci identified\n"))
  fwrite(all_loci, out, sep="\t", quote=F, na=NA)
  #  return(all_loci)
}
