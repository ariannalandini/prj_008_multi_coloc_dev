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
  if(!is.null(trait.res)){
    trait.res = as.data.frame(trait.res, stringsAsFactors = FALSE)
    trait.res = trait.res[, -(which(names(trait.res) == chr.label))]
    names(trait.res)[1:3] = c("chr", "start", "end")
    rownames(trait.res) <- NULL
  }
  return(trait.res)
}


### locus.lister
locus.lister <- function(
  my_path_loci=NULL,
  out="./panlocus_table.tsv"){
  
#  require(data.table)
#  require(GenomicRanges)
#  require(dplyr)

  loci_list <- list.files(path=my_path_loci, pattern="*_loci.tsv", full.names=T)
  all_loci <- as.data.frame(rbindlist(lapply(loci_list, function(x) fread(x, data.table=F)), fill=TRUE))
  cat("\nIdentify overlapping loci across traits\n")
  
  pan_loci <- reduce(GRanges(seqnames = all_loci$chr, IRanges(as.numeric(all_loci$start), as.numeric(all_loci$end)))) 
  pan_loci_non_reduced <- GRanges(seqnames = all_loci$chr, IRanges(as.numeric(all_loci$start), as.numeric(all_loci$end)))
  
  overlapping <- findOverlaps(pan_loci_non_reduced, pan_loci) # find overlaps between all the loci and use it as an index of the unique non overlapping loci
  all_loci$pan_locus <- rep(0, nrow(all_loci)) # allocate the column 
  all_loci[overlapping@from,]$pan_locus <- overlapping@to  # assigning the number as index of which macro loci is overlapping 
  all_loci$pan_locus_name <- rep(0, nrow(all_loci))
  
  # assign a name refereed to the position of each pan_locus
  for(k in 1:length(unique(all_loci$pan_locus))){
    all_loci[which(all_loci$pan_locus==k), ]$pan_locus_name <- paste0(all_loci[which(all_loci$pan_locus==k), ]$chr,'_',  min(all_loci[which(all_loci$pan_locus==k), ]$start), '_',  max(all_loci[which(all_loci$pan_locus==k), ]$end))
  }
  rownames(all_loci) <- NULL
  cat(paste0("\nTotal of ", length(unique(all_loci$pan_locus)), " pan loci identified\n"))
  fwrite(all_loci %>% arrange(pan_locus, trait), out, sep="\t", quote=F, na=NA)
  #  return(all_loci)
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
