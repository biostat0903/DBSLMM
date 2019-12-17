library(data.table)
library(optparse)

## Parameter setting
args_list <- list(
  make_option("--dbslmm", type = "character", default = NULL,
              help = "INPUT:  the result of dbslmm", metavar = "character"),
  make_option("--external", type = "character", default = NULL,
              help = "INPUT:  the external summary statistics", metavar = "character"),
  make_option("--ref", type = "character", default = NULL,
              help = "INPUT:  the perfix of reference panel", metavar = "character"),
  make_option("--valid", type = "character", default = NULL,
              help = "INPUT:  the prefix of valid", metavar = "character"),
  make_option("--block", type = "character", default = NULL,
              help = "INPUT:  the block information (Berisa and Pickrell 2015)", metavar = "character"),
  make_option("--chr", type = "character", default = NULL,
              help = "INPUT:  the chromosome number", metavar = "character"),
  make_option("--mafMax", type = "character", default = NULL,
              help = "INPUT:  the maximium of the difference between reference panel and external data", metavar = "character"),
  make_option("--outPath", type="character", default=NULL,
              help = "OUTPUT: the output path", metavar="character")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

## output the options
cat("Acccept Options: \n")
cat("--dbslmm:   ", opt$dbslmm, "\n")
cat("--external: ", opt$external, "\n")
cat("--ref:      ", opt$ref, "\n")
cat("--block:    ", opt$block, "\n")
cat("--chr:      ", opt$chr, "\n")
mafMax <- ifelse(is.null(opt$mafMax), 0.2, opt$mafMax)
cat("--mafMax:   ", mafMax, "\n")
cat("--outPath:  ", opt$outPath, "\n")

## check the options
if (!file.exists(opt$external)){
  cat(paste0("ERROR: ", opt$external, " does not exist!\n"))
  q()
}
if(opt$chr != "all"){
  dbslmm_chr <- paste0(opt$dbslmm, opt$chr, ".assoc.dbslmm.txt")
  if(!file.exists(dbslmm_chr)){
    cat(paste0("ERROR: ", dbslmm_chr, " does not exist!\n"))
    q()
  }
  ref_chr1 <- paste0(opt$ref, opt$chr, ".bed")
  ref_chr2 <- paste0(opt$ref, opt$chr, ".bim")
  ref_chr3 <- paste0(opt$ref, opt$chr, ".fam")
  if(!file.exists(ref_chr1) || !file.exists(ref_chr2) || !file.exists(ref_chr3)){
    cat(paste0("ERROR: ", opt$ref, opt$chr, " does not exist!\n"))
    q()
  }
  block_chr <- paste0(opt$block, opt$chr, ".bed")
  if(!file.exists(block_chr)){
    cat(paste0("ERROR: ", block_chr, " does not exist!\n"))
    q()
  }
} else {
  for (chr in 1: 22){
    dbslmm_chr <- paste0(opt$dbslmm, chr, ".assoc.dbslmm.txt")
    if(!file.exists(dbslmm_chr)){
      cat(paste0("ERROR: ", dbslmm_chr, " does not exist!\n"))
      q()
    }
  }
  for (chr in 1: 22){
    ref_chr1 <- paste0(opt$ref, chr, ".bed")
    ref_chr2 <- paste0(opt$ref, chr, ".bim")
    ref_chr3 <- paste0(opt$ref, chr, ".fam")
    if(!file.exists(ref_chr1) || !file.exists(ref_chr2) || !file.exists(ref_chr3)){
     cat(paste0("ERROR: ", opt$ref, chr, " does not exist!\n"))
      q()
    }
  }
  for (chr in 1: 22){
    block_chr <- paste0(opt$block, chr, ".bed")
    if(!file.exists(block_chr)){
      cat(paste0("ERROR: ", block_chr, " does not exist!\n"))
      q()
    }
  }
}

## valid function
if (opt$chr != "all"){
  valid_cmd <- paste0(opt$valid, " -d ", dbslmm_chr, " -s ", opt$external, 
                      " -r ", opt$ref, opt$chr, " -b ", opt$block, opt$chr, ".bed -mafMax ", mafMax, 
                      " -r2 ", opt$outPath, "/r2_chr", opt$chr)
  system(valid_cmd)
  r2_dat <- read.table(paste0(opt$outPath, "/r2_chr", opt$chr, ".txt"))
  r <- sum(r2_dat[, 1]) / sqrt(sum(r2_dat[, 2]))
  cat ("R square between internal and external data on Chromosome ", opt$chr, " is ", r^2, "\n")
} else {
  r2_dat_tot <- vector()
  for (chr in 1:22){
    valid_cmd <- paste0(opt$valid, " -d ", opt$dbslmm, chr, ".assoc.sbslmm.txt -s ", opt$external, 
                        " -r ", opt$ref, chr, " -b ", opt$block, opt$chr, ".bed -mafMax ", mafMax, 
                        " -r2 ", opt$outPath, "/r2_chr", opt$chr)
    system(valid_cmd)
    r2_dat <- read.table(paste0(opt$outPath, "/r2_chr", opt$chr, ".txt"))
    r2_dat_tot <- rbind(r2_dat_tot, r2_dat)
    r <- sum(r2_dat_tot[, 1]) / sqrt(sum(r2_dat_tot[, 2]))
    cat ("R square between internal and external data on whole genome is ", r^2, "\n")
  }
}

