#! /usr/bin/env Rscript

########################################################################
# Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)            #
# Copyright (C) 2019  Sheng Yang and Xiang Zhou                        #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
########################################################################

library(bigreadr)
library(optparse)

## Parameter setting
args_list <- list(
  make_option("--summary", type = "character", default = NULL,
              help = "INPUT: the summary statistics (gemma output format)", metavar = "character"),
  make_option("--plink", type = "character", default = NULL,
              help = "INPUT: the perfix of Plink software", metavar = "character"),
  make_option("--dbslmm", type = "character", default = NULL,
              help = "INPUT: the perfix of dbslmm software", metavar = "character"),
  make_option("--ref", type = "character", default = NULL,
              help = "INPUT: the perfix of reference panel", metavar = "character"),
  make_option("--block", type = "character", default = NULL,
              help = "INPUT: the block information (Berisa and Pickrell 2015)", metavar = "character"),
  make_option("--outPath", type="character", default=NULL,
              help="INPUT: the output path", metavar="character"),
  make_option("--n", type = "character", default = NULL,
              help = "INPUT: the sample size of summary data", metavar = "character"),
  make_option("--nsnp", type = "character", default = NULL,
              help = "INPUT: the number of SNPs in whole genome", metavar = "character"),
  make_option("--h2", type = "character", default = NULL,
              help = "INPUT: the heritability of trait", metavar = "character"),
  # make_option("--h2f", type = "character", default = NULL,
  #             help = "INPUT: the heritability of trait", metavar = "character"),
  make_option("--type", type="character", default="d",
              help="INPUT: type of DBSLMM (default: default version)", 
              metavar="character"),
  make_option("--r2", type = "character", default = "0.1",
              help = "INPUT: the cutoff of SNPs clumping (default: 0.1)", 
              metavar = "character"),
  make_option("--pv", type = "character", default = "1e-6",
              help = "INPUT: the cutoff of SNPs pruning (default: 1e-6)", 
              metavar = "character"),
  make_option("--mafMax", type = "character", default = "0.2",
              help = "INPUT: the maximium of the difference between reference panel and summary data (default:0.2)", 
              metavar = "character"),
  make_option("--thread", type = "character", default = "5",
              help = "INPUT: the number of threads (default: 5)", 
              metavar = "character")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

## output the options
cat("Acccept Options: \n")
cat("--summary:  ", opt$summary, "\n")
cat("--plink:    ", opt$plink, "\n")
cat("--dbslmm:   ", opt$dbslmm, "\n")
cat("--ref:      ", opt$ref, "\n")
cat("--block:    ", opt$block, "\n")
cat("--outPath:  ", opt$outPath, "\n")
cat("--n:        ", opt$n, "\n")
cat("--nsnp:     ", opt$nsnp, "\n")
cat("--h2:       ", opt$h2, "\n")
cat("--type:     ", opt$type, "\n")
cat("--r2:       ", opt$r2, "\n")
cat("--p:        ", opt$pv, "\n")
cat("--mafMax:   ", opt$mafMax, "\n")
cat("--thread:   ", opt$thread, "\n")

## check the options
if (!file.exists(opt$summary)){
  cat(paste0("ERROR: ", opt$summary, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$outPath)){
  cat(paste0("ERROR: ", opt$outPath, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$plink)){
  cat(paste0("ERROR: ", opt$plink, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$dbslmm)){
  cat(paste0("ERROR: ", opt$dbslmm, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$plink)){
  cat(paste0("ERROR: ", opt$plink, " does not exist! Please check!\n"))
  q()
}
if (!file.exists(opt$block)){
  cat(paste0("ERROR: ", opt$block, " does not exist! Please check!\n"))
  q()
}

start <- proc.time()
## gemma format to plink format
s <- strsplit(opt$summary, "\\/")[[1]]
s <- s[length(s)]
prefix_file <- strsplit(s, "\\.")[[1]]
len_prefix_file <- length(prefix_file)
prefix_file <- paste(prefix_file[-c((len_prefix_file-1):len_prefix_file)], collapse = ".")

if (opt$type != "d"){
  prefix_file <- paste0(prefix_file, "_pv", opt$pv, "_r", opt$r2)
}

summstats <- fread2(opt$summ)
t <- summstats[, 9]/summstats[, 10]
p_val <- ifelse(t < 0, pnorm(t), pnorm(t, lower.tail = F))*2
p_val_z <- ifelse(p_val == 0, min(p_val[-which(p_val==0)]), p_val)
plink_summstats <- data.frame(SNP = summstats[, 2], P = p_val_z)
write.table(plink_summstats, file = paste0(opt$outPath, "plink_", prefix_file, ".txt"),
            row.names = F, quote = F)

## large effect clumping
l_clumping_cmd <- paste0(opt$plink, " --bfile ", opt$ref, 
                         " --silent  --clump ", opt$outPath, "plink_", prefix_file, ".txt",
                         " --clump-r2 ", opt$r2, " --clump-p1 ", opt$pv,
                         " --clump-kb 10000", 
                         " --out ", opt$outPath, "/l_", prefix_file)
system(l_clumping_cmd)

if (!file.exists(paste0(opt$outPath, "l_", prefix_file, ".clumped"))){
  cat ("Using", opt$pv, ", no significant SNPs. The most significant SNPs are regarded as large effect SNPs.\n")
  sort_cmd <- paste0("sort -k 11 -n ", opt$summary, " | head -n 1 > ",
                     opt$outPath, "/l_snp_", prefix_file, ".txt")
  system(sort_cmd)
  lsnp_p <- read.table(paste0(opt$outPath, "/l_snp_", prefix_file, ".txt"))[1, 2]
} else {
  tosnp_cmd <- paste0("awk '{print $3}' ", opt$outPath, "l_", prefix_file, ".clumped", " > ",
                      opt$outPath, "/l_snp_", prefix_file, ".txt")
  system(tosnp_cmd)
  lsnp <- read.table(paste0(opt$outPath, "/l_snp_", prefix_file, ".txt"), header = T,
                     stringsAsFactors = F)[, 1]
  cat ("Using", opt$pv, ",", length(lsnp), "SNPs are regarded as large effect.\n")
  lsnp_p <- paste(lsnp, collapse="|")
}
l_cmd <- paste0("grep -w -E -a '", lsnp_p, "' ", opt$summary, " > ",
                opt$outPath, "/l_", prefix_file, ".txt")
system(l_cmd)

## small effect clumping
s_cmd1 <- paste0("grep -w -E -a -v '", lsnp_p, "' ", opt$outPath, "plink_", prefix_file, ".txt", " > ",
                 opt$outPath, "/s_plink_", prefix_file, ".txt")
s_clumping_cmd <- paste0(opt$plink, " --bfile ", opt$ref, 
                         " --silent --clump ", opt$outPath, "s_plink_", prefix_file, ".txt",
                         " --clump-r2 ", 0.9, " --clump-p1 ", 0.9,
                         " --clump-kb 10000", 
                         " --out ", opt$outPath, "/s_", prefix_file)
system(s_cmd1)
system(s_clumping_cmd)
tosnp_cmd <- paste0("awk '{print $3}' ", opt$outPath, "s_", prefix_file, ".clumped", " > ",
                    opt$outPath, "/s_snp_", prefix_file, ".txt")
system(tosnp_cmd)
ssnp <- read.table(paste0(opt$outPath, "/s_snp_", prefix_file, ".txt"), header = T,
                   stringsAsFactors = F)[, 1]
s_summstats <- summstats[summstats[, 2] %in% ssnp, ]
write.table(s_summstats, file = paste0(opt$outPath, "/s_", prefix_file, ".txt"), 
            col.names = F, row.names = F, quote = F, sep = "\t")
cat(paste0(dim(s_summstats)[1], " SNPs are regarded as small effect.\n"))
end <- proc.time()
cat("Clumping time: ", end[3]-start[3], "s.\n")

## dbslmm
# h2d <- as.numeric(opt$h2f) * as.numeric(opt$h2)
system(paste0(opt$dbslmm,
              " -s ",      opt$outPath, "s_", prefix_file, ".txt",
              " -l ",      opt$outPath, "l_", prefix_file, ".txt",
              " -r ",      opt$ref,
              " -n ",      opt$n,
              " -mafMax ", opt$mafMax,
              " -nsnp ",   opt$nsnp,
              " -b ",      opt$block,
              " -h ",      opt$h2,
              " -t ",      opt$thread,
              " -eff ",    opt$outPath, prefix_file, ".dbslmm"))
system(paste0("rm ", opt$outPath, "plink_", prefix_file, ".txt"))
system(paste0("rm ", opt$outPath, "l_", prefix_file, "*"))
system(paste0("rm ", opt$outPath, "s_", prefix_file, "*"))
system(paste0("rm ", opt$outPath, "s_plink_", prefix_file, ".txt"))
system(paste0("rm ", opt$outPath, "l_snp_", prefix_file, ".txt"))
system(paste0("rm ", opt$outPath, "s_snp_", prefix_file, ".txt"))
