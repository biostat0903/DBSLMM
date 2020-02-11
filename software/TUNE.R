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

library(data.table)
library(Metrics)
require(scoring)
require(pROC)
library(optparse)

## Parameter setting
args_list <- list(
  make_option("--phenoPred", type = "character", default = NULL,
              help = "INPUT: the predicted phenotype (PLINK output)", metavar = "character"),
  make_option("--phenoVal", type = "character", default = NULL,
              help = "INPUT: the validation phenotype file and column number", metavar = "character"),
  make_option("--index", type="character", default="r2",
              help="INPUT: four different indexes (MSE and R2 for continuous trait, Brier score and AUC for binary trait)", 
              metavar="character"), 
  make_option("--pvRange", type = "character", default = NULL,
              help = "INPUT: the range of p-value threshold", metavar = "character"),
  make_option("--ldRange", type = "character", default = NULL,
              help = "INPUT: the range of LD threshold", metavar = "character"),
  make_option("--cov", type = "character", default = NULL,
              help = "INPUT: covariate", metavar = "character"),
  make_option("--chr", type = "character", default = NULL,
              help = "INPUT: chromosome", metavar = "character")
)

opt_parser <- OptionParser(option_list=args_list)
opt <- parse_args(opt_parser)

## output the options
cat("Acccept Options: \n")
cat("--phenoPred: ", opt$phenoPred, "\n")
phenoVal <- unlist(strsplit(opt$phenoVal, ","))[1]
phenoNum <- unlist(strsplit(opt$phenoVal, ","))[2]
cat("--phenoVal:  ", phenoVal, "\n")
cat("--phenoNum:  ", phenoNum, "\n")
cat("--index:     ", opt$index, "\n")
if (is.null(opt$pvRange) == F){
  cat("--pvRange:   ", opt$pvRange, "\n")
  pv_vec <- unlist(strsplit(opt$pvRange, ","))
}
if (is.null(opt$ldRange) == F){
  cat("--ldRange:   ", opt$ldRange, "\n")
  ld_vec <- unlist(strsplit(opt$ldRange, ","))
}
if (length(pv_vec) == 1 & length(ld_vec) == 1){
  cat("ERROR: no tuning parameter!\n")
  q()
}

if(is.null(opt$chr) == F){
  if (as.numeric(opt$chr) %in% c(1:22)){
    cat("--chr:       ", opt$chr, "\n")
    cat(paste0("WARNNING: Only chromosome ", opt$chr, "is included!\n"))
  }
} else {
  cat(paste0("Validate for whole genome!\n"))
}
if (is.null(opt$cov) == F){
  cat("--cov:       ", opt$cov, "\n")
} else{
  cat("Please check: phenotypes are removed covariates!\n")
}

## check file
if (opt$index != "r2" & opt$index != "mse" & opt$index != "bs" & opt$index == "auc"){
  cat("ERROR: The index is wrong!\n")
  q()
}

pheno <- data.frame(fread(phenoVal))[, as.numeric(phenoNum)]
ld_num <- length(ld_vec)
pv_num <- length(pv_vec)
r2_res <- vector("numeric", pv_num*ld_num)
mse_res <- vector("numeric", pv_num*ld_num)
bs_res <- vector("numeric", pv_num*ld_num)
auc_res <- vector("numeric", pv_num*ld_num)
for (pp in 1: pv_num){
  for (rr in 1: ld_num){
    if (is.null(opt$chr) == F){
      pheno_chr1_str <- paste0(opt$phenoPred, opt$chr, "_pv", pv_vec[pp],
                               "_r", ld_vec[rr], ".profile")
      geno_pheno <- data.frame(fread(pheno_chr1_str, header = T))[, 6]
    } else {
      pheno_chr1_str <- paste0(opt$phenoPred, 1, "_pv", pv_vec[pp],
                               "_r", ld_vec[rr], ".profile")
      geno_pheno <- data.frame(fread(pheno_chr1_str, header = T))[, 6]
      for (chr in 2: 22){
        pheno_chr_str <- paste0(opt$phenoPred, chr, "_pv", pv_vec[pp],
                                "_r", ld_vec[rr], ".profile")
        pred_chr <- data.frame(fread(pheno_chr_str, header = T))[, 6]
        geno_pheno <- geno_pheno + pred_chr
      }
    }
    if (opt$index == "auc" || opt$index == "bs"){
      cov <- as.matrix(fread(opt$cov))
      coef_logit <- coef(glm(pheno~cov, family = "binomial"))
      pred_pheno <- geno_pheno + cbind(1, cov) %*% coef_logit
      bs_res[c((pp-1)*ld_num+rr)] <- mean(brierscore(pheno ~ pred_pheno))
      auc_res[c((pp-1)*ld_num+rr)] <- as.numeric(auc(roc(pheno, pred_pheno, levels = c(0, 1))))
    }
    if ( (opt$index == "r2" || opt$index == "mse")){
      if (is.null(opt$cov) == F){
        cov <- as.matrix(fread(opt$cov))
        coef_lm <- coef(lm(pheno~cov))
        geno_pheno <- geno_pheno + cbind(1, cov) %*% coef_lm
      } 
      r2_res[c((pp-1)*ld_num+rr)] <- cor(geno_pheno, pheno)^2
      mse_res[c((pp-1)*ld_num+rr)] <- mse(geno_pheno, pheno)
    }
  }
}

## selection
if (opt$index == "r2"){
  out_file <- data.frame(pv = as.numeric(rep(pv_vec, each=ld_num)),
                         r2 = as.numeric(ld_vec),
                         R2 = r2_res)
  out_best <- out_file[which.max(out_file$R2), c(1, 2, 3)]
  cat("Using", opt$index, "the best combination: p-threshold: ", out_best[1, 1],
      " LD threshold: ", out_best[1, 2], ".\n")
}
if (opt$index == "mse"){
  out_file <- data.frame(pv = as.numeric(rep(pv_vec, each=ld_num)),
                         r2 = as.numeric(ld_vec),
                         mse = mse_res)
  out_best <- out_file[which.min(out_file$mse), ]
  cat("Using", opt$index, "the best combination: p-threshold: ", out_best[1, 1],
      " LD threshold: ", out_best[1, 2], ".\n")
}
if (opt$index == "auc"){
  out_file <- data.frame(pv = as.numeric(rep(pv_vec, each=ld_num)),
                         r2 = as.numeric(ld_vec),
                         auc = auc_res)
  out_best <- out_file[which.max(out_file$auc), ]
  cat("Using", opt$index, "the best combination: p-threshold: ", out_best[1, 1],
      " LD threshold: ", out_best[1, 2], ".\n")
}
if (opt$index == "bs"){
  out_file <- data.frame(pv = as.numeric(rep(pv_vec, each=ld_num)),
                         r2 = as.numeric(ld_vec),
                         bs = bs_res)
  out_best <- out_file[which.min(out_file$bs), ]
  cat("Using", opt$index, "the best combination: p-threshold: ", out_best[1, 1],
      " LD threshold: ", out_best[1, 2], ".\n")
}

## output
write.table(out_best[1], paste0(opt$phenoPred, "_pbest.", opt$index), 
            row.names = F, col.names = F, quote = F)
write.table(out_best[2], paste0(opt$phenoPred, "_rbest.", opt$index), 
            row.names = F, col.names = F, quote = F)
