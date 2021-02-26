#!/bin/bash

while getopts "D:p:B:s:H:m:n:G:P:l:T:c:i:t:o:" opt; do
  case $opt in
    D) software_path="$OPTARG"
    ;;
    p) plink="$OPTARG"
    ;;
    B) block_prefix="$OPTARG"
    ;;
    s) summary_file_prefix="$OPTARG"
    ;;
    H) herit="$OPTARG"
    ;;
    m) nsnp="$OPTARG"
    ;;
    n) n="$OPTARG"
    ;;
    G) val_geno="$OPTARG"
    ;;
    P) val_pheno="$OPTARG"
    ;;
    l) col="$OPTARG"
    ;;
    T) type="$OPTARG"
    ;;
    c) cov="$OPTARG"
    ;;
    i) index="$OPTARG"
    ;;
    t) thread="$OPTARG"
    ;;
    o) outpath="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

printf "\033[33mArgument software_path is %s  \033[0m\n" "$software_path"
printf "\033[33mArgument plink is %s  \033[0m\n" "$plink"
printf "\033[33mArgument block_prefix is %s  \033[0m\n" "$block_prefix"
printf "\033[33mArgument summary_file_prefix is %s  \033[0m\n" "$summary_file_prefix"
printf "\033[33mArgument herit is %s  \033[0m\n" "$herit"
printf "\033[33mArgument val_geno is %s  \033[0m\n" "$val_geno"
printf "\033[33mArgument valid_pheno is %s  \033[0m\n" "$val_pheno"
printf "\033[33mArgument col is %s  \033[0m\n" "$col"
printf "\033[33mArgument type is %s  \033[0m\n" "$type"
if [ -z "${cov}" ]; then 
cov='N'
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
else 
printf "\033[33mArgument cov is %s  \033[0m\n" "$cov"
fi
printf "\033[33mArgument index is %s  \033[0m\n" "$index"
printf "\033[33mArgument thread is %s  \033[0m\n" "$thread"
printf "\033[33mArgument outpath is %s  \033[0m\n" "$outpath"

# DBSLMM software
DBSLMM=${software_path}DBSLMM/software/DBSLMM.R
TUNE=${software_path}DBSLMM/software/TUNE.R
dbslmm=${software_path}/DBSLMM/software/dbslmm

# split to chromosome
Rscript ${SPLIT} --summary ${summary_file_prefix}.assoc.txt

# DBSLMM: tuning version
if [[ "$type" == "t" ]]
then
for chr in `seq 1 22`
do

BLOCK=${block_prefix}${chr}
summchr=${summary_file_prefix}_chr${chr}
Rscript ${DBSLMM} --summary ${summchr}.assoc.txt --outPath ${outpath} --plink ${plink}\
                  --dbslmm ${dbslmm} --ref ${val_geno} --n ${n} --type ${type} --nsnp ${nsnp} --block ${BLOCK}.bed\
                  --h2 ${herit} --h2f 0.8,1,1.2 --thread ${thread}
for h2f in 0.8 1 1.2
do
if [ -f "${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.badsnps" ];then
rm ${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.badsnps
fi
summchr_prefix=`echo ${summchr##*/}`
${plink}  --silent --bfile ${val_geno} --score ${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.txt 1 2 4 sum\
          --out ${outpath}${summchr_prefix}_h2f${h2f}
rm ${outpath}${summchr_prefix}_h2f${h2f}.log
if [ -f "${outpath}${summchr_prefix}_h2f${h2f}.nopred" ];then
rm ${outpath}${summchr_prefix}_h2f${h2f}.nopred
fi

done
done

summchr_prefix2=`echo ${summchr_prefix%_*}`
if [[ "${cov}" == "N" ]]
then 
Rscript ${TUNE} --phenoPred ${outpath}${summchr_prefix2} --phenoVal ${val_pheno},${col} \
       --h2Range 0.8,1,1.2 --index ${index}
else 
Rscript ${TUNE} --phenoPred ${outpath}${summchr_prefix2} --phenoVal ${val_pheno},${col} \
       --h2Range 0.8,1,1.2 --index ${index} --cov ${cov}
fi

hbest=`cat ${outpath}${summchr_prefix2}_hbest.${index}`
for chr in `seq 1 22`
do
mv ${outpath}${summchr_prefix2}_chr${chr}_h2f${hbest}.dbslmm.txt ${outpath}${summchr_prefix2}_chr${chr}_best.dbslmm.txt
rm ${outpath}${summchr_prefix2}_chr${chr}_h2f*
done

fi

## DBSLMM default version
if [[ "$type" == "d" ]]
then
for chr in `seq 1 22` 
do
BLOCK=${block_prefix}${chr}
summchr=${summary_file_prefix}${chr}
val_geno=${val_geno}${chr}
Rscript ${DBSLMM} --summary ${summchr}.assoc.txt --outPath ${outpath} --plink ${plink}\
                  --dbslmm ${dbslmm} --ref ${val_geno} --n ${n} --nsnp ${nsnp} --block ${BLOCK}.bed\
                  --h2 ${herit} --thread ${thread}
summchr_prefix=`echo ${summchr##*/}`
rm ${outpath}${summchr_prefix}.dbslmm.badsnps

done
fi
