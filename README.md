# DBSLMM
Deterministic Bayesian Sparse Linear Mixed Model

## Update log
### v0.21
* fits DBSLMM in different heritabiliies with one clumping
* uses a script to make the usage of DBSLMM-tuning version more easier.
### v0.2 
* validates the flods of heritability: 0.8, 1 and 1.2
* fits the LMM, when the chromosome without any large effect SNPs
* fits the external validation all by R code. 

## Tutorial for DBSLMM (v0.21)
In this version, we use `DBSLMM_script.sh` to make the usage of DBSLMM more easier. Different to past version, v0.21 needs `valg` to be the whole genome plink file.
````bash
PACK_DIR=/your/path/
SPLITCHR=${PACK_DIR}DBSLMM/software/SPLITCHR.R
PLINK=${PACK_DIR}plink_1.9/plink
DBSLMM=${PACK_DIR}DBSLMM/DBSLMM_script.sh

# Split into chromosome
DATADIR=/your/path/example_data/
summ=${DATADIR}all/summary
Rscript ${SPLITCHR} --summary ${summ}.assoc.txt

# Set parameters
valg=${DATADIR}val/valid
valp=${DATADIR}val/valid_pheno.txt
BLOCK=${PACK_DIR}DBSLMM/LDblock/chr
herit=/your/path/herit.log
index=r2
thread=1
outpath=${DATADIR}output/

# DBSLMM tuning version (without covariates)
type=t
sh ${DBSLMM} -D ${PACK_DIR} -p ${PLINK} -B ${BLOCK} -s ${summ} -H ${herit} -n 300 -G ${valg} -P ${valp} -l ${col} -T ${type}  -i ${index} -t ${thread} -o ${outpath}

# DBSLMM determinitic version (without covariates)
type=d
sh ${DBSLMM} -D ${PACK_DIR} -p ${PLINK} -B ${BLOCK} -s ${summ} -H ${herit} -n 300 -G ${valg} -P ${valp} -l ${col} -T ${type}  -i ${index} -t ${thread} -o ${outpath}
````
If the user wants to change the fold of heritability, you can revise the row 78, 79, 99 and 102.
You should use the output file of `ldsc` as `-H` parameter of `DBSLMM`.
The download link of `dbslmm` is <https://drive.google.com/file/d/1TRqPozXtenDW9buQgzFQ2KSxxWhl7dNn/view?usp=sharing>.

## Tutorial for DBSLMM (v0.2)
In this version, we treat the heritability as the tuning parameter. We give the bash script for the DBSLMM-tuning as following:
````bash
PACK_DIR=/your/path/
SPLITCHR=${PACK_DIR}DBSLMM/software/SPLITCHR.R
DBSLMM=${PACK_DIR}DBSLMM/software/DBSLMM.R
dbslmm=${PACK_DIR}DBSLMM/software/dbslmm
TUNE=${PACK_DIR}DBSLMM/software/TUNE.R
PLINK=${PACK_DIR}plink_1.9/plink

## SNP number
m=`cat ${summary_file_prefix}.assoc.txt | wc -l`

## sample size of GWAS
nobs=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summary_file_prefix}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)

## DBSLMM validation
Rscript ${SPLITCHR} --summary ${summary_file_prefix}.assoc.txt
for chr in `seq 1 22`
do
BLOCK=${PACK_DIR}DBSLMM/LDblock/chr${chr}
summchr=${summary_file_prefix}_chr${chr}

for h2f in 0.8 1 1.2
do
Rscript ${DBSLMM} --summary ${summchr}.assoc.txt --outPath ${outpath} --plink ${PLINK} --dbslmm ${dbslmm} --ref ${val_geno} \
--n ${n} --type ${type} --nsnp ${m} --block ${BLOCK}.bed --h2 ${herit} --h2f ${h2f} --thread ${thread}

summchr_prefix=`echo ${summchr##*/}`
${PLINK} --silent --bfile ${val_geno} --score ${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.txt 1 2 4 sum\
          --out ${outpath}${summchr_prefix}_h2f${h2f}

if [ -f "${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.badsnps" ];then
rm ${outpath}${summchr_prefix}_h2f${h2f}.dbslmm.badsnps
fi
rm ${outpath}${summchr_prefix}_h2f${h2f}.log
if [ -f "${outpath}${summchr_prefix}_h2f${h2f}.nopred" ];then
rm ${outpath}${summchr_prefix}_h2f${h2f}.nopred
fi
done
done

## DBSLMM selection 
summchr_prefix2=`echo ${summchr_prefix%_*}`
if [ ${cov} -eq 0 ]
then 
Rscript ${TUNE} --phenoPred ${outpath}${summchr_prefix2} --phenoVal ${val_pheno},${col} \
       --h2Range 0.8,1,1.2 --index ${index}
else 
Rscript ${TUNE} --phenoPred ${InterPred} --phenoVal ${val_pheno},${col} \
       --h2Range 0.8,1,1.2 --index ${index} --cov ${cov}
fi

## 
hbest=`cat ${outpath}${summchr_prefix2}_hbest.${index}`
for chr in `seq 1 22`
do
mv ${outpath}${summchr_prefix2}_chr${chr}_h2f${hbest}.dbslmm.txt ${outpath}${summchr_prefix2}_chr${chr}_best.dbslmm.txt
rm ${outpath}${summchr_prefix2}_chr${chr}_h2f*
done
````
We recommend the folds of heritability are set as 0.8, 1 and 1. The setting of them is flexiable as your data. If you have other setting, you can e-mail me. I will change them. The deault version is as following: 
````bash
Rscript ${DBSLMM} --summary ${summ}.assoc.txt --outPath ${outPath} --plink ${plink}\
                  --dbslmm ${dbslmm} --ref ${ref} --n ${n} --nsnp ${m} --block ${blockf}.bed\
                  --h2 ${h2} --thread ${thread}
````

## Tutorial for external test only using summary statistics
### Make SNP correlation matrix
We should use reference panel to construct SNP correlation matrix. You can treat [1000 Genome Project](https://www.internationalgenome.org/data#download) as reference panel. We also need the [block information](http://bitbucket.org/nygcresearch/ldetect-data). The code is `MAKELD.R`. 
````bash
mkld=/your/path/MAKELD.R
# construct SNP corrlation matrix for each chromosome
for chr in `seq 1 22`
do
bfile=/your/path/chr${chr}
blockfile=/your/path/chr${chr}.bed
plink=/your/path/plink-1.9
outpath=/your/path/outpath/
Rscript ${mkld} --bfile ${bfile} --blockfile ${blockfile} --plink ${plink} --outpath ${outpath} --chr ${chr}
done
````
### Estimate R<sup>2
The format of `extsumm` is the same as that of `LDSC`. 
````bash
SNP N Z A1 A2
rs1983865 253135 3.79310344827586 T C
rs1983864 251364 -4.51612903225806 T G
rs12411954 253213 0.413793103448276 T C
rs7077266 250092 -0.92 T G
````
We should use the standardized SNP effect size. You have two options. First, the variant ID is read from column 1, an allele code is read from the following column, the standardized effect size associated with the named allele is read from the column after the allele column. E.g.
````bash
esteff=/your/path/esteff.txt,1,2,3
````
Second, the variant ID is read from column 1, an allele code is read from the following column, the standardized effect size associated with the named allele is read from the column after the allele column, the MAF with the named allele is read from the last column. E.g.
````bash
esteff=/your/path/esteff.txt,1,2,3,4
````
The example of `EXTERNAL.R` is as following:
````bash
external=/your/path/EXTERNAL.R
extsumm=/your/path/extsumm.ldsc.gz
esteff=/your/path/esteff.txt,1,2,3
LDpath=/your/path/LDpath/
outpath=/your/path/outpath/
Rscript ${external} --extsumm ${extsumm} --esteff ${esteff} --LDpath ${LDpath} --outpath ${outpath}
````
