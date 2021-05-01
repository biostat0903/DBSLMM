# DBSLMM
Deterministic Bayesian Sparse Linear Mixed Model

## Update log
### v0.3
* fixes the bug of model fiting
### v0.21
* fits DBSLMM in different heritabiliies with one clumping
* uses a script to make the usage of DBSLMM-tuning version more easier.
### v0.2 
* validates the flods of heritability: 0.8, 1 and 1.2
* fits the LMM, when the chromosome without any large effect SNPs
* fits the external validation all by R code. 

## Tutorial for DBSLMM (v0.3)
In this version, we use `DBSLMM_script.sh` to make the usage of DBSLMM more easier. v0.3 needs `valg` to be the whole genome plink file. By parameter `-m`, we give two different model assumptions, including DBSLMM and LMM. Specially, LMM model assumption only have the default version, because it do not need to tune `-h2f`.
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
model=DBSLMM
sh ${DBSLMM} -D ${PACK_DIR} -p ${PLINK} -B ${BLOCK} -s ${summ} -H ${herit} -m ${model} -G ${valg} -P ${valp} -l ${col} -T ${type}  -i ${index} -t ${thread} -o ${outpath}

# DBSLMM determinitic version (without covariates)
type=d
model=DBSLMM
sh ${DBSLMM} -D ${PACK_DIR} -p ${PLINK} -B ${BLOCK} -s ${summ} -H ${herit} -m ${model} -G ${valg} -P ${valp} -l ${col} -T ${type}  -i ${index} -t ${thread} -o ${outpath}

# LMM version
type=d
model=LMM
sh ${DBSLMM} -D ${PACK_DIR} -p ${PLINK} -B ${BLOCK} -s ${summ} -H ${herit} -m ${model} -G ${valg} -P ${valp} -l ${col} -T ${type}  -i ${index} -t ${thread} -o ${outpath}
````
If the user wants to change the fold of heritability, you can revise the setting in `DBSLMM_script.sh`.
You should use the output file of `ldsc` as `-H` parameter of `DBSLMM`.
The download link of `dbslmm` is <https://drive.google.com/file/d/1Rpo2pNv4bXgWYiHnDZQfYeD0GUSJFBAM/view?usp=sharing>.

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
