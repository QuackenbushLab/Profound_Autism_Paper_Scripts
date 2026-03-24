#!/bin/bash
#SBATCH --job-name=eqtlRand
#SBATCH --time=20:00:00 
#SBATCH --mem=10GB 
#SBATCH --output=eqtlRand.log 
#SBATCH --error=eqtlRand_error.log 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=teicher@hsph.harvard.edu
#SBATCH --gres=gpu:1

conda activate eQTL_env
for i in 1 2 3 4 5 6 7 8 9 10; do
	for j in 1 2; do
		python3 -m tensorqtl --covariates ../eQTL/covarNewProfoundBoth_${i}_${j}.txt --mode cis_nominal ../eQTL/profoundAutismBoth_rand_${i}_${j} ../eQTL/expression.pheno.profoundBoth_${i}_${j}_filt.bed ../eQTL/profoundAutismBoth_rand_${i}_${j}
		python3 -m tensorqtl --covariates ../eQTL/covarNewProfoundBoth_${i}_${j}.txt --mode trans ../eQTL/profoundAutismBoth_rand_${i}_${j} ../eQTL/expression.pheno.profoundBoth_${i}_${j}_filt.bed ../eQTL/profoundAutismBoth_rand_${i}_${j}_trans
		python3 -m tensorqtl --covariates ../eQTL/covarNewProfoundBoth_${i}_${j}_diff.txt --mode cis_nominal ../eQTL/profoundAutismBoth_rand_${i}_${j} ../eQTL/expression.pheno.profoundBoth_${i}_${j}_diff_filt.bed ../eQTL/profoundAutismBoth_rand_${i}_${j}_diff
		python3 -m tensorqtl --covariates ../eQTL/covarNewProfoundBoth_${i}_${j}_diff.txt --mode trans ../eQTL/profoundAutismBoth_rand_${i}_${j} ../eQTL/expression.pheno.profoundBoth_${i}_${j}_diff_filt.bed ../eQTL/profoundAutismBoth_rand_${i}_${j}_diff_trans
	done
done
