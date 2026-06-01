# Replicating Profound Autism Genomics, Transcriptomics, and Regulomics Analysis

## Basic Setup
1. Install [R](https://www.r-project.org/).
2. Install [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/).
3. Install [Globus Connect Personal](https://www.globus.org/globus-connect-personal).
4. If you are not working in a Unix environment or do not have **awk** pre-installed, install [awk](https://www.gnu.org/software/gawk/manual/html_node/Installation.html).
5. Install [BEDOPS](https://bedops.readthedocs.io/en/latest/content/installation.html).
6. Install [Python](https://www.python.org/downloads/).
7. Install and setup [tensorqtl and Pgenlib](https://github.com/broadinstitute/tensorqtl).

## Accessing Data
1.  Request access to the Simons Simplex Collection RNAseq and phenotype data (large data) using [SFARI Base](https://base.sfari.org/).
2.  Register for a Globus account.
3.  Request and download the phenotype data directly from SFARI Base.
4.  Request and download the RNASeq_Large data using Globus Connect Personal. If at Harvard, these can be transferred directly to [FASRC](https://docs.rc.fas.harvard.edu/kb/globus-file-transfer/).
5.  Request and download the Structural Variants Dataset (SSC, SPARK, RARE) using Globus Connect Personal. You will only need the PLINK-formatted data, not the FASTQ or BAM files.

## Preprocessing Data
1. To preprocess the phenotypic data, run **phenotypeProcessing.R**, modifying the following parameters:
    - *phenotypeDir*: The directory where the phenotype data are stored.
    - *outDir*: The directory where the output groups should be stored.
2. To filter phenotypic data by age and generate the demographic table (Table 1), run **StatSummaryAbove8.R**.
    - *phenotypeDir*: The directory where the phenotype data are stored.
    - *outDir*: The directory where the phenotype groups are stored.
3. To compile and normalize gene expression data, run **compileGeneExpressionData.R**.
    - *expressionDir*: The directory where the expression data files are stored.
4. To compute logCPM_PSE, run **diffExpressionAgainstSiblingsLarge.R**.
    - *sourceDir*: The directory where the expression data files are stored.
5. To subset the logCPM gene expression data for each PA and ASD subgroup, run **splitExpressionBySubgroup.R**.
    - *phenoGrp*: The directory where the phenotype groups (generated in step 1) are stored.
    - *sourceDirGenomics*: The directory where the logCPM gene expression matrix is stored.
6. To subset the logCPM_PSE for each PA and ASD subgroup, run **splitDiffExpressionBySubgroup.R**.
    - *phenoGrp*: The directory where the phenotype groups (generated in step 1) are stored.
    - *sourceDirGenomics*: The directory where the standard difference in gene expression matrix is stored.

## Computing Separability Along Eigengenes
1. To run PCA on the logCPM data and perform a Wilcoxon rank-sum test on each PC, run **genomicsPCA.R**.
    - *outDir*: The directory where the phenotype groups are stored.
    - *genomicsDir*: The directory where the gene expression matrix is stored.
2. To run PCA on the logCPM_PSE data and perform a Wilcoxon rank-sum test on each PC, run **genomicsDiffPCA.R**.
    - *outDir*: The directory where the phenotype groups are stored.
    - *genomicsDir*: The directory where the standard difference in gene expression matrix is stored.

## Running GWAS
1. To calculate the first 50 PCs from the SNP data, run the following command with the files in the angle brackets specified accordingly: ```plink2 --pca 50 --bfile <plink-file-path> --out <pca-file-path> --maf 0.05```.
2. To generate the PCA elbow plot, run **summarizeEigenvaluesGeneVariantPCA.R**.
    - *inDir*: The directory where the PLINK PCA was saved (pca-file-path).
3. To remove the header from the eigenvector file, run the following command, filling in the items in angle brackets accordingly: (**pca-eigenvector-file** should have been generated in step 1): ```awk 'NR > 1' <pca-eigenvector-file> <pca-eigenvector-file-without-header>```
4. To format the PCs as covariates for input to PLINK, run the following command, filling in the items in angle brackets accordingly: ```awk '{print $1, $2, $3, $4, $5}' <pca-eigenvector-file-without-header> > <covariate-file>```
5. To run the PLINK logistic regression models, run ```sbatch profoundAutismRunPLINK_pcs.sh``` in a SLURM environment. If you do not have access to a SLURM environment, run the commands sequentially.
6. To filter for significant SNPs and calculate SNP rankings (Table 2), run **summarizePLINKProfoundAutismResultsPCA.R**. 
  - *sourceDir*: The directory where the PLINK model results were stored (see step 5).
7. Download the basic gene annotation from [GENCODE Human Release 49](https://www.gencodegenes.org/human/).
8. To map the SNPs to genes, run the following command, filling in the items in angle brackets accordingly: ```gtf2bed < <gtf-gencode-file> > <bed-gencode-file>```.
7. To generate an UpSet plot of the overlapping SNPs (Figure 3) and map the SNPs to genes, run **GWAS_figure_pcs.R**.
   - *sourceDir*: The directory where the PLINK results are stored.
   - *refFile*: Path to the GENCODE reference file.

## Running eQTL Analysis
1. To filter samples for eQTL analysis, run the following command for each PA / ASD subgroup to generate a PLINK 2 binary genotype table, filling in the items in angle brackets accordingly: ```plink2 --bfile <plink-file-path> --keep <samples-in-subgroup-file> --make-pgen --out <subgroup-file>```. Note that samples-in-subgroup-file corresponds to the 
2. Remove the family ID from the file for compatibility, filling in the items in angle brackets accordingly: ```awk '{sub(/^[^[:space:]]+[[:space:]]*/, ""); print}' <subgroup-file>.psam > <subgroup-file-without-family-id>.psam```
3. To convert the logCPM gene expression data to BED format, run **convertExpressionToBED.R**.
    - *expr_csv*: Path to the CSV file containing logCPM expression data.
    - *gtf_path*: Path to the GENCODE reference file.
    - *geno_prefix*: Path (with prefix) to the PLINK results.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *out_tsv*: Path to the directory where you wish to save the gene expression BED file.
4. To convert the logCPM_PSE expression data to BED format, run **convertDiffExpressionToBED.R** for each phenotype group.
    - *expr_csv*: Path to the CSV file containing standard difference expression data.
    - *gtf_path*: Path to the GENCODE reference file.
    - *geno_prefix*: Path (with prefix) to the PLINK results for the phenotype group.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *out_tsv*: Path to the directory where you wish to store the gene expression BED file for the phenotype group.
5. To format the covariates for eQTL analysis, including only samples with logCPM data, run **formatCovariatesForEQTL.R**.
    - *clinDir*: Directory where phenotype data are stored.
    - *covarFile*: Path to the file where covariates (the same ones used in the GWAS) are stored.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *expressionBed*: Path to the BED file where the expression data are stored.
6. To format the covariates for eQTL analysis, including only samples with standard difference data, run **formatDiffCovariatesForEQTL.R**, for each phenotype group.
    - *clinDir*: Directory where phenotype data are stored.
    - *covarFile*: File where covariates (the same ones used in the GWAS) are stored.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *expressionBed*: Path to the BED file where the expression data are stored.
    - *covarNewOut*: File where the newly formatted covariates should be stored.
    - *expressionBedNewOut*: File where the newly formatted expression BED file should be stored.
7. To compute cis-eQTLs for each PA / ASD subgroup, run the following command, filling in the items in angle brackets accordingly: ```python3 -m tensorqtl --covariates <eqtl-covariate-file> --mode cis_nominal <subgroup-file-path> <gene-expression-bed-file> <subgroup-result-file>```
8. To compute trans-eQTLs for each PA / ASD subgroup, run the following command, filling in the items in angle brackets accordingly: ```python3 -m tensorqtl --covariates <eqtl-covariate-file> --mode trans <subgroup-file-path> <gene-expression-bed-file> <subgroup-result-file>```
9. To combine all cis-and-trans-eQTL results from logCPM data, run **processParquetResults.R**.
    - *eQTLDir*: Path to directory where eQTL results are stored (subgroup-result-file).
    - *cutoff*: Cutoff for significance.
9. To combine all cis-and-trans-eQTL results from logCPM_PSE data, run **processParquetResultsDiff.R**.
    - *eQTLDir*: Path to directory where eQTL results are stored (subgroup-result-file).
    - *cutoff*: Cutoff for significance.
10. To generate combination matrices for eQTL UpSet plots, run **eQTLMatrices.R**.
    - *eQTLDir*: Path to eQTL results.
11. To generate the eQTL UpSet plots, run **eQTLFigures.R**.
    - *eQTLDir*: Path to eQTL results.

## Running eQTL Analysis on Randomized Subsets
1. Run **randomSubsets.R**, changing the following variables accordingly:
    - *paDir*: Path to the directory where the PA / ASD subgroup files are stored.
    - *eqtlDir*: Path to the directory where eQTL inputs and results are stored.
2. Run ```./makePgenForPLINKRand.sh <paDir> <eqtlDir> <plinkFilePath>```
3. Run ```./removeFamilyID.sh <eqtlDir>```.
4. Run **convertExpressionToBED_rand.R** and **convertExpressionToBED_rand_diff.R**, changing the following paths as needed:
    - *expr_csv*: Path to the CSV file containing logCPM / pro band-specific expression data.
    - *gtf_path*: Path to the GENCODE reference file.
    - *geno_prefix_dir*: Path (with prefix) to the PLINK results for a given PA / ASD subgroup.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *out_tsv_dir*: Path to the BED file where you wish to save the gene expression data.
5. Run **formatCovariatesForEqtl_rand.R** and **formatDiffCovariatesForEqtl_rand.R* to format the covariates.
6. Run **sbatch runEqtlAnalysisRand.sh** in a SLURM environment.
7. Run **processParquetResults_rand.R** and **processParquetResults_diff_rand.R**, modifying the following:
    - *eQTL*: Path to the eQTL results.
8. Run **calculateJaccard_rand.R**.
8. Run **jaccardSummaryEqtl.R**.


# Running Differential Gene Expression Analysis
1. To run the linear models for PA / ASD subgroup comparison using standard difference gene expression data, run **linearModelsDiff.R**.
    - *sourceDirGenomics*: Source directory where gene expression data are stored.
    - *phenoGrp*: Source directory where PA /ASD subgroup files are stored.
    - *outDir*: Directory where you wish to store the linear model results.
2. To run the linear models for PA / ASD subgroup comparison using logCPM gene expression data, run **linearModels.R**.
    - *sourceDirGenomics*: Source directory where gene expression data are stored.
    - *phenoGrp*: Source directory where PA /ASD subgroup files are stored.
    - *outDir*: Directory where you wish to store the linear model results.

# Running Differential Gene Variance Analysis
1. To run the Levene's Test for subgroup comparison using standard difference gene expression data, run **leveneDiff.R**.
    - *sourceDirGenomics*: Source directory where gene expression data are stored.
    - *phenoGrp*: Source directory where PA /ASD subgroup files are stored.
    - *outDirFinal*: Directory where you wish to store the linear model results.
2. To run the Levene's Test for subgroup comparison using logCPM gene expression data, run **levene.R**.
    - *sourceDirGenomics*: Source directory where gene expression data are stored.
    - *phenoGrp*: Source directory where PA /ASD subgroup files are stored.
    - *outDirFinal*: Directory where you wish to store the linear model results.
3. To run the GSEA for the nonspeaking group, run **GSEA.R**.
    - *msigdbPath*: Path to GMT file from MSigDB.
    - *diffExpressionPath*: Path to standard difference in expression files for each PA / ASD subgroup.

# Running Litman Gene Set Enrichment Analysis
1. To run the linear models for IQ, Speech, and Combination analyses using standard difference gene expression data, run **linearModelsLitmanDiff.R**.
    - *sourceDirGenomics*: Source directory where gene expression data are stored.
    - *phenoGrp*: Source directory where PA /ASD subgroup files are stored.
    - *outDir*: Directory where you wish to store the linear model results.
2. To run the linear models for PA / ASD subgroup comparison using logCPM gene expression data, run **linearModelsLitman.R**.
    - *sourceDirGenomics*: Source directory where gene expression data are stored.
    - *phenoGrp*: Source directory where PA /ASD subgroup files are stored.
    - *outDir*: Directory where you wish to store the linear model results.
3. To generate the heat maps, run **LitmanSubgroupAnalysisAggregated.R**.
    - *sourceDir*: Source directory where the Litman gene sets are stored.
# Consolidating Results for Supplementary Tables
1. Run **MakeCompositeSuppTables_v2.R**.
    - *diffExprDiffDir*: Path to the directory where the linear models for the standard difference analyses are stored.
    - *diffExprRawDir*: Path to the directory where the linear models for the logCPM gene expression analyses are stored.
    - *diffVarDiffDir*: Path to the directory where the Levene's test results for the standard difference data are stored.
    - *diffVarRawDir*: Path to the directory where the Levene's test results for the logCPM gne expression data are stored.
    - *litmanDir*: Path to the directory where the GSEA results for the Litman sets are stored.
2. Run **profoundAutism_LitmanBinaryTable.R**.
    
    
