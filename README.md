# Replicating Profound Autism Genomics, Transcriptomics, and Regulomics analysis.

## Basic Setup
1. Install [R](https://www.r-project.org/).
2. Install [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/).
3. Install [Globus Connect Personal](https://www.globus.org/globus-connect-personal).
4. If you are not working in a Unix environment or do not have **awk** pre-installed, install [awk](https://www.gnu.org/software/gawk/manual/html_node/Installation.html).
5. Install [BEDOPS](https://bedops.readthedocs.io/en/latest/content/installation.html).
6. Install [Python](https://www.python.org/downloads/).
7. Install and setup [tensorqtl and Pgenlib](https://github.com/broadinstitute/tensorqtl).

## Accessing Data
1.  Request access to the Simons Simplex Collection RNAseq and phenotype data (pilot and large data set) using [SFARI Base](https://base.sfari.org/).
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
4. To compute the standard difference in gene expression between probands and non-autistic siblings, run **diffExpressionAgainstSiblingsLarge.R**.
    - *sourceDir*: The directory where the expression data files are stored.
5. To subset the gene expression data for each PA and ASD subgroup, run **splitExpressionBySubgroup.R**.
    - *phenoGrp*: The directory where the phenotype groups (generated in step 1) are stored.
    - *sourceDirGenomics*: The directory where the logCPM gene expression matrix is stored.
6. To subset the standard difference in gene expression for each PA and ASD subgroup, run **splitDiffExpressionBySubgroup.R**.
    - *phenoGrp*: The directory where the phenotype groups (generated in step 1) are stored.
    - *sourceDirGenomics*: The directory where the standard difference in gene expression matrix is stored.

## Computing Separability Along Axes of Variation
1. To run PCA on the logCPM data and perform a Wilcoxon rank-sum test on each PC, run **genomicsPCA.R**. This will also generate Supplementary Figure 1.
    - *outDir*: The directory where the phenotype groups (generated in step 1) are stored.
    - *genomicsDir*: The directory where the gene expression matrix is stored.
2. To run PCA on the standard difference data and perform a Wilcoxon rank-sum test on each PC, run **genomicsDiffPCA.R**. This will also generate Figure 2.
    - *outDir*: The directory where the phenotype groups (generated in step 1) are stored.
    - *genomicsDir*: The directory where the standard difference in gene expression matrix is stored.

## Running GWAS
1. To calculate the first 50 PCs from the SNP data, run the following command with the files in the angle brackets specified accordingly: ```plink2 --pca 50 --bfile <plink-file-path> --out <pca-file-path> --maf 0.05```.
2. To generate the PCA elbow plot (Supplementary Figure 2), run **summarizeEigenvaluesGeneVariantPCA.R**.
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
   - *ref*: Path to the GENCODE reference file.

## Running eQTL Analysis
1. To filter samples for eQTL analysis, run the following command for each PA / ASD subgroup to generate a PLINK 2 binary genotype table, filling in the items in angle brackets accordingly: ```plink2 --bfile <plink-file-path> --keep <samples-in-subgroup-file> --make-pgen --out <subgroup-file>```.
2. Remove the family ID from the file for compatibility, filling in the items in angle brackets accordingly: ```awk '{sub(/^[^[:space:]]+[[:space:]]*/, ""); print}' <subgroup-file>.psam > <subgroup-file-without-family-id>.psam```
3. To convert the logCPM gene expression data to BED format, run **convertExpressionToBED.R**.
    - *expr_csv*: Path to the CSV file containing logCPM expression data.
    - *gtf_path*: Path to the GENCODE reference file.
    - *geno_prefix*: Path (with prefix) to the PLINK results for a given PA / ASD subgroup.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *out_tsv*: Path to the BED file where you wish to save the gene expression data.
4. To convert the standard difference gene expression data to BED format, run **convertDiffExpressionToBED.R**.
    - *expr_csv*: Path to the CSV file containing logCPM expression data.
    - *gtf_path*: Path to the GENCODE reference file.
    - *geno_prefix*: Path (with prefix) to the PLINK results for a given PA / ASD subgroup.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *out_tsv*: Path to the BED file where you wish to save the gene expression data.
5. To format the covariates for eQTL analysis, including only samples with logCPM data, run **formatCovariatesForEQTL.R**.
    - *clinDir*: Directory where phenotype data are stored.
    - *covar*: File where covariates (the same ones used in the GWAS) are stored.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *expressionBed*: Path to the BED file where the expression data are stored.
6. To format the covariates for eQTL analysis, including only samples with standard difference data, run **formatDiffCovariatesForEQTL.R**.
    - *clinDir*: Directory where phenotype data are stored.
    - *covar*: File where covariates (the same ones used in the GWAS) are stored.
    - *mergedFile*: Path to the file containing the expression sample to SNP sample mapping.
    - *expressionBed*: Path to the BED file where the expression data are stored.
7. To compute cis-eQTLs for each PA / ASD subgroup, run the following command, filling in the items in angle brackets accordingly: ```python3 -m tensorqtl --covariates <eqtl-covariate-file> --mode cis_nominal <subgroup-file-path> <gene-expression-bed-file> <subgroup-result-file>```
8. To compute trans-eQTLs for each PA / ASD subgroup, run the following command, filling in the items in angle brackets accordingly: ```python3 -m tensorqtl --covariates <eqtl-covariate-file> --mode trans <subgroup-file-path> <gene-expression-bed-file> <subgroup-result-file>```
9. To combine all cis-and-trans-eQTL results, run **processParquetResults.R**.
    - *eQTLDir*: Path to directory where eQTL results are stored (subgroup-result-file).
    - *cutoff*: Cutoff for significance.
10. To generate combination matrices for eQTL UpSet plots, run **eQTLMatrices.R**.
    - *eQTLDir*: Path to eQTL results.
11. To generate the eQTL UpSet plots, run **eQTLFigures.R**.
    - *eQTLDir*: Path to eQTL results.

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
Run **MakeCompositeSuppTables_v2.R**.
    - *diffExprDiffDir*: Path to the directory where the linear models for the standard difference analyses are stored.
    - *diffExprRawDir*: Path to the directory where the linear models for the logCPM gene expression analyses are stored.
    - *diffVarDiffDir*: Path to the directory where the Levene's test results for the standard difference data are stored.
    - *diffVarRawDir*: Path to the directory where the Levene's test results for the logCPM gne expression data are stored.
    - *litmanDir*: Path to the directory where the GSEA results for the Litman sets are stored.
    
    
