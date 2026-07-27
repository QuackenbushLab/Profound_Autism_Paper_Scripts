#!/bin/bash
#SBATCH --job-name=profoundAutismPLINK 
#SBATCH --time=24:00:00 
#SBATCH --mem=5GB 
#SBATCH --output=profoundAutismRunPLINK.log 
#SBATCH --error=profoundAutismRunPLINK_error.log 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=teicher@hsph.harvard.edu

#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_profoundBoth.txt --out ../PLINK/resultFiles/profoundBoth_noSex_pcs
plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2_missingFirst.txt --keep ../PLINK/generatedFiles/Omni2.5_profoundModerateIDOnly.txt --out ../PLINK/resultFiles/profoundModerateIDOnly_noSex_pcs
plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2_missingFirst.txt --keep ../PLINK/generatedFiles/Omni2.5_profoundNonverbalOnly.txt --out ../PLINK/resultFiles/profoundNonverbalOnly_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_verbalMildID.txt --out ../PLINK/resultFiles/verbalMildID_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_verbalNoID.txt --out ../PLINK/resultFiles/verbalNoID_noSex_pcs
plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2_missingFirst.txt --keep ../PLINK/generatedFiles/Omni2.5_verbalGifted.txt --out ../PLINK/resultFiles/verbalGifted_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_notProfoundNoGifted.txt --out ../PLINK/resultFiles/notProfoundNoGifted_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_profoundBothNonverbal.txt --out ../PLINK/resultFiles/profoundBothNonverbal_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_notProfoundNoID.txt --out ../PLINK/resultFiles/notProfoundNoID_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_notProfoundAll.txt --out ../PLINK/resultFiles/notProfoundAll_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_profoundBothID.txt --out ../PLINK/resultFiles/profoundBothID_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_verbalNoID.txt --out ../PLINK/resultFiles/verbalNoID_noSex_pcs
#plink2 --bfile ../PLINK/Omni2.5/SSC_Omni2.5 --logistic --covar ../PLINK/Omni_covar_2.txt --keep ../PLINK/generatedFiles/Omni2.5_profoundAll.txt --out ../PLINK/resultFiles/profoundAll_noSex_pcs

