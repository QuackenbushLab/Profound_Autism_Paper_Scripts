phenotypeDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/SSC\ Version\ 15.3\ Phenotype\ Dataset/Proband\ Data/"
outDir <- "/Users/tae771/Library/CloudStorage/OneDrive-HarvardUniversity/Documents/postdoc/SFARI/profoundAutism"
ssc <- read.csv(paste0(phenotypeDir, "ssc_core_descriptive.csv"), row.names = 1)
ados1 <- read.csv(paste0(phenotypeDir, "ados_1.csv"), row.names = 1)
ados2 <- read.csv(paste0(phenotypeDir, "ados_2.csv"), row.names = 1)
ados3 <- read.csv(paste0(phenotypeDir, "ados_3.csv"), row.names = 1)
ados4 <- read.csv(paste0(phenotypeDir, "ados_4.csv"), row.names = 1)

# Breakdown by verbal abilities and ID.
moderateToSevereID <- rownames(ssc)[which(ssc$ssc_diagnosis_full_scale_iq < 50)]
mildID <- rownames(ssc)[intersect(which(ssc$ssc_diagnosis_full_scale_iq >= 50),
                                  which(ssc$ssc_diagnosis_full_scale_iq < 70))]
noID <- rownames(ssc)[intersect(which(ssc$ssc_diagnosis_full_scale_iq >= 70),
                                  which(ssc$ssc_diagnosis_full_scale_iq < 115))]
gifted <- rownames(ssc)[which(ssc$ssc_diagnosis_full_scale_iq >= 115)]
profoundAutism <- intersect(moderateToSevereID, rownames(ados1))
moderateToSevereIDWithPhraseSpeech <- intersect(moderateToSevereID, rownames(ados2))
moderateToSevereIDWithVerbalFluency <- intersect(moderateToSevereID, rownames(ados3))
moderateToSevereIDWithVerbalFluencyOlderAdolescent <- intersect(moderateToSevereID, rownames(ados4))
nonverbalMildID <- intersect(mildID, rownames(ados1))
phraseSpeechMildID <- intersect(mildID, rownames(ados2))
verbalFluencyMildID <- intersect(mildID, rownames(ados3))
verbalFluencyMildIDOlderAdolescent <- intersect(mildID, rownames(ados4))
nonverbalNoID <- intersect(noID, rownames(ados1))
phraseSpeechNoID <- intersect(noID, rownames(ados2))
verbalFluencyNoID <- intersect(noID, rownames(ados3))
verbalFluencyNoIDOlderAdolescent <- intersect(noID, rownames(ados4))
nonverbalGifted <- intersect(gifted, rownames(ados1))
phraseSpeechGifted <- intersect(gifted, rownames(ados2))
verbalFluencyGifted <- intersect(gifted, rownames(ados3))
verbalFluencyGiftedOlderAdolescent <- intersect(gifted, rownames(ados4))

# Save SSC data for the four groups.
write.csv(ssc[profoundAutism,], paste0(outDir, "/profoundAutism_SSC.csv"))
write.csv(ssc[c(moderateToSevereIDWithPhraseSpeech,moderateToSevereIDWithVerbalFluency,
                moderateToSevereIDWithVerbalFluencyOlderAdolescent),], 
          paste0(outDir, "/moderateToSevereIDVerbal_SSC.csv"))
write.csv(ssc[nonverbalMildID,], 
          paste0(outDir, "/mildIDnonverbal_SSC.csv"))
write.csv(ssc[c(phraseSpeechMildID,verbalFluencyMildID,
                verbalFluencyMildIDOlderAdolescent),], 
          paste0(outDir, "/mildIDVerbal_SSC.csv"))
write.csv(ssc[nonverbalNoID,], 
          paste0(outDir, "/noIDnonverbal_SSC.csv"))
write.csv(ssc[c(phraseSpeechNoID,verbalFluencyNoID,
                verbalFluencyNoIDOlderAdolescent),], 
          paste0(outDir, "/noIDVerbal_SSC.csv"))
write.csv(ssc[nonverbalGifted,], 
          paste0(outDir, "/giftedNonverbal_SSC.csv"))
write.csv(ssc[c(phraseSpeechGifted,verbalFluencyGifted,
                verbalFluencyGiftedOlderAdolescent),], 
          paste0(outDir, "/giftedVerbal_SSC.csv"))



