diffExprDiffDir <- NULL
diffExprRawDir <- NULL
diffVarRawDir <- NULL
diffVarDiffDir <- NULL
litmanDir <- NULL

bothGiftedRaw <- read.csv(paste0(diffExprRawDir, "profoundBoth_GiftedVerbal.csv"), row.names = 1)
bothMildRaw <- read.csv(paste0(diffExprRawDir, "profoundBoth_MildIDVerbal.csv"), row.names = 1)
bothNoRaw <- read.csv(paste0(diffExprRawDir, "profoundBoth_NoIDVerbal.csv"), row.names = 1)
idGiftedRaw <- read.csv(paste0(diffExprRawDir, "profoundModerateIDOnly_GiftedVerbal.csv"), row.names = 1)
idMildRaw <- read.csv(paste0(diffExprRawDir, "profoundModerateIDOnly_MildIDVerbal.csv"), row.names = 1)
idNoRaw <- read.csv(paste0(diffExprRawDir, "profoundModerateIDOnly_NoIDVerbal.csv"), row.names = 1)
idNonverbalRaw <- read.csv(paste0(diffExprRawDir, "profoundModerateIDOnly_NonverbalOnly.csv"), row.names = 1)
idBothRaw <- read.csv(paste0(diffExprRawDir, "profoundModerateIDOnly_profoundBoth.csv"), row.names = 1)
nonverbalGiftedRaw <- read.csv(paste0(diffExprRawDir, "profoundNonverbalOnly_GiftedVerbal.csv"), row.names = 1)
nonverbalMildRaw <- read.csv(paste0(diffExprRawDir, "profoundNonverbalOnly_MildIDVerbal.csv"), row.names = 1)
nonverbalNoRaw <- read.csv(paste0(diffExprRawDir, "profoundNonverbalOnly_NoIDVerbal.csv"), row.names = 1)
nonverbalBothRaw <- read.csv(paste0(diffExprRawDir, "profoundNonverbalOnly_profoundBoth.csv"), row.names = 1)
consolidatedDERaw <- data.frame(gene = bothGiftedRaw$gene, paBoth_asdGifted = bothGiftedRaw$padj,
                                paBoth_asdMildID = bothMildRaw$padj,
                                paBoth_asdNoID = bothNoRaw$padj,
                                paID_asdGifted = idGiftedRaw$padj,
                                paID_asdMild = idMildRaw$padj,
                                paID_asdNo = idNoRaw$padj,
                                paNonspeaking_asdGifted = nonverbalGiftedRaw$padj,
                                paNonspeaking_asdMild = nonverbalMildRaw$padj,
                                paNonspeaking_asdNo = nonverbalNoRaw$padj,
                                paID_paBoth = idBothRaw$padj,
                                paNonspeaking_paID = idNonverbalRaw$padj)
write.csv(consolidatedDERaw, paste0(diffExprRawDir, "consolidatedPvals.csv"))

bothGiftedDiff <- read.csv(paste0(diffExprDiffDir, "profoundAutismBoth_GiftedVerbal.csv"), row.names = 1)
bothMildDiff <- read.csv(paste0(diffExprDiffDir, "profoundAutismBoth_MildIDVerbal.csv"), row.names = 1)
bothNoDiff <- read.csv(paste0(diffExprDiffDir, "profoundAutismBoth_NoIDVerbal.csv"), row.names = 1)
idGiftedDiff <- read.csv(paste0(diffExprDiffDir, "profoundModerateIDOnly_GiftedVerbal.csv"), row.names = 1)
idMildDiff <- read.csv(paste0(diffExprDiffDir, "profoundModerateIDOnly_MildIDVerbal.csv"), row.names = 1)
idNoDiff <- read.csv(paste0(diffExprDiffDir, "profoundModerateIDOnly_NoIDVerbal.csv"), row.names = 1)
idNonverbalDiff <- read.csv(paste0(diffExprDiffDir, "profoundModerateIDOnly_NonverbalOnly.csv"), row.names = 1)
idBothDiff <- read.csv(paste0(diffExprDiffDir, "profoundModerateIDOnly_profoundBoth.csv"), row.names = 1)
nonverbalGiftedDiff <- read.csv(paste0(diffExprDiffDir, "profoundNonverbalOnly_GiftedVerbal.csv"), row.names = 1)
nonverbalMildDiff <- read.csv(paste0(diffExprDiffDir, "profoundNonverbalOnly_MildIDVerbal.csv"), row.names = 1)
nonverbalNoDiff <- read.csv(paste0(diffExprDiffDir, "profoundNonverbalOnly_NoIDVerbal.csv"), row.names = 1)
nonverbalBothDiff <- read.csv(paste0(diffExprDiffDir, "profoundNonverbalOnly_profoundBoth.csv"), row.names = 1)
consolidatedDEDiff <- data.frame(gene = bothGiftedDiff$gene, paBoth_asdGifted = bothGiftedDiff$padj,
                                paBoth_asdMildID = bothMildDiff$padj,
                                paBoth_asdNoID = bothNoDiff$padj,
                                paID_asdGifted = idGiftedDiff$padj,
                                paID_asdMild = idMildDiff$padj,
                                paID_asdNo = idNoDiff$padj,
                                paNonspeaking_asdGifted = nonverbalGiftedDiff$padj,
                                paNonspeaking_asdMild = nonverbalMildDiff$padj,
                                paNonspeaking_asdNo = nonverbalNoDiff$padj,
                                paID_paBoth = idBothDiff$padj,
                                paNonspeaking_paID = idNonverbalDiff$padj)
write.csv(consolidatedDEDiff, paste0(diffExprDiffDir, "consolidatedPvals.csv"))

# Do the differential variance.
bothRawVar <- read.csv(paste0(diffVarRawDir, "profoundAutismBoth_NotProfound.csv"), row.names = 1)
idRawVar <- read.csv(paste0(diffVarRawDir, "profoundModerateIDOnly_NotProfound.csv"), row.names = 1)
nonverbalRawVar <- read.csv(paste0(diffVarRawDir, "profoundNonverbalOnly_NotProfound.csv"), row.names = 1)
idNonverbalVar <- read.csv(paste0(diffVarRawDir, "profoundModerateIDOnly_NonverbalOnly.csv"), row.names = 1)
idBothVar <- read.csv(paste0(diffVarRawDir, "profoundModerateIDOnly_ProfoundBoth.csv"), row.names = 1)
nonverbalBothVar <- read.csv(paste0(diffVarRawDir, "profoundNonverbalOnly_ProfoundBoth.csv"), row.names = 1)
consolidatedRawVar <- data.frame(gene = bothRawVar$gene,
                                 paBoth_asd = bothRawVar$padj,
                                 paID_asd = idRawVar$padj,
                                 paNonspeaking_asd = nonverbalRawVar$padj,
                                 paID_paNonspeaking = idNonverbalVar$padj,
                                 paID_paBoth = idBothVar$padj,
                                 paNonspeaking_paBoth = nonverbalBothVar$padj)
write.csv(consolidatedRawVar, paste0(diffVarRawDir, "consolidatedPvals.csv"))

bothDiffVar <- read.csv(paste0(diffVarDiffDir, "profoundAutismBoth_NotProfound.csv"), row.names = 1)
idDiffVar <- read.csv(paste0(diffVarDiffDir, "profoundModerateIDOnly_NotProfound.csv"), row.names = 1)
nonverbalDiffVar <- read.csv(paste0(diffVarDiffDir, "profoundNonverbalOnly_NotProfound.csv"), row.names = 1)
idNonverbalVar <- read.csv(paste0(diffVarDiffDir, "profoundModerateIDOnly_NonverbalOnly.csv"), row.names = 1)
idBothVar <- read.csv(paste0(diffVarDiffDir, "profoundModerateIDOnly_ProfoundBoth.csv"), row.names = 1)
nonverbalBothVar <- read.csv(paste0(diffVarDiffDir, "profoundNonverbalOnly_ProfoundBoth.csv"), row.names = 1)
consolidatedDiffVar <- data.frame(gene = bothDiffVar$gene,
                                 paBoth_asd = bothDiffVar$padj,
                                 paID_asd = idDiffVar$padj,
                                 paNonspeaking_asd = nonverbalDiffVar$padj,
                                 paID_paNonspeaking = idNonverbalVar$padj,
                                 paID_paBoth = idBothVar$padj,
                                 paNonspeaking_paBoth = nonverbalBothVar$padj)
write.csv(consolidatedDiffVar, paste0(diffVarDiffDir, "consolidatedPvals.csv"))


# Do the Litman sets.
full <- rbind(read.csv(paste0(litmanDir, "mean_geneSet_pvals_fullpw_pos.csv"), row.names = 1),
              read.csv(paste0(litmanDir, "mean_geneSet_pvals_fullpw_neg.csv"), row.names = 1))
full$direction <- c(rep("pos", nrow(full) / 2),
                    rep("neg", nrow(full) / 2))
write.csv(full, paste0(litmanDir, "consolidatedPvalsFull.csv"))

upw <- rbind(read.csv(paste0(litmanDir, "mean_geneSet_pvals_upw_pos.csv"), row.names = 1),
              read.csv(paste0(litmanDir, "mean_geneSet_pvals_upw_neg.csv"), row.names = 1))
upw$direction <- c(rep("pos", nrow(upw) / 2),
                    rep("neg", nrow(upw) / 2))
write.csv(upw, paste0(litmanDir, "consolidatedPvalsUpw.csv"))
