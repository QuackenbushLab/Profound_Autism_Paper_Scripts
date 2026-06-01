# Read SSC data.
phenotypeDir <- NULL
ssc <- read.csv(paste0(phenotypeDir, "ssc_core_descriptive.csv"), row.names = 1)

# Read subgroups.
outDir <- NULL
profoundAutismSSC <- read.csv(paste0(outDir, "/profoundAutism_SSC.csv"), row.names = 1, stringsAsFactors = TRUE)
moderateToSevereVerbalSSC <- read.csv(paste0(outDir, "/moderateToSevereIDVerbal_SSC.csv"), row.names = 1, stringsAsFactors = TRUE)
nonverbalMildIDSSC <- read.csv(paste0(outDir, "/mildIDnonverbal_SSC.csv"), row.names = 1, stringsAsFactors = TRUE)
verbalMildIDSSC <- read.csv(paste0(outDir, "/mildIDVerbal_SSC.csv"), row.names = 1, stringsAsFactors = TRUE)
nonverbalNoIDSSC <- read.csv(paste0(outDir, "/noIDnonverbal_SSC.csv"), row.names = 1, stringsAsFactors = TRUE)
verbalNoIDSSC <- read.csv(paste0(outDir, "/noIDVerbal_SSC.csv"), row.names = 1, stringsAsFactors = TRUE)
nonverbalGiftedSSC <- read.csv(paste0(outDir, "/giftedNonverbal_SSC.csv"), row.names = 1, stringsAsFactors = TRUE)
verbalGiftedSSC <- read.csv(paste0(outDir, "/giftedVerbal_SSC.csv"), row.names = 1, stringsAsFactors = TRUE)

# Split into groups of interest.
profoundAutismModerateIDOnly <- moderateToSevereVerbalSSC
profoundAutismNonverbalOnly <- do.call(rbind, list(nonverbalMildIDSSC, nonverbalNoIDSSC,
                                                   nonverbalGiftedSSC))
profoundAutismBoth <- profoundAutismSSC
verbalMildID <- verbalMildIDSSC
verbalNoID <- verbalNoIDSSC
verbalGifted <- verbalGiftedSSC

# Subset to above 8.
profoundAutismModerateIDOnly <- profoundAutismModerateIDOnly[which(profoundAutismModerateIDOnly$age_at_ados / 12 > 8),]
write.csv(profoundAutismModerateIDOnly, paste0(outDir, "/profoundAutismModerateIDOnly_above8.csv"))
profoundAutismNonverbalOnly <- profoundAutismNonverbalOnly[which(profoundAutismNonverbalOnly$age_at_ados / 12 > 8),]
write.csv(profoundAutismNonverbalOnly, paste0(outDir, "/profoundAutismNonverbalOnly_above8.csv"))
profoundAutismBoth <- profoundAutismBoth[which(profoundAutismBoth$age_at_ados / 12 > 8),]
write.csv(profoundAutismBoth, paste0(outDir, "/profoundAutismBoth_above8.csv"))
verbalMildID <- verbalMildID[which(verbalMildID$age_at_ados / 12 > 8),]
write.csv(verbalMildID, paste0(outDir, "/verbalMildID_above8.csv"))
verbalNoID <- verbalNoID[which(verbalNoID$age_at_ados / 12 > 8),]
write.csv(verbalNoID, paste0(outDir, "/verbalNoID_above8.csv"))
verbalGifted <- verbalGifted[which(verbalGifted$age_at_ados / 12 > 8),]
write.csv(verbalGifted, paste0(outDir, "/verbalGifted_above8.csv"))

# Breakdown by age, sex, race, and ethnicity.
listToIterate <- list(rownames(profoundAutismModerateIDOnly), 
                      rownames(profoundAutismNonverbalOnly), 
                      rownames(profoundAutismBoth),
                      rownames(verbalMildID), 
                      rownames(verbalNoID), 
                      rownames(verbalGifted))
sexBreakdown <- do.call(rbind, lapply(listToIterate, function(group){
  sex <- data.frame(male = length(intersect(group, rownames(ssc)[which(ssc$sex == "male")])),
                    female = length(intersect(group, rownames(ssc)[which(ssc$sex == "female")])))
}))
ageBreakdown <- do.call(rbind, lapply(listToIterate, function(group){
  age <- data.frame(ageMean = mean(ssc[group, "age_at_ados"] / 12),
                    ageSD = sd(ssc[group, "age_at_ados"] / 12))
}))
raceBreakdown <- do.call(rbind, lapply(listToIterate, function(group){
  race <- data.frame(africanAmerican = length(intersect(group, rownames(ssc)[which(ssc$race == "african-amer")])),
                     asian = length(intersect(group, rownames(ssc)[which(ssc$race == "asian")])),
                     multiracial = length(intersect(group, rownames(ssc)[which(ssc$race == "more-than-one-race")])),
                     nativeAmerican = length(intersect(group, rownames(ssc)[which(ssc$race == "native-american")])),
                     nativeHawaiian = length(intersect(group, rownames(ssc)[which(ssc$race == "native-hawaiian")])),
                     notSpecified = length(intersect(group, rownames(ssc)[which(ssc$race == "not-specified")])),
                     other = length(intersect(group, rownames(ssc)[which(ssc$race == "other")])),
                     white = length(intersect(group, rownames(ssc)[which(ssc$race == "white")])))
}))
ethnicityBreakdown <- do.call(rbind, lapply(listToIterate, function(group){
  ethnicity <- data.frame(hispanic = length(intersect(group, rownames(ssc)[which(ssc$ethnicity == "hispanic")])),
                          nonhispanic = length(intersect(group, rownames(ssc)[which(ssc$ethnicity == "non-hispanic")])))
}))

# Make counts into table.
counts <- data.frame(iqRange = c("0-49", "50-160", "0-49", "50-69", "70-114", "115-160"),
                     expressiveLanguage = c("phrase speech", "no phrase speech", "no phrase speech",
                                            rep("phrase speech", 3)))
allData <- do.call(cbind, list(counts, sexBreakdown, ageBreakdown, raceBreakdown, ethnicityBreakdown))
write.csv(allData, paste0(outDir, "/statSummaryAbove8.csv"), row.names = FALSE)

# Test for statistical significance.
n <- allData$male + allData$female
cutoff <- 0.05
allFacStats <- do.call(rbind, lapply(c("africanAmerican", "asian", "multiracial", "white", "hispanic", "nonhispanic",
                "male", "female"), function(factor){
  facStats <- do.call(rbind, lapply(1:6, function(i){
    dat <- data.frame(stat_yes = c(allData[i, factor],
                                   sum(allData[,factor]) - allData[i, factor]),
                      stat_no = c(n[i] - allData[i, factor], 
                                  (sum(n) - n[i]) - (sum(allData[,factor]) - allData[i, factor])),
                      row.names = c(paste0("yes_", allData[i, "iqRange"]), 
                                    paste0("no_", allData[i, "iqRange"])))
    colnames(dat) <- c(paste0("yes_", factor), paste0("no_", factor))
    fishLess <- stats::fisher.test(dat, alternative = "less")
    fishMore <- stats::fisher.test(dat, alternative = "greater")
    return(data.frame(fishersLess = fishLess$p.value, fishersMore = fishMore$p.value,
                      factor = factor, iqRange = allData[i, "iqRange"],
                      expressiveLanguage = allData[i, "expressiveLanguage"]))
  }))
  return(facStats)
}))
write.csv(allFacStats, paste0(outDir, "/fishersTestsAbove8.csv"), row.names = FALSE)
