# High-Low Comparison

MaxQuantPeaks <- read_tsv(paste0(MaxQuantResultsFolders$M37_HeLa, "/RecoveredSpectra.psmtsv")) %>%
  mutate(Condition = AssignConditionCombined(`File Name`, ConditionPatterns[[i]])) %>%
  filter(Condition %in% singleCellConditions)
MaxQuantPeaks$ApexMs2Offset <- MaxQuantPeaks$`Peak Apex RT (min)` - MaxQuantPeaks$`Scan Retention Time`

LowAnglePeaks <- GetZeroPeaks(MaxQuantPeaks, threshold = 0.05) %>% 
  mutate(AngleGroup = "Low")
HighAnglePeaks <- GetHighPeaks(MaxQuantPeaks) %>% 
  mutate(AngleGroup = "High")

classifiedPeaks <- bind_rows(LowAnglePeaks, HighAnglePeaks)

ggplot(classifiedPeaks, aes(`Precursor Isotopic Envelope Score`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2)  + ggtitle("Precursor Isotopic Envelope Score")

ggplot(classifiedPeaks, aes(`Mass Diff (ppm)`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + xlim(c(-50,50)) + ggtitle("Mass Diff (ppm)")

ggplot(classifiedPeaks, aes(`Peak Width (min)`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + xlim(c(0,2)) + ggtitle("Peak Width (Minutes)")

ggplot(classifiedPeaks, aes(`RT Shift`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + ggtitle("Retention Time Shift")

ggplot(classifiedPeaks, aes(`Precursor m/z - Isolation Center Distance (Th)`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + xlim(-1, 1) + ggtitle("Isolation Window Center - Precursor m/z Distance")

ggplot(classifiedPeaks, aes(`ApexMs2Offset`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + ggtitle("Time between peak apex and MS2 collection")
