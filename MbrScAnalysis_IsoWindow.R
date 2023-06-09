# Single Cell Data Analysis
# V2, used for IsoWindow Spectral Recover Analysis
library(tidyverse)
library(viridis)
library(magrittr)

MetaResultsFolders <- list(
  M37_HeLa = "D:/SingleCellDataSets/MetaMorpheusResults/MSV000085937_HeLa",
  M37_Lung = "D:/SingleCellDataSets/MetaMorpheusResults/MSV000085937_Lung",
  M24 = "D:/SingleCellDataSets/MetaMorpheusResults/MSV000087524",
  P15 = "D:/SingleCellDataSets/MetaMorpheusResults/PXD019515",
  P17 = "D:/SingleCellDataSets/MetaMorpheusResults/PXD024017",
  P55 = "D:/SingleCellDataSets/MetaMorpheusResults/PXD031955"
)

MaxQuantResultsFolders <- list(
  "D:/SingleCellDataSets/MaxQuantResults/MSV000085937_HeLa",
  "D:/SingleCellDataSets/MaxQuantResults/MSV000085937_Lung",
  "D:/SingleCellDataSets/MaxQuantResults/MSV000087524",
  "D:/SingleCellDataSets/MaxQuantResults/PXD019515",
  "D:/SingleCellDataSets/MaxQuantResults/PXD024017",
  "D:/SingleCellDataSets/MaxQuantResults/PXD031955"
)

names(MaxQuantResultsFolders) <- names(MetaResultsFolders)

# Condition Patterns ----
M24Conditions <- c("2ng", "0.2ng")
M24ConditionPatterns <- c("_2ng", "_02ng")
names(M24ConditionPatterns) <- M24Conditions

P17Conditions <- c("5ng", "2.5ng", "1ng", "500pg", "250pg")
P17ConditionPatterns <- c("_250", "_500pg", "_1ng",  "_2c5ng", "_5ng") %>% rev()
names(P17ConditionPatterns) <- P17Conditions


P55Conditions <- c("500pg", "250pg", "125pg")
P55ConditionPatterns <- c("neu0_5", "neu0_25", "neu0_125")
names(P55ConditionPatterns) <- P55Conditions

ConditionPatterns <- list(
  M37_HeLa = NULL,
  M37_Lung = NULL,
  M24 = rev(M24ConditionPatterns),
  P15 = NULL,
  P17 = rev(P17ConditionPatterns),
  P55 = rev(P55ConditionPatterns)
)

singleCellConditions <- c("0.2ng", "500pg", "250pg", "125pg", "Single Cell")

# Constants
scoreThreshold = 0.6

# File Analysis ----
RunAnalysis <- function(MetaFolder, MaxQuantFolder) {
  
  # MetaMorpheus Analysis ----
  
  mbrTableMeta <- read_tsv(paste0(MetaFolder, "/AllQuantifiedPeaks.tsv" )) %>%
    group_by(`File Name`) %>%
    summarise(
      TotalPeaks = n(),
      MsmsPeaks = sum(`Peak Detection Type` == "MSMS"),
      MbrPeaks = sum(`Peak Detection Type` == "MBR"),
      SpectraRecovered = sum(!(`Spectral Contrast Angle` %in% c("Spectrum Not Found", "-1"))),
      ZeroScoreSpectra = sum(`Spectral Contrast Angle` == "0"),
      HighScoreSpectra = `Spectral Contrast Angle` [`Spectral Contrast Angle` != "Spectrum Not Found"] %>% 
        as.numeric() %>% .[. > scoreThreshold] %>% length(),
      PercentMbr = 100 * MbrPeaks / TotalPeaks,
      PercentRecovery = 100 * SpectraRecovered / MbrPeaks,
      PercentZeroSpectra = 100 * ZeroScoreSpectra / SpectraRecovered,
      # AvgSpectralAngle = `Spectral Contrast Angle`[`Spectral Contrast Angle` != "Spectrum Not Found"] %>% 
      #   as.numeric() %>% mean(),
      # AvgAbsRtShift_Min = `Retention Time Shift (min)`[!is.na(`Retention Time Shift (min)`)] %>% 
      #   as.numeric() %>% abs() %>% mean(),
      # AvgAbsRtShift_Zscore = `Retention Time Z-Score`[!is.na(`Retention Time Z-Score`)] %>% 
      #   as.numeric() %>% abs() %>% mean(),
      # AvgAbsPpmError = `Peak Apex Mass Error (ppm)`[`Spectral Contrast Angle` != "Spectrum Not Found" & `Peak Apex Mass Error (ppm)` != "NaN"] %>% 
      #   as.numeric() %>% abs() %>% mean()
      )
  
  metaSummary <- read_tsv(paste0(MetaFolder, "/AllPSMs.psmtsv" )) %>%
    filter(Contaminant == "N" & Decoy == "N") %>%
    filter(as.numeric(QValue) <= 0.01) %>%
    group_by(`File Name`) %>%
    summarise(
      UniquePeptides = unique(`Full Sequence`) %>% length()
    )
  
  metaResults <- merge(mbrTableMeta, metaSummary, by = "File Name")
  metaResults$`File Name` <- str_replace(metaResults$`File Name`, "-calib", "")
    
  
  # MaxQuant Analysis ----
  
  maxQuantMbrCount <- read_tsv(paste0(MaxQuantFolder, "/evidence.txt" )) %>%
    group_by(`Raw file`) %>%
    summarise(MbrPeaks = sum(!is.na(`Match m/z difference`)))
  
  maxQuantPeptideCount <- read_tsv(paste0(MaxQuantFolder, "/evidence.txt" )) %>%
    filter(is.na(Reverse) & is.na(`Potential contaminant`)) %>%
    filter(`MS/MS count` > 0) %>%
    filter(PEP <= 0.01) %>%
    group_by(`Raw file`) %>%
    summarise(
      UniquePeptides = unique(`Modified sequence`) %>% length()
    )
  
  maxQuantCount <- merge(maxQuantMbrCount, maxQuantPeptideCount, by = "Raw file")
  
  maxQuantSummary <- merge(maxQuantCount,
                           read_tsv(paste0(MaxQuantFolder, "/summary.txt" )), 
                           by = "Raw file") %>%
    mutate(
      `File Name` = `Raw file`,
      MbrPeaks = MbrPeaks,
      MsmsPeaks = `MS/MS identified`,
      PercentMbr = 100 * MbrPeaks / (MsmsPeaks + MbrPeaks),
      UniquePeptides = `Peptide sequences identified`,
      MsmsScans = `MS/MS`,
      TotalPeaks = Peaks,
      FragmentedPeaks = `Peaks sequenced`,
      IsotopeEnvelopes = `Isotope patterns`,
      FragmentedIsotopeEnvelopes = `Isotope patterns sequenced`,
      PercentIsotopeEnvelopesFragmented = 100 * FragmentedIsotopeEnvelopes / IsotopeEnvelopes,
      .keep = "none"
    )
  
  mbrTableMaxQuant <- read_tsv(paste0(MaxQuantFolder, "/RecoveredSpectra.psmtsv"  )) %>%
    group_by(`File Name`) %>%
    summarise(
      #MbrPeaksRead = sum(`Peak Detection Type` == "MBR"), # Should be equal to n()
      SpectraRecovered = sum(`Normalized Spectral Angle` != "Spectrum Not Found") - sum(`Normalized Spectral Angle` == "-1"),
      ZeroScoreSpectra = sum(`Normalized Spectral Angle` == "0"),
      HighScoreSpectra = `Normalized Spectral Angle` [`Normalized Spectral Angle` != "Spectrum Not Found"] %>% 
        as.numeric() %>% .[. > scoreThreshold] %>% length(),
      PercentZeroSpectra = 100 * ZeroScoreSpectra / SpectraRecovered,
      AvgSpectralAngle = `Normalized Spectral Angle`[`Normalized Spectral Angle` != "Spectrum Not Found"] %>% 
        as.numeric() %>% mean(),
      AvgAbsRtShift_Min = `RT Shift`[!is.na(`RT Shift`)] %>% 
        as.numeric() %>% abs() %>% mean(),
      AvgAbsPpmError = `Mass Diff (ppm)`[`Normalized Spectral Angle` != "Spectrum Not Found"] %>% 
        as.numeric() %>% abs() %>% mean()
    )
  
  maxQuantResults <- merge(mbrTableMaxQuant, maxQuantSummary, by = "File Name") %>%
    mutate(PercentRecovery = 100 * SpectraRecovered/MbrPeaks)
  
  metaResults$SearchEngine <- "MetaMorpheus"
  maxQuantResults$SearchEngine <- "MaxQuant"
  
  return( bind_rows(metaResults, maxQuantResults) )
  
}


# Combined File Analysis ----

AssignConditionCombined <- function(fileNames, conditionPatterns){
  conditionAssignment <- character(length(fileNames))
  if(is.null(conditionPatterns)){
    conditionAssignment <- rep("Single Cell", length(conditionAssignment))
  } else{
    for(i in 1:length(conditionPatterns))
    {
      conditionAssignment[grep(conditionPatterns[i], fileNames)] <- names(conditionPatterns)[i]
    }
  }
  return(conditionAssignment)
}


for (i in 1:6){
  singleAnalysis <- RunAnalysis(MetaResultsFolders[[i]], MaxQuantResultsFolders[[i]]) %>%
    mutate(Condition = AssignConditionCombined(`File Name`, ConditionPatterns[[i]]),
           Dataset = names(MetaResultsFolders)[i]) %>%
    filter(Condition %in% singleCellConditions)
  if (i == 1) {
    combinedAnalysis <- singleAnalysis
  } else {
    combinedAnalysis %<>% add_row(singleAnalysis)
  }
}
combinedAnalysis$Dataset <- factor(combinedAnalysis$Dataset, 
                                   c("M37_HeLa", "M37_Lung", "P15", "M24", "P17", "P55") )

test <- combinedAnalysis %>%
  group_by(Dataset, SearchEngine) %>%
  summarise(
    PercentRecovery = mean(PercentRecovery),
    PercentZeroSpectra = mean(PercentZeroSpectra),
    SpectraRecovered = sum(SpectraRecovered)
  )
test2 <- data.frame(
  x_label = c( rep("MBR Matches with\nCorresponding Spectra", 12), rep("Recovered Spectra with\nSpectral Angle of Zero", 12) ),
  values = c(test$PercentRecovery, test$PercentZeroSpectra),
  Dataset = rep(test$Dataset,2),
  SearchEngine = rep(test$SearchEngine, 2)
)

p <- ggplot(test2, aes(x = x_label,
                  y = values)) +
  geom_boxplot(lwd = 0.7,
               aes(color = SearchEngine),
               fill = "gray90",
               outlier.shape = NA) + 
  ylab("Fraction of Events") + 
  xlab("") +
  scale_color_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred")) +
  geom_point(position=position_dodge(width = 0.75),
    #dodge.width = 0.75, jitter.width = 0.35),
             aes(fill = Dataset, group = SearchEngine),
             shape = 21,
             size = 2.25,
             stroke = 0.05,
             color = "black") +
  scale_fill_manual( values = c(
    "M37_HeLa" = "black",
    "M37_Lung" = "#440154FF",
    "P15" = "#443A83FF" ,
    "M24" = "#31688EFF",
    "P17" = "#21908CFF",
    "P55" = "#35B779FF" ),
    labels = c(
      "M37_HeLa" = "HeLa, 2022, PNNL",
      "M37_Lung" = "HeLa, 2022, PNNL",
      "P15" = "HeLa + Neurons, 2021, Kelly" ,
      "M24" = "Hela, 2021, Kelly",
      "P17" = "HeLa, 2021, Mechtler",
      "P55" = "Neurons, 2022, Nemes")
  ) +
  theme_classic() + 
  theme(panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="darkgray" ),
        panel.grid.minor.y = element_line( size=.1, color="gray" )) +
  # Force the x axis to 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80))


png("EventFractions.png", units = "in", width = 6, height = 5, res = 300)
p
dev.off() 

d <- ggplot(combinedAnalysis, aes(x = Dataset,
                             y = PercentMbr)) +
  geom_boxplot(
    lwd = 0.7,
     aes(color = SearchEngine),
     fill = "gray90",
     outlier.shape = NA) +
  ylab("Percentage of Identifications Made via MBR") + 
  xlab("Dataset") + 
  scale_color_manual(
    values = c(
      "MaxQuant" = "darkblue", 
      "MetaMorpheus" = "darkred")) + 
  geom_point(
    position = position_dodge(width = 0.75),
     aes(x = Dataset, group = SearchEngine),
     shape = 21,
     size = 2,
     stroke = 0.05,
     color = "black",
     fill = "black") +
  theme(legend.position = "bottom") +
  theme_classic() + 
  theme(panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="darkgray" ),
        panel.grid.minor.y = element_line( size=.1, color="gray" ),
        axis.text.x = element_text(angle = 0)) +
  # Force the x axis to 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_x_discrete(labels = c(
    "M37_HeLa" = "HeLa\n2022, PNNL",
    "M37_Lung" = "HeLa\n2022, PNNL",
    "P15" = "HeLa + Neurons\n2021, Kelly" ,
    "M24" = "Hela\n2021, Kelly",
    "P17" = "HeLa\n2021,Mechtler",
    "P55" = "Neurons\n2022, Nemes"))

png("PercentMBR_boxplot.png", units = "in", width = 8, height = 6, res = 300)
d
dev.off()

ggplot(combinedAnalysis, aes(x = Dataset,
                              y = PercentMbr)) +
  geom_col(
    stat="identity",
    position = "dodge",
    aes(color = SearchEngine, fill = SearchEngine),
    lwd = 1) +
  ylab("Percentage of Identifications Made via MBR") + 
  xlab("Dataset") + 
  scale_color_manual(
    values = c(
      "MaxQuant" = "darkblue", 
      "MetaMorpheus" = "darkred")) + 
  scale_fill_manual(
    values = c(
      "MaxQuant" = adjustcolor("darkblue", alpha.f = 0.75),
      "MetaMorpheus" = adjustcolor("darkred", alpha.f = 0.75))) + 
  geom_point(
    position = position_dodge(width = 0.75),
    aes(x = Dataset, group = SearchEngine),
    shape = 21,
    size = 2,
    stroke = 0.05,
    color = "black",
    fill = "black") +
  theme(legend.position = "bottom") +
  theme_classic() + 
  theme(panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="darkgray" ),
        panel.grid.minor.y = element_line( size=.1, color="gray" ),
        axis.text.x = element_text(angle = 0)) +
  # Force the x axis to 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_x_discrete(labels = c(
    "M37_HeLa" = "HeLa\n2022, PNNL",
    "M37_Lung" = "HeLa\n2022, PNNL",
    "P15" = "HeLa + Neurons\n2021, Kelly" ,
    "M24" = "Hela\n2021, Kelly",
    "P17" = "HeLa\n2021,Mechtler",
    "P55" = "Neurons\n2022, Nemes"))


# Summary Plotting ----

for (i in 1:6){
  singleAnalysis <- RunAnalysis(MetaResultsFolders[[i]], MaxQuantResultsFolders[[i]]) %>%
    mutate(Condition = AssignConditionCombined(`File Name`, ConditionPatterns[[i]])) %>%
    filter(Condition %in% singleCellConditions) %>%
    group_by(SearchEngine) %>%
    summarise(
      `File Name` = names(MetaResultsFolders)[i],
      MbrPeaksMean = mean(MbrPeaks),
      MbrPeaksSd = sd(MbrPeaks),
      MsmsIdsMean = mean(MsmsPeaks),
      MsmsIdsSd = sd(MsmsPeaks),
      PercentMbrMean = mean(PercentMbr),
      PercentMbrSd = sd(PercentMbr),
      SpectraRecoveredMean = mean(SpectraRecovered),
      SpectraRecoveredSd = sd(SpectraRecovered),
      PercentRecoveryMean = mean(PercentRecovery),
      PercentRecoverySd = sd(PercentRecovery),
      UniquePeptidesMean = mean(UniquePeptides),
      UniquePeptidesSd = sd(UniquePeptides),
      IsotopeEnvelopesMean = mean(IsotopeEnvelopes),
      IsotopeEnvelopesSd = sd(IsotopeEnvelopes),
      FragmentedIsotopeEnvelopesMean = mean(FragmentedIsotopeEnvelopes),
      FragmentedIsotopeEnvelopesSd = sd(FragmentedIsotopeEnvelopes),
      PercentIsotopeEnvelopesFragmentedMean = mean(PercentIsotopeEnvelopesFragmented),
      PercentIsotopeEnvelopesFragmentedSd = sd(PercentIsotopeEnvelopesFragmented)
    )
  if (i == 1) {
    combinedAnalysis <- singleAnalysis
  } else {
    combinedAnalysis %<>% add_row(singleAnalysis)
  }
}
combinedAnalysis$Dataset <- factor(combinedAnalysis$`File Name`, 
                                   c("M37_HeLa", "M37_Lung", "P15", "M24", "P17", "P55") )


d <- ggplot(combinedAnalysis, aes(x = Dataset,
                             y = PercentMbrMean,
                             fill = SearchEngine)) +
  geom_bar(
    stat="identity",
    position = position_dodge(),
    aes(color = SearchEngine),
    lwd = 0.7) +
  geom_errorbar(
    aes(ymin = PercentMbrMean-PercentMbrSd, ymax = PercentMbrMean+PercentMbrSd),
    width = 0.4, color = "darkorange", size = 1, position = position_dodge(.9)) +
  ylab("Percentage of Identifications Made via MBR") + 
  xlab("Dataset") + 
  scale_color_manual(
    values = c(
      "MaxQuant" = "darkblue", 
      "MetaMorpheus" = "darkred")) + 
  scale_fill_manual(
    values = c(
      "MaxQuant" = adjustcolor("darkblue", alpha.f = 0.85),
      "MetaMorpheus" = adjustcolor("darkred", alpha.f = 0.85))) + 
  # geom_point(
  #   position = position_dodge(width = 0.75),
  #   aes(x = Dataset, group = SearchEngine),
  #   shape = 21,
  #   size = 2,
  #   stroke = 0.05,
  #   color = "black",
  #   fill = "black") +
  theme_classic() + 
  theme(panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="darkgray" ),
        panel.grid.minor.y = element_line( size=.1, color="gray" ),
        axis.text.x = element_text(angle = 0),
        legend.position = "bottom") +
  # Force the x axis to 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_x_discrete(labels = c(
    "M37_HeLa" = "HeLa\n2022, PNNL",
    "M37_Lung" = "HeLa\n2022, PNNL",
    "P15" = "HeLa + Neurons\n2021, Kelly" ,
    "M24" = "Hela\n2021, Kelly",
    "P17" = "HeLa\n2021,Mechtler",
    "P55" = "Neurons\n2022, Nemes"))


d <- ggplot(combinedAnalysis, aes(x = Dataset,
                                  y = PercentMbr)) +
  geom_boxplot(
    lwd = 0.7,
    aes(color = SearchEngine),
    fill = "gray90",
    outlier.shape = NA) +
  ylab("Percentage of Identifications Made via MBR") + 
  xlab("Dataset") + 
  scale_color_manual(
    values = c(
      "MaxQuant" = "darkblue", 
      "MetaMorpheus" = "darkred")) + 
  geom_point(
    position = position_dodge(width = 0.75),
    aes(x = Dataset, group = SearchEngine),
    shape = 21,
    size = 2,
    stroke = 0.05,
    color = "black",
    fill = "black") +
  theme(legend.position = "bottom") +
  theme_classic() + 
  theme(panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="darkgray" ),
        panel.grid.minor.y = element_line( size=.1, color="gray" ),
        axis.text.x = element_text(angle = 0)) +
  # Force the x axis to 0
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  scale_x_discrete(labels = c(
    "M37_HeLa" = "HeLa\n2022, PNNL",
    "M37_Lung" = "HeLa\n2022, PNNL",
    "P15" = "HeLa + Neurons\n2021, Kelly" ,
    "M24" = "Hela\n2021, Kelly",
    "P17" = "HeLa\n2021,Mechtler",
    "P55" = "Neurons\n2022, Nemes"))

png("PercentMBR_BocPlot.png", units = "in", width = 8, height = 6, res = 300)
d
dev.off()

png("PercentMBR.png", units = "in", width = 6, height = 4.2, res = 300)
d
dev.off()

ggplot(combinedAnalysis, aes(x = Dataset, y = PercentRecoveryMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = PercentRecoveryMean-PercentRecoverySd, ymax = PercentRecoveryMean+PercentRecoverySd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Percent MBR Events with Associated Spectra") + 
  xlab("Dataset") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))

ggplot(combinedAnalysis, aes(x = Dataset, y = MbrPeaksMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = MbrPeaksMean-MbrPeaksSd, ymax = MbrPeaksMean+MbrPeaksSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Number of MBR Peaks") + 
  xlab("Dataset") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))

p <- ggplot(combinedAnalysis, aes(x = Dataset,
                                  y = PercentMbrMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = PercentMbrMean-PercentMbrSd, ymax = PercentMbrMean+PercentMbrSd),
    width = 0.4, color = "darkorange2", size = 1, position = position_dodge(.9)) +
  ylab("Percentage of Identifications Made Via MBR") + 
  xlab("Dataset") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"),
    labels = c(
      "M37_HeLa" = "HeLa, 2022, PNNL",
      "M37_Lung" = "HeLa, 2022, PNNL",
      "P15" = "HeLa + Neurons, 2021, Kelly" ,
      "M24" = "Hela, 2021, Kelly",
      "P17" = "HeLa, 2021, Mechtler",
      "P55" = "Neurons, 2022, Nemes")) + 
  theme(legend.position = "bottom")





ggplot(combinedAnalysis, aes(x = Dataset, y = MsmsIdsMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = MsmsIdsMean-MsmsIdsSd, ymax = MsmsIdsMean+MsmsIdsSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Number of MSMS IDs") + 
  xlab("Dataset") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))

ggplot(combinedAnalysis, aes(x = Dataset, y = SpectraRecoveredMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = SpectraRecoveredMean-SpectraRecoveredSd, ymax = SpectraRecoveredMean+SpectraRecoveredSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Number of RecoveredSpectra") + 
  xlab("Dataset") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))


# Plotting ----
ggplot(test, aes(fill = SearchEngine, x = `File Name`, y = UniquePeptides)) +
  geom_bar(position = "dodge", stat = "identity")

ggplot(test, aes(fill = SearchEngine, x = `File Name`, y = MbrPeaks)) +
  geom_bar(position = "dodge", stat = "identity")

ggplot(test, aes(fill = SearchEngine, x = `File Name`, y = SpectraRecovered)) +
  geom_bar(position = "dodge", stat = "identity")

ggplot(test, aes(fill = SearchEngine, x = `File Name`, y = PercentMbr)) +
  geom_bar(position = "dodge", stat = "identity")

ggplot(test, aes(fill = SearchEngine, x = `File Name`, y = PercentRecovery)) +
  geom_bar(position = "dodge", stat = "identity")

# ConditionBasedAnalysis ----

AssignCondition <- function(fileNames, conditionPatterns){
  conditionAssignment <- character(length(fileNames))
  for(i in 1:length(conditionPatterns))
  {
    conditionAssignment[grep(conditionPatterns[i], fileNames)] <- names(conditionPatterns)[i]
  }
  return(conditionAssignment)
}

AnalyzeByCondition <- function(analysisData, conditionPatterns){
  if (is.null(conditionPatterns)){
    conditionAnalysis <- analysisData %>%
      mutate(Condition = AssignCondition(`File Name`, conditionPatterns)) %>%
      group_by(SearchEngine) %>%
      summarise(
        MbrPeaksMean = mean(MbrPeaks),
        MbrPeaksSd = sd(MbrPeaks),
        MsmsIdsMean = mean(MsmsPeaks),
        MsmsIdsSd = sd(MsmsPeaks),
        PercentMbrMean = mean(PercentMbr),
        PercentMbrSd = sd(PercentMbr),
        SpectraRecoveredMean = mean(SpectraRecovered),
        SpectraRecoveredSd = sd(SpectraRecovered),
        PercentRecoveryMean = mean(PercentRecovery),
        PercentRecoverySd = sd(PercentRecovery),
        UniquePeptidesMean = mean(UniquePeptides),
        UniquePeptidesSd = sd(UniquePeptides),
        IsotopeEnvelopesMean = mean(IsotopeEnvelopes),
        IsotopeEnvelopesSd = sd(IsotopeEnvelopes),
        FragmentedIsotopeEnvelopesMean = mean(FragmentedIsotopeEnvelopes),
        FragmentedIsotopeEnvelopesSd = sd(FragmentedIsotopeEnvelopes),
        PercentIsotopeEnvelopesFragmentedMean = mean(PercentIsotopeEnvelopesFragmented),
        PercentIsotopeEnvelopesFragmentedSd = sd(PercentIsotopeEnvelopesFragmented)
      )
    return(conditionAnalysis)
  }
  conditionAnalysis <- analysisData %>%
    mutate(Condition = AssignCondition(`File Name`, conditionPatterns)) %>%
    group_by(Condition, SearchEngine) %>%
    summarise(
      MbrPeaksMean = mean(MbrPeaks),
      MbrPeaksSd = sd(MbrPeaks),
      MsmsIdsMean = mean(MsmsPeaks),
      MsmsIdsSd = sd(MsmsPeaks),
      PercentMbrMean = mean(PercentMbr),
      PercentMbrSd = sd(PercentMbr),
      SpectraRecoveredMean = mean(SpectraRecovered),
      SpectraRecoveredSd = sd(SpectraRecovered),
      PercentRecoveryMean = mean(PercentRecovery),
      PercentRecoverySd = sd(PercentRecovery),
      UniquePeptidesMean = mean(UniquePeptides),
      UniquePeptidesSd = sd(UniquePeptides),
      IsotopeEnvelopesMean = mean(IsotopeEnvelopes),
      IsotopeEnvelopesSd = sd(IsotopeEnvelopes),
      FragmentedIsotopeEnvelopesMean = mean(FragmentedIsotopeEnvelopes),
      FragmentedIsotopeEnvelopesSd = sd(FragmentedIsotopeEnvelopes),
      PercentIsotopeEnvelopesFragmentedMean = mean(PercentIsotopeEnvelopesFragmented),
      PercentIsotopeEnvelopesFragmentedSd = sd(PercentIsotopeEnvelopesFragmented)
    )
  conditionAnalysis$Condition <- factor(conditionAnalysis$Condition, levels = names(conditionPatterns) )
  return(conditionAnalysis)
}

fileAnalysis <- RunAnalysis(MetaResultsFolders$M24, MaxQuantResultsFolders$M24)
conditionAnalysis <- AnalyzeByCondition(fileAnalysis, ConditionPatterns$M24)

fileAnalysis <- RunAnalysis(MetaResultsFolders$P17, MaxQuantResultsFolders$P17)
conditionAnalysis <- AnalyzeByCondition(fileAnalysis, ConditionPatterns$P17)

fileAnalysis <- RunAnalysis(MetaResultsFolders$P55, MaxQuantResultsFolders$P55)
conditionAnalysis <- AnalyzeByCondition(fileAnalysis, ConditionPatterns$P55)

# Condition Based Plotting ----

ggplot(conditionAnalysis, aes(x = Condition, y = MbrPeaksMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = MbrPeaksMean-MbrPeaksSd, ymax = MbrPeaksMean+MbrPeaksSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Number of MBR Peaks") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))

ggplot(conditionAnalysis, aes(x = Condition, y = MsmsIdsMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = MsmsIdsMean-MsmsIdsSd, ymax = MsmsIdsMean+MsmsIdsSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Number of MSMS Ids") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))

ggplot(conditionAnalysis, aes(x = Condition, y = PercentMbrMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = PercentMbrMean-PercentMbrSd, ymax = PercentMbrMean+PercentMbrSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))

ggplot(conditionAnalysis, aes(x = Condition, y = SpectraRecoveredMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = SpectraRecoveredMean-SpectraRecoveredSd, ymax = SpectraRecoveredMean+SpectraRecoveredSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Spectra Recovered") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))

ggplot(conditionAnalysis, aes(x = Condition, y = PercentRecoveryMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = PercentRecoveryMean-PercentRecoverySd, ymax = PercentRecoveryMean+PercentRecoverySd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Percent MBR Events with Associated Spectra") + 
  scale_fill_manual(values = c(
    "MaxQuant" = "darkblue", 
    "MetaMorpheus" = "darkred"))

# MaxQuant Stats ----
ggplot(conditionAnalysis %>% filter(SearchEngine == "MaxQuant"),
       aes(x = Condition, y = IsotopeEnvelopesMean)) +
  geom_bar(position = position_dodge(), stat = "identity", fill = "royalblue") +
  geom_errorbar(
    aes(ymin = IsotopeEnvelopesMean-IsotopeEnvelopesSd, ymax = IsotopeEnvelopesMean+IsotopeEnvelopesSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9))  +
  ylab("Isotope Envelopes Detected")

ggplot(conditionAnalysis %>% filter(SearchEngine == "MaxQuant"),
       aes(x = Condition, y = PercentIsotopeEnvelopesFragmentedMean)) +
  geom_bar(position = position_dodge(), stat = "identity", fill = "darkblue") +
  geom_errorbar(
    aes(ymin = PercentIsotopeEnvelopesFragmentedMean-PercentIsotopeEnvelopesFragmentedSd,
        ymax = PercentIsotopeEnvelopesFragmentedMean+PercentIsotopeEnvelopesFragmentedSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  ylab("Percent of Isotope Envelopes Fragmented")

ggplot(conditionAnalysis %>% filter(SearchEngine == "MaxQuant"),
       aes(x = Condition, y = MbrPeaksMean, fill = SearchEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = MbrPeaksMean-MbrPeaksSd, ymax = MbrPeaksMean+MbrPeaksSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9))


# Spectral Recovery Stats Functions ----

GetZeroPeaks <- function(peakQuant) {
  zeroAnglePeaks <- peakQuant %>%
    filter(`Spectral Contrast Angle` == "0")
  return(zeroAnglePeaks)
}

GetHighPeaks <- function(peakQuant, highAngleThreshold = 0.6) {
  highAnglePeaks <- peakQuant %>%
    filter(`Peak Detection Type` == "MBR") %>% 
    filter(`Spectral Contrast Angle` != "Spectrum Not Found") %>%
    mutate(Angle  = as.numeric(`Spectral Contrast Angle`)) %>%
    filter(Angle > highAngleThreshold)
  return(highAnglePeaks)
}

std.error <- function(data, na.rm = TRUE) {
  return(sd(data, na.rm=na.rm)/sum(!is.na(data)))
}

SummarizePeaksMeta <- function(peakTable, conditionPattern){
  peakTable$`File Name` <- str_replace(peakTable$`File Name`, "-calib", "")
  peakTable$Condition <- AssignCondition(peakTable$`File Name`, conditionPattern)
  peakSummary <- peakTable %>%
    group_by(Condition) %>%
    summarise(
      AbsPpmErrorMean = mean(abs(`Peak Apex Mass Error (ppm)`), na.rm = TRUE),
      AbsPpmErrorSd = std.error(abs(`Peak Apex Mass Error (ppm)`), na.rm = TRUE),
      IntensityMean = mean(`Peak intensity`),
      IntensitySd = std.error(`Peak intensity`),
      AbsRtShiftMean = mean(abs(`Retention Time Shift (min)`), na.rm = TRUE),
      AbsRtShiftSd = std.error(abs(`Retention Time Shift (min)`), na.rm = TRUE),
      AbsRtZScoreMean = mean(abs(`Retention Time Z-Score`), na.rm = TRUE),
      AbsRtZScoreSd = std.error(abs(`Retention Time Z-Score`), na.rm = TRUE),
    )
  peakSummary$Condition <- factor(peakSummary$Condition, levels = names(conditionPattern) )
  return(peakSummary)
}

SummarizePeaksMaxQuant <- function(peakTable, conditionPattern){
  peakTable$Condition <- AssignCondition(peakTable$`File Name`, conditionPattern)
  peakSummary <- peakTable %>%
    group_by(Condition) %>%
    summarise(
      AbsPpmErrorMean = mean(abs(`Ppm Error`), na.rm = TRUE),
      AbsPpmErrorSd = std.error(abs(`Ppm Error`), na.rm = TRUE),
      IntensityMean = mean(`Peak intensity`),
      IntensitySd = std.error(`Peak intensity`),
      AbsRtShiftMean = mean(abs(`Retention Time Shift (min)`), na.rm = TRUE),
      AbsRtShiftSd = std.error(abs(`Retention Time Shift (min)`), na.rm = TRUE)
    )
  peakSummary$Condition <- factor(peakSummary$Condition, levels = names(conditionPattern) )
  return(peakSummary)
}

# Spectral Recovery Analysis ---- 

# QuantifiedPeaks
MetaPeaks <- read_tsv(paste0(MetaResultsFolders$P17, "/AllQuantifiedPeaks.tsv"))
MaxQuantPeaks <- read_tsv(paste0(MaxQuantResultsFolders$P17, "/PeakQuant_NoArtifact.tsv"))

# Angle Groups
ZeroAngleSummary <- SummarizePeaksMeta(GetZeroPeaks(MetaPeaks), ConditionPatterns$P17) %>% 
  mutate(AngleGroup = "Zero", SearchEngine = "MetaMorpheus")
HighAngleSummary <- SummarizePeaksMeta(GetHighPeaks(MetaPeaks), ConditionPatterns$P17) %>% 
  mutate(AngleGroup = "High", SearchEngine = "MetaMorpheus")

ZeroAngleSummaryQ <- SummarizePeaksMaxQuant(GetZeroPeaks(MaxQuantPeaks), ConditionPatterns$P17) %>% 
  mutate(AngleGroup = "Zero", SearchEngine = "MaxQuant")
HighAngleSummaryQ <- SummarizePeaksMaxQuant(GetHighPeaks(MaxQuantPeaks), ConditionPatterns$P17) %>% 
  mutate(AngleGroup = "High", SearchEngine = "MaxQuant")

test <- bind_rows(ZeroAngleSummary, HighAngleSummary, ZeroAngleSummaryQ, HighAngleSummaryQ) %>%
  mutate(AngleEngine = paste0(AngleGroup, "-", SearchEngine))
test$AngleEngine <- factor(test$AngleEngine, c("Zero-MaxQuant", "High-MaxQuant", "Zero-MetaMorpheus", "High-MetaMorpheus") )


ggplot(test, aes(x = Condition, y = AbsPpmErrorMean, fill = AngleEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = AbsPpmErrorMean-AbsPpmErrorSd, ymax = AbsPpmErrorMean+AbsPpmErrorSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  scale_fill_manual(values = c(
    "Zero-MaxQuant" = "royalblue", "High-MaxQuant" = "darkblue", 
    "Zero-MetaMorpheus"= "coral2", "High-MetaMorpheus" = "darkred")) + 
  ylab("Absolute Ppm Error")

ggplot(test, aes(x = Condition, y = IntensityMean, fill = AngleEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = IntensityMean-IntensitySd, ymax = IntensityMean+IntensitySd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9)) +
  scale_fill_manual(values = c(
    "Zero-MaxQuant" = "royalblue", "High-MaxQuant" = "darkblue", 
    "Zero-MetaMorpheus"= "coral2", "High-MetaMorpheus" = "darkred")) + 
  ylab("Intensity")

ggplot(test, aes(x = Condition, y = AbsRtShiftMean, fill = AngleEngine)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(
    aes(ymin = AbsRtShiftMean-AbsRtShiftSd, ymax = AbsRtShiftMean+AbsRtShiftSd),
    width = 0.4, color = "orange", size = 1, position = position_dodge(.9))  +
  scale_fill_manual(values = c(
    "Zero-MaxQuant" = "royalblue", "High-MaxQuant" = "darkblue", 
    "Zero-MetaMorpheus"= "coral2", "High-MetaMorpheus" = "darkred")) + 
  ylab("Retention Time Shift (Min)")


t.test(abs(GetZeroPeaks(MaxQuantPeaks)$`Retention Time Shift (min)`),
       abs(GetHighPeaks(MaxQuantPeaks)$`Retention Time Shift (min)`))

t.test(abs(GetZeroPeaks(MaxQuantPeaks)$`Ppm Error`),
       abs(GetHighPeaks(MaxQuantPeaks)$`Ppm Error`))

# Spectral Recovery Analysis 2 ---- 

GetZeroPeaks <- function(peakQuant) {
  zeroAnglePeaks <- peakQuant %>%
    mutate(Angle  = as.numeric(`Normalized Spectral Angle`)) %>%
    filter(Angle < 0.42 & Angle > -1)
  return(zeroAnglePeaks)
}

GetHighPeaks <- function(peakQuant, highAngleThreshold = 0.54) {
  highAnglePeaks <- peakQuant %>%
    mutate(Angle  = as.numeric(`Normalized Spectral Angle`)) %>%
    filter(Angle > highAngleThreshold)
  return(highAnglePeaks)
}

# QuantifiedPeaks
MaxQuantPeaks <- read_tsv(paste0(MaxQuantResultsFolders$M24, "/RecoveredSpectra.psmtsv"))
MaxQuantPeaks$ApexMs2Offset <- MaxQuantPeaks$`Peak Apex RT (min)` - MaxQuantPeaks$`Scan Retention Time`

ZeroAnglePeaks <- GetZeroPeaks(MaxQuantPeaks) %>% 
  mutate(AngleGroup = "Zero")
HighAnglePeaks <- GetHighPeaks(MaxQuantPeaks) %>% 
  mutate(AngleGroup = "High")

classifiedPeaks <- bind_rows(ZeroAnglePeaks, HighAnglePeaks)


classifiedPeaks <- classifiedPeaks[classifiedPeaks$ApexMs2Offset < 0.15 & classifiedPeaks$ApexMs2Offset > -0.15,]

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



# MetaMorpheus, same dataset
MaxQuantPeaks <- read_tsv(paste0(MetaResultsFolders$M24, "/RecoveredSpectra.psmtsv"))

MaxQuantPeaks$`Normalized Spectral Angle`

ZeroAnglePeaks <- GetZeroPeaks(MaxQuantPeaks) %>% 
  mutate(AngleGroup = "Zero")
HighAnglePeaks <- GetHighPeaks(MaxQuantPeaks) %>% 
  mutate(AngleGroup = "High")

classifiedPeaks <- bind_rows(ZeroAnglePeaks, HighAnglePeaks)

classifiedPeaks$ApexMs2Offset <- classifiedPeaks$`Peak Apex RT (min)` - classifiedPeaks$`Scan Retention Time`

ggplot(classifiedPeaks, aes(`Precursor Isotopic Envelope Score`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + ggtitle("Precursor Isotopic Envelope Score")

ggplot(classifiedPeaks, aes(`Mass Diff (ppm)`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + xlim(c(-50,50)) + ggtitle("Mass Diff (ppm)")

ggplot(classifiedPeaks, aes(`Peak Width (min)`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + xlim(c(0,2))

ggplot(classifiedPeaks, aes(`RT Shift`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2)

ggplot(classifiedPeaks, aes(`Precursor m/z - Isolation Center Distance (Th)`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2) + xlim(-1, 1)

ggplot(classifiedPeaks, aes(`ApexMs2Offset`, fill = AngleGroup)) + 
  geom_density(alpha = 0.2)
