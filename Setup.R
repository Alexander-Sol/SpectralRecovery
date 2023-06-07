#Single Cell Data Analysis

library(tidyverse)
library(viridis)
library(magrittr)

# Path to Local Folders ----
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

# Vectors used to determine conditions based on file names ----
M24Conditions <- c("2ng", "0.2ng")
M24ConditionPatterns <- c("_2ng", "_02ng")
names(M24ConditionPatterns) <- M24Conditions

P17Conditions <- c("5ng", "2.5ng", "1ng", "500pg", "250pg")
P17ConditionPatterns <- c("_250", "_500pg", "_1ng",  "_2c5ng", "_5ng") %>% rev()
names(P17ConditionPatterns) <- P17Conditions

P55Conditions <- c("500pg", "250pg", "125pg")
P55ConditionPatterns <- c("neu0_5", "neu0_25", "neu0_125")
names(P55ConditionPatterns) <- P55Conditions

M37_HelaConditions <- c("Library", "Single Cell")
M37_HelaConditionPatterns <- c("Lib_", "SC_")
names(M37_HelaConditionPatterns) <- M37_HelaConditions

ConditionPatterns <- list(
  M37_HeLa = M37_HelaConditionPatterns,
  M37_Lung = NULL,
  M24 = rev(M24ConditionPatterns),
  P15 = NULL,
  P17 = rev(P17ConditionPatterns),
  P55 = rev(P55ConditionPatterns)
)

singleCellConditions <- c("0.2ng", "500pg", "250pg", "125pg", "Single Cell")
# Combined File Analysis ----

SummarizeAllQuantifiedPeaks <- function(filePath, scoreThreshold) 
{
  read_tsv(filePath) %>%
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
    )
}

SummarizeRecoveredSpectra <- function(filePath, scoreThreshold)
{
  read_tsv(filePath) %>%
    group_by(`File Name`) %>%
    summarise(
      SpectraRecovered = sum(`Normalized Spectral Angle` != "Spectrum Not Found") - sum(`Normalized Spectral Angle` == "-1"),
      ZeroScoreSpectra = sum(`Normalized Spectral Angle` == "0"),
      HighScoreSpectra = `Normalized Spectral Angle` [`Normalized Spectral Angle` != "Spectrum Not Found"] %>% 
        as.numeric() %>% .[. > scoreThreshold] %>% length(),
      PercentZeroSpectra = 100 * ZeroScoreSpectra / SpectraRecovered,
    )
}

SummarizeAllPsms <- function(filePath)
{
  read_tsv(filePath) %>%
    filter(Contaminant == "N" & Decoy == "N") %>% 
    filter(as.numeric(QValue) <= 0.01) %>%
    group_by(`File Name`) %>%
      summarise(
        UniquePeptides = unique(`Full Sequence`) %>% length()
      )
}

SummarizeMqEvidence <- function(filePath)
{
  maxQuantMbrCount <- read_tsv(filePath) %>%
    group_by(`Raw file`) %>%
    summarise(MbrPeaks = sum(!is.na(`Match m/z difference`)))
  
  maxQuantPeptideCount <- read_tsv(filePath) %>%
    filter(is.na(Reverse) & is.na(`Potential contaminant`)) %>%
    filter(`MS/MS count` > 0) %>%
    filter(PEP <= 0.01) %>%
    group_by(`Raw file`) %>%
    summarise(
      UniquePeptides = unique(`Modified sequence`) %>% length()
    )
  
  return(merge(maxQuantMbrCount, maxQuantPeptideCount, by = "Raw file"))
}

MergeMqResults <- function(folder)
{
  evidenceResults <- SummarizeMqEvidence(paste0(folder, "/evidence.txt"))
  summaryFile <- read_tsv(paste0(folder, "/summary.txt"))
  
  mergedResults <- merge(evidenceResults, summaryFile, by = "Raw file") %>%
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
      .keep = "none")
    
  return(mergedResults)
}

#Analyzes multiple output files, then combines the information into one table
RunAnalysis <- function(MetaFolder, MaxQuantFolder, HighScoreThreshold = 0.6) {
  
  # MetaMorpheus Analysis
  mbrTableMeta <- 
    SummarizeAllQuantifiedPeaks(paste0(MetaFolder, "/AllQuantifiedPeaks.tsv"), HighScoreThreshold)
  metaSummary <- SummarizeAllPsms(paste0(MetaFolder, "/AllPSMs.psmtsv" ))
  metaResults <- merge(mbrTableMeta, metaSummary, by = "File Name")
  
  # MaxQuant Analysis
  maxQuantSummary <- MergeMqResults(MaxQuantFolder)
  maxQuantRecoveredSpectra <- 
    SummarizeRecoveredSpectra(paste0(MaxQuantFolder, "/RecoveredSpectra.psmtsv"), HighScoreThreshold)
  
  # Merge results from MaxQuant evidence and summary files 
  maxQuantResults <- 
    merge(maxQuantSummary, maxQuantRecoveredSpectra, by = "File Name") %>%
    mutate(PercentRecovery = 100 * SpectraRecovered/MbrPeaks)

  # Removing the -calib apppended to file names to match the maxQuant results
  metaResults$`File Name` <- str_replace(metaResults$`File Name`, "-calib", "")
  
  metaResults$SearchEngine <- "MetaMorpheus"
  maxQuantResults$SearchEngine <- "MaxQuant"
  
  # Merge the results of each search engine
  return( bind_rows(metaResults, maxQuantResults) )
}

# Function used to write the condition (e.g., sample mass) for each file in
# a given dataset
AssignConditionCombined <- function(fileNames, conditionPatterns){
  
  conditionAssignment <- character(length(fileNames))
  
  if(is.null(conditionPatterns))
  {
    conditionAssignment <- rep("Single Cell", length(conditionAssignment))
  } 
  else
  {
    for(i in 1:length(conditionPatterns))
    {
      conditionAssignment[grep(conditionPatterns[i], fileNames)] <- names(conditionPatterns)[i]
    }
  }
  
  return(conditionAssignment)
}

GetAnalysisResults <- function(singleCellEquivalent = TRUE)
{
  
  for (i in 1:6)
  {
    singleAnalysis <- RunAnalysis(MetaResultsFolders[[i]], MaxQuantResultsFolders[[i]]) %>%
      mutate(
        Condition = AssignConditionCombined(`File Name`, ConditionPatterns[[i]]),
        Dataset = names(MetaResultsFolders)[i]
        )
    
    if(singleCellEquivalent)
    {
      singleAnalysis %<>% filter(Condition %in% singleCellConditions)
    }

    if (i == 1) 
    {
      combinedAnalysis <- singleAnalysis
    } 
    else 
    {
      combinedAnalysis %<>% add_row(singleAnalysis)
    }
  }
  
  combinedAnalysis$Dataset <- factor(combinedAnalysis$Dataset, 
                                     c("M37_HeLa", "M37_Lung", "P15", "M24", "P17", "P55") )
  return(combinedAnalysis)
}

GetSummaryResults <- function(singleCellEquivalent = TRUE)
{
  for (i in 1:6){
    singleAnalysis <- RunAnalysis(MetaResultsFolders[[i]], MaxQuantResultsFolders[[i]]) %>%
      mutate(Condition = AssignConditionCombined(`File Name`, ConditionPatterns[[i]])) %>%
      filter(Condition %in% singleCellConditions) %>% # Need to actually make this a check
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
  return(combinedAnalysis)
}

# Individual file analysis functions ----
GetZeroPeaks <- function(peakQuant, threshold = 0.45) 
{
  zeroAnglePeaks <- peakQuant %>%
    mutate(Angle  = as.numeric(`Normalized Spectral Angle`)) %>%
    filter(Angle <= threshold & Angle > -1)
  return(zeroAnglePeaks)
}

GetHighPeaks <- function(peakQuant, threshold = 0.54) 
{
  highAnglePeaks <- peakQuant %>%
    mutate(Angle  = as.numeric(`Normalized Spectral Angle`)) %>%
    filter(Angle > threshold)
  return(highAnglePeaks)
}












