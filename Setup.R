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

ConditionPatterns <- list(
  M37_HeLa = NULL,
  M37_Lung = NULL,
  M24 = rev(M24ConditionPatterns),
  P15 = NULL,
  P17 = rev(P17ConditionPatterns),
  P55 = rev(P55ConditionPatterns)
)

singleCellConditions <- c("0.2ng", "500pg", "250pg", "125pg", "Single Cell")