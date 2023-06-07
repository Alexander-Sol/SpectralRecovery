# Percent Recovered Spectra

summaryResults <- GetSummaryResults()

ggplot(summaryResults, aes(x = Dataset,
                             y = PercentRecoveryMean,
                             fill = SearchEngine)) +
  geom_bar(
    stat="identity",
    position = position_dodge(),
    aes(color = SearchEngine),
    lwd = 0.7) +
  geom_errorbar(
    aes(ymin = PercentRecoveryMean-PercentRecoverySd, ymax = PercentRecoveryMean+PercentRecoverySd),
    width = 0.4, color = "darkorange", size = 1, position = position_dodge(.9)) +
  ylab("Percentage of Transferred IDs with Missed Spectra") + 
  xlab("Dataset") + 
  scale_color_manual(
    values = c(
      "MaxQuant" = "darkblue", 
      "MetaMorpheus" = "darkred")) + 
  scale_fill_manual(
    values = c(
      "MaxQuant" = adjustcolor("darkblue", alpha.f = 0.85),
      "MetaMorpheus" = adjustcolor("darkred", alpha.f = 0.85))) + 
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
    "M37_Lung" = "Lung\n2022, PNNL",
    "P15" = "HeLa + Neurons\n2021, Kelly" ,
    "M24" = "Hela\n2021, Kelly",
    "P17" = "HeLa\n2021,Mechtler",
    "P55" = "Neurons\n2022, Nemes"))


