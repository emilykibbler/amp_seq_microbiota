# Read in metadata, not sure if I need this?
meta <- readRDS("meta.rds")
# Cleaned and rarified data
# All taxa
phylo_rar <- readRDS("clean_phylo_rarified.rds")
# Cleaned and rarified data
# Eliminated Streptococcaceae which is the target of the dosed antibiotics
phylo_clean_no_strep_rar <- readRDS('phylo_clean_no_strep_rar.rds') 
# Cleaned and rarified data
# Selected for only Streptococcaceae
phylo_clean_strep_rar <- readRDS('phylo_clean_strep_rar.rds') 


## Considering syndrome -------------

plot_richness(phylo_rar, 
              x = "Syndrome", #CHANGE ME, A is whatever factor you want on x-axis
              measures = "Observed", #CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title.x = element_blank()) + #example, get rid of x axis labels
  geom_violin(trim = TRUE, aes(fill = Syndrome)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(color = Syndrome)) + #optional. CHANGE ME, A is whatever factor to group box plots
  scale_color_hue(guide = "none") +
  facet_grid(.~Group, switch = "x", space = "free") + # facet_grid(.~Week, space = "free") +
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") +
  # xlab("Week") +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  labs(title = "Nasal aspirate alpha diversity",
       subtitle = "Rarefied<br>No species removed")
ggsave("observed_richness_plot.png")

plot_richness(phylo_rar, 
              x = "Syndrome", #CHANGE ME, A is whatever factor you want on x-axis
              measures = "Shannon", #CHANGE ME, whatever alpha diversity measure you want. to have multiple, use: = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA)) + #example, get rid of x axis labels
  geom_violin(trim = TRUE, aes(fill = Syndrome)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(color = Syndrome)) + #optional. CHANGE ME, A is whatever factor to group box plots
  scale_color_hue(guide = "none") +
  facet_grid(.~Group, switch = "x", space = "free") + # facet_grid(.~Week, space = "free") +
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Shannon Bacterial Richness (SVs)") +
  # xlab("Week") +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_markdown()
  ) +
  labs(title = "Nasal aspirate alpha diversity",
       subtitle = "Rarefied<br>No species removed")

# two different richness plots grouped

plot1 <- plot_richness(phylo_rar, 
                       x = "Syndrome", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Observed", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim = TRUE, aes(fill = Syndrome)) + #optional. CHANGE ME, A is whatever factor to color violins
  # geom_boxplot(width = 0.1, aes(group = Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
        axis.title.x = element_blank()) +
  ylab("Observed Bacterial Richness (SVs)")
plot1

plot2 <- plot_richness(phylo_rar, 
                       x = "Syndrome", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Shannon", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  geom_violin(trim = TRUE, aes(fill = Syndrome)) + #optional. CHANGE ME, A is whatever factor to color violins
  # geom_boxplot(width = 0.1, aes(group = Treatment)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
        axis.title.x = element_blank()) +
  ylab("Shannon Bacterial Richness (SVs)")


# plot together
ggarrange(plot1,
          ggarrange(plot2, nrow = 1, labels = c("B")),
          ncol = 2, labels = "A", common.legend = TRUE)
ggsave("different_richness_plots_panels.png")

## Splitting by strepto y-n ------------

# phylo_rar:
# 64 samples, 1414 taxa

plot1 <- plot_richness(phylo_rar, 
              x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
              measures = "Observed", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Observed Bacterial Richness (SVs)") +
  ggtitle("All taxa")

plot2 <- plot_richness(phylo_clean_no_strep_rar, 
              x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
              measures = "Observed", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Observed Bacterial Richness (SVs)") +
  ggtitle("No Streptococcaceae")

plot3 <- plot_richness(phylo_clean_strep_rar, 
              x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
              measures = "Observed", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Observed Bacterial Richness (SVs)") +
  ggtitle("Only Streptococcaceae")


ggarrange(plotlist = list(plot1, plot2, plot3),
          labels = c("A", "B", "C"),
          nrow = 1,
          ncol = 3)
ggsave("observed_richness.png")

plot4 <- plot_richness(phylo_rar,
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Shannon", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) +
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Shannon Bacterial Richness (SVs)") +
  ggtitle("All taxa")

plot5 <- plot_richness(phylo_clean_no_strep_rar,
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Shannon", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) +
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Shannon Bacterial Richness (SVs)") +
  ggtitle("No Streptococcaceae")

plot6 <- plot_richness(phylo_clean_strep_rar,
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = "Shannon", #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) +
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Shannon Bacterial Richness (SVs)") +
  ggtitle("Only Streptococcaceae")


ggarrange(plotlist = list(plot4, plot5, plot6),
          labels = c("A", "B", "C"),
          nrow = 1,
          ncol = 3)
ggsave("shannon_richness.png")

plot7 <- plot_richness(phylo_rar, 
              x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
              measures = c("Observed", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Bacterial Richness") +
  ggtitle("All taxa")

plot8 <- plot_richness(phylo_clean_no_strep_rar,
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = c("Observed", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) +
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Bacterial Richness") +
  ggtitle("No Streptococcaceae")

plot9 <- plot_richness(phylo_clean_strep_rar,
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = c("Observed", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) +
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank()) +
  ylab("Bacterial Richness") +
  ggtitle("Only Streptococcaceae")


ggarrange(plotlist = list(plot7, plot8, plot9),
          labels = c("A", "B", "C"),
          nrow = 3,
          ncol = 1)

ggsave("shannon_and_observed.png")
