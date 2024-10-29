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

## Richness with significance ---------------





# richness plot observed SVs with lines to fit the view screen
plot_richness(phylo_rar, 
              x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
              measures = c("Observed", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + 
  # geom_violin(trim = TRUE, aes(fill = Diet)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = "Group")) + #optional. CHANGE ME, A is whatever factor to group boxplots
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") + 
  ylim(0,1500) + #define the y axis min/max
  # geom_segment(aes(x = 1, y = 1200, xend = 2, yend = 1200)) +  
  geom_text(x = 1.5, y = 1250, label = "***") # add a drawn in line and significance tag, adjusting the x and y coordinates till it fits where you want in the window.  Add another for each line to add. As written, this will fit the view window you have, if you adjust that your segments will not adjust with it.

phylo_decontam_rar <- readRDS("phylo_decontam_rar.rds")
phylo_decontam_no_strep_rar <- readRDS("phylo_decontam_no_strep_rar.rds")
phylo_decontam_strep_rar <- readRDS("phylo_decontam_strep_rar.rds")


plot10 <- plot_richness(phylo_decontam_rar, 
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = c("Observed", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 14)) +
  ylab("Bacterial Richness") +
  ggtitle("All taxa")

plot11 <- plot_richness(phylo_decontam_no_strep_rar,
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = c("Observed", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) +
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 14)) +
  ylab("Bacterial Richness") +
  ggtitle("No Streptococcaceae")

plot12 <- plot_richness(phylo_decontam_strep_rar,
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = c("Observed", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) +
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 14)) +
  ylab("Bacterial Richness") +
  ggtitle("Only Streptococcaceae")


combo <- ggarrange(plotlist = list(plot10, plot11, plot12),
          labels = c("A", "B", "C"),
          nrow = 3,
          ncol = 1)

annotate_figure(combo, top = text_grob("Decontam method alpha diversity plots", size = 16))

ggsave("decontam_shannon_and_observed.png")

source("/Users/emilykibbler/Desktop/projects/R/AVS_554/lab5_functions.R")

phylo <- readRDS("phylo.rds")


phylo_atb <- subset_and_trim(phylo, "Group", "IPD_ATB")
phylo_no_atb <- subset_and_trim(phylo, "Group", "IPD")
phylo_NC <- subset_and_trim(phylo, "Group", "controlneg")

dirty_atb_SVs <- colnames(as.data.frame(phylo_atb@otu_table))
dirty_no_atb_SVs <- colnames(as.data.frame(phylo_no_atb@otu_table))
dirty_NC_SVs <- colnames(as.data.frame(phylo_NC@otu_table))

list(
  Antibiotics = dirty_atb_SVs,
  No_Antibiotics = dirty_no_atb_SVs,
  Neg_CT = dirty_NC_SVs
) %>%
  ggVennDiagram() +
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(plot.title = element_text(face = "bold", size = 16)) +
    ggtitle("SVs -- before decontam and rarification")
  ggsave("dirty_venn_diagram_SVS_atb_no_atb.png")

phylo_decontam_rar <- readRDS("phylo_decontam_rar.rds")

phylo_decontam_rar_atb <- subset_and_trim(phylo_decontam_rar, "Group", "Antibiotics")
phylo_decontam_rar_no_atb <- subset_and_trim(phylo_decontam_rar, "Group", "No antibiotics")
atb_SVs <- row.names(as.data.frame(phylo_decontam_rar_atb@tax_table))
no_atb_SVs <- row.names(as.data.frame(phylo_decontam_rar_no_atb@tax_table))

list(
  Antibiotics = atb_SVs,
  No_Antibiotics = no_atb_SVs
) %>%
  ggVennDiagram() +
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(plot.title = element_text(face = "bold", size = 16)) +
    ggtitle("SVs -- after decontam and rarification")
  ggsave("venn_diagram_SVS_atb_no_atb.png")

clean_phylo_rarified <- readRDS("clean_phylo_rarified.rds")

phylo_clean_rar_atb <- subset_and_trim(clean_phylo_rarified, "Group", "Antibiotics")
phylo_clean_rar_no_atb <- subset_and_trim(clean_phylo_rarified, "Group", "No antibiotics")
clean_atb_SVs <- row.names(as.data.frame(phylo_clean_rar_atb@tax_table))
clean_no_atb_SVs <- row.names(as.data.frame(phylo_clean_rar_no_atb@tax_table))

list(
  Antibiotics = clean_atb_SVs,
  No_Antibiotics = clean_no_atb_SVs
) %>%
  ggVennDiagram( set_size = 3) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(plot.title = element_text(face = "bold", size = 16)) +
    ggtitle("SVs -- after clean and rarification")
  ggsave("clean_venn_diagram_SVS_atb_no_atb.png")

list(
  ATB_clean = clean_atb_SVs,
  ATB_decon = atb_SVs,
  No_ATB_clean = clean_no_atb_SVs,
  No_ATB_decon = no_atb_SVs
) %>%
  ggVennDiagram() +
    scale_y_continuous(expand = expansion(mult = .2)) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
    ggplot2::scale_fill_gradient(low = "blue",high = "yellow") +
    ggtitle("SVs -- after decontam/clean and rarification")
ggsave("venn_diagram_clean_decon_SVs.png")

# heat maps
top_300 <- prune_taxa(names(sort(taxa_sums(phylo_decontam_rar),TRUE)[1:300]), phylo_decontam_rar)

p <- plot_heatmap(top_300)
#FIXME there is a theme issue that is causing me to not be able to remove y axis labels
  # plot_heatmap(top_300) +
  #   theme(axis.y.text = element_blank())


atb_300 <- subset_samples(top_300, Group == "Antibiotics") # CHANGE ME to the column name that holds the variables associated with being a negative control
atb_p <- plot_heatmap(atb_300)
no_atb_300 <- subset_samples(top_300, Group == "No antibiotics")
no_atb_p <- plot_heatmap(no_atb_300)
ggarrange(plotlist = list(p, atb_p, no_atb_p),
          labels = c("All", "Antibiotics", "No antibiotics"),
          vjust = -0.5,
          ncol = 3,
          nrow = 1,
          align = "hv",
          common.legend = TRUE)
ggsave("heatmaps_top300_by_group.png")


plot_heatmap(phylo_decontam_rar, fill = "Phylum") +   
  facet_grid(~Group, space = "free", scales = "free") + 
  theme(legend.position = "bottom", axis.text.y = element_blank()) +
  ggtitle("Heat maps, grouped by Phylum")
ggsave("heat_map_phylum.png")

plot_heatmap(phylo_decontam_rar, fill = "Class") +   
  facet_grid(~Group, space = "free", scales = "free") + 
  theme(legend.position = "bottom", axis.text.y = element_blank()) +
  ggtitle("Heat maps, grouped by Class")
ggsave("heat_map_class.png")

plot_heatmap(top_300, fill = "Class") +   
  facet_grid(~Group, space = "free", scales = "free") + 
  theme(legend.position = "bottom", axis.text.y = element_blank()) +
  ggtitle("Heat maps, top 300, grouped by Class")
ggsave("heat_map_300_class.png")

## Dirty and clean taxa analysis --------------
phylo_dirty_with_species <- readRDS("phylo_dirty_with_species.rds")

df_dirty <- as.data.frame(phylo_dirty_with_species@tax_table)
table(df_dirty$Kingdom)
table(df_dirty$Phylum)
table(df_dirty$Class)
table(df_dirty$Order)
table(df_dirty$Family)

Ishaq_clean <- readRDS("Ishaq_phylo_clean_with_species.rds")
df_clean <- as.data.frame(Ishaq_clean@tax_table)

decontam_with_species <- readRDS("phylo_decontam_with_species.rds")
df_decontam <- as.data.frame(decontam_with_species@tax_table)


# kingdoms <- list(
#   dirty = df_dirty$Kingdom,
#   Ishaq_clean = df_clean$Kingdom,
#   decontam = df_decontam$Kingdom
# )
# ggVennDiagram(kingdoms) +
#   ggtitle("kingdoms")

classifications <- colnames(df_dirty)

for (i in 1:length(classifications)) {
   temp <- list(
     Dirty = df_dirty[,i],
     Ishaq_clean = df_clean[,i],
     Decontam = df_decontam[,i])
   p <- ggVennDiagram(temp) +
     ggtitle(classifications[i])
   print(p)
   ggsave(paste0(classifications[i], "_venn_diagram.png"), plot = p)
 }

plot_list <- list()

# combo <- ggarrange(plotlist = list(plot10, plot11, plot12),
#                    labels = c("A", "B", "C"),
#                    nrow = 3,
#                    ncol = 1)

