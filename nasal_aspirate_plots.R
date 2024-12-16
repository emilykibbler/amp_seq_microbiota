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

## Patient metadata ------
patient_data_correlation_summary <- readRDS("numeric_patient_data_correlation_summary.rds")
patient_data_correlation_summary$Variable <- str_replace_all(patient_data_correlation_summary$Variable, "_", " ")
patient_data_correlation_summary$Variable <- str_to_title(patient_data_correlation_summary$Variable) 
head(patient_data_correlation_summary)

p1 <-  subset(patient_data_correlation_summary, Variable != "Group" & Variable != "SampleID" ) %>%
        ggplot(aes(y = Variable, x = as.numeric(p.value), size = 2)) +
        geom_point(aes(color = as.numeric(estimate.cor))) +
        scale_size(guide = "none") +
        scale_color_gradient2(low = "blue",high = "red", midpoint = 0) +
        theme(panel.background = element_rect(fill = "gray"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 9)) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        labs(color = "Correlation coefficient",
             x = "p-value",
             title = "Correlation of Numeric Patient Metrics",
             subtitle = "Treatment Group vs Control Group") +
        geom_vline(xintercept = 0.05, color = "red") +
        annotate("text", x = 0.03, y = 10, label = "p = 0.05", angle = 90)
p1
# ggsave("numeric_correlation_patient_metrics.png")

chisq_summary <- readRDS("chisq_summary.rds")
chisq_summary$variable <- str_replace_all(chisq_summary$variable, "_", " ")
chisq_summary$variable <- str_to_title(chisq_summary$variable)
head(chisq_summary)

p2 <-  subset(chisq_summary, variable != "Group" & variable != "SampleID" ) %>%
        ggplot(aes(y = variable, x = as.numeric(p.value), size = 2)) +
        geom_point(aes(color = as.numeric(`statistic.X-squared`))) +
        scale_size(guide = "none") +
        scale_color_gradient2(low = "blue",high = "red", midpoint = 0) +
        theme(panel.background = element_rect(fill = "gray"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 9)) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        labs(color = "X-squared",
             x = "p-value",
             title = "Correlation of Categorical Patient Metrics",
             subtitle = "Treatment Group vs Control Group") +
        geom_vline(xintercept = 0.05, color = "red") +
        annotate("text", x = 0.03, y = 25, label = "p = 0.05", angle = 90)
# ggsave("categorical_correlation_patient_metrics.png")

ggarrange(plotlist = list(p1, p2), 
          labels = c("A", "B"), 
          nrow = 2, 
          ncol = 1, 
          heights = c(1, 1.25),
          align = "v")
ggsave("correlation_patient_metrics.png")

## Plot QC steps ------------------
file_info <- readRDS("file_info.rds")
# Filtered
# F
temp <- list.files(paste0(getwd(), "/", file_info$data_subset[1], "_filtered"), full.names = TRUE)
temp <- temp[grep("fastq.gz", temp)] # take out non-fastq files
p1 <- plotQualityProfile(temp, aggregate = TRUE) +
  ggtitle("Forward Reads")
print(p1)

# R
temp <- list.files(paste0(getwd(), "/", file_info$data_subset[2], "_filtered"), full.names = TRUE)
temp <- temp[grep("fastq.gz", temp)] # take out non-fastq files
p2 <- plotQualityProfile(temp, aggregate = TRUE) +
  ggtitle("Reverse Reads")

print(p2)
ggarrange(plotlist = list(p1, p2), labels = c("A", "B")) %>%
  annotate_figure(top = text_grob("Filtered Data Quality Plots", size = 16))
ggsave("filtered_quality_plots.png")

# Raw
# F
temp <- file_info[1,]$file_names_full[[1]] # the files we want to analyze for this step
raw_p1 <- plotQualityProfile(temp, aggregate = TRUE) +
  ggtitle("Forward Reads")
print(raw_p1)

# R

temp <- file_info[2,]$file_names_full[[1]] # the files we want to analyze for this step
raw_p2 <- plotQualityProfile(temp, aggregate = TRUE) +
  ggtitle("Reverse Reads")
print(raw_p2)

ggarrange(plotlist = list(raw_p1, raw_p2), labels = c("A", "B")) %>%
  annotate_figure(top = text_grob("Raw Data Quality Plots", size = 16))
ggsave("raw_data_qual_plots.png")


# Look at read numbers in and out of QC steps
track <- readRDS('tracked.rds')

plotData <- as.data.frame(track) %>% gather(type, totals, reads.in, filtered, nonchimeras)

#order the types from raw reads to cleanest
plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"), labels = c("Unfiltered reads", "Filtered and trimmed reads", "Nonchimeric reads"))



# plot with Sample_type along the X axis
# ggplot(plotData,aes(x = Sample_type, y = as.numeric(totals))) + geom_jitter(aes(color = type)) + 
#   ylab("Sequences") + 
#   xlab("Sample_type") +
#   # theme(axis.text.x = element_text(angle = 0, size = 12, vjust = -0.5)) +
#   ggtitle("Reads by sample type")
# ggsave("reads_filt_chim_sample_type.png")


# or, plot with QA stage along the X axis
plotData$Sample_type <- factor(plotData$Sample_type, levels = c("experimental", "negative"), labels = c("Patient", "Negative control"))

ggplot(plotData, aes(x = type, y = as.numeric(totals))) +
  geom_jitter(aes(color = Sample_type)) + 
  geom_boxplot(aes(color = Sample_type)) +
  scale_color_hue(name = "Sample type") +
  ylab("Sequences") + 
  xlab("QA stage") +
  theme(axis.text.x = element_text(angle = 0, size = 12),
        axis.title.x = element_text(size = 14),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        legend.position = "top",
        plot.title = element_text(size = 16, face = "bold")) +
  ggtitle("Reads by filtering step")
ggsave("reads_sample_type_QC_status.png")



## Dirty alpha diversity plot --------

phylo <- readRDS("phylo.rds")
# See below. The plot_richness function has some quirks
# I found it easier to make a df with estimate_richness and do my own plot from scratch
plot_richness(phylo, x = "Group",
              measures = c("Observed","Chao1", "Shannon"),
              color = "Group") +
  geom_violin() +
  # geom_boxplot() +
  scale_color_hue(labels = c("Negative control", "No antibiotics", "Antibiotics")) +
  theme_bw() + # a ggplot theme to make the graph look nice
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top",
        plot.title = element_text(size = 16, face = "bold")) +
  # scale_x_discrete(labels = c("Negative control", "No antibiotics", "Antibiotics")) +
  ggtitle("Alpha Diversity; Before Data Cleaning")
ggsave("initial_alpha_diversity_plot.png")

plot_richness(phylo, x = "Group",
              measures = c("Observed","Chao1", "Shannon")) +
  geom_boxplot(aes(color = "Group")) +
  scale_color_hue(labels = c("Negative control", "No antibiotics", "Antibiotics")) +
  theme_bw() + # a ggplot theme to make the graph look nice
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 16, face = "bold")) +
  # scale_x_discrete(labels = c("Negative control", "No antibiotics", "Antibiotics")) +
  ggtitle("Alpha Diversity; Before Data Cleaning")

# ycl6 in the github forum for the phyloseq package recommends using estimate_richness instead of plot_richness
# for more control over the plot
# let's try that
dirty_rich_df <- estimate_richness(phylo, measures = c("Observed","Chao1", "Shannon"))
# row names get a little messed up, function must not tolerate a mix of numeric and character
row.names(dirty_rich_df) <- row.names(phylo@sam_data)
dirty_rich_df <- cbind(dirty_rich_df, phylo@sam_data)
saveRDS(dirty_rich_df, "dirty_rich_df.rds")

esk_plot_richness <- function(input, the_method) { # some error messages might be nice to check inputs
  p <- ggplot(input, aes(x = Group)) +
    geom_violin(aes(y = !!sym(the_method), fill = Group)) +
    scale_fill_discrete(labels = c("Negative control", "No antibiotics", "Antibiotics")) +
    geom_boxplot(aes(y = !!sym(the_method)), width = 0.1) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "top")
  return(p)
}
p1 <- esk_plot_richness(dirty_rich_df, "Observed")
p2 <- esk_plot_richness(dirty_rich_df, "Chao1")
p3 <- esk_plot_richness(dirty_rich_df, "Shannon")
ggarrange(plotlist = list(p1, p2, p3),
          # labels = c("A", "B", "C"),
          common.legend = TRUE,
          nrow = 3,
          ncol = 1) %>%
  annotate_figure(top = text_grob("Before cleaning alpha diversity plots", size = 16))

ggsave("dirty_diversity_panel_plot.png")

## Dirty ordination plot ----
phylo_ord <- ordinate(phylo, #calculate similarities
                      method = "PCoA", #ordination type
                      "jaccard", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(phylo, phylo_ord, type = "samples", color = "Group") +
  geom_point(size = 2.5) +
  # scale_color_hue(name = "Sample group", labels = c("Lab negative control", "No antibiotics", "Antibiotics")) +
  theme(#panel.border = element_rect(color = "black", size = 1),
        panel.background = element_rect(fill = "gray84", color = "black"),
        legend.position = "top") +
  ggtitle("Ordination Plot, Before Cleaning", subtitle = "Jaccard Distances")
ggsave("jaccard_ordination_before_plot.png")

phylo_ord <- ordinate(phylo, #calculate similarities
                      method = "PCoA", #ordination type
                      "bray", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)
plot_ordination(phylo, phylo_ord, type = "samples", color = "Group") +
  geom_point(size = 2.5) +
  # scale_color_hue(name = "Sample group", labels = c("Lab negative control", "No antibiotics", "Antibiotics")) +
  theme(#panel.border = element_rect(color = "black", size = 1),
    panel.background = element_rect(fill = "gray90", color = "black"),
    legend.position = "top") +
  ggtitle("Ordination Plot, Before Cleaning", subtitle = "Bray-curtis Distances")
ggsave("bray_ordination_before_plot.png") 

## Clean ordination plot -----
phylo_ord <- ordinate(phylo_decontam_rar, #calculate similarities
                      method = "PCoA", #ordination type
                      "bray", binary = TRUE) #similarity type. Jaccard is binary, Bray can be binary (unweighted) or not (weighted)

plot_ordination(phylo_decontam_rar, phylo_ord, type = "samples", color = "Group") +
  geom_point(size = 2.5) +
  # scale_color_hue(name = "Sample group", labels = c("Lab negative control", "No antibiotics", "Antibiotics")) +
  theme(#panel.border = element_rect(color = "black", size = 1),
    panel.background = element_rect(fill = "gray90", color = "black"),
    legend.position = "top") +
  ggtitle("Ordination Plot, After Cleaning and Rarefaction", subtitle = "Bray-curtis Distances")
ggsave("bray_ordination_after_clean_plot.png") 

## Richness considering syndrome -------------
phylo_decontam_rar <- readRDS("phylo_decontam_rar.rds")
plot_richness(phylo_decontam_rar, 
              x = "Syndrome", 
              measures = c("Observed"), #whatever alpha diversity measure you want
              title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) + #make it look pretty
  theme(axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title.x = element_blank()) + 
  geom_violin(trim = TRUE, aes(fill = Syndrome)) + 
  geom_boxplot(width = 0.1, aes(color = Syndrome)) + 
  scale_color_hue(guide = "none") +
  facet_grid(.~Group, switch = "x", space = "free", scales = "free") + # facet_grid(.~Week, space = "free") +
  # theme(legend.position = "none") + #use to get rid of your legend
  ylab("Observed Bacterial Richness (SVs)") +
  # xlab("Week") +
  theme(
    plot.title = element_markdown(),
    plot.subtitle = element_markdown(),
    legend.position = "top"
  ) +
  labs(title = "Nasal aspirate alpha diversity",
       subtitle = "Decontaminated<br>Rarefied")
ggsave("syndrome_plus_observed_richness_plot.png")

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


## Clean alpha diversity ---------

sample_data(phylo_decontam_rar)$Group <- factor(sample_data(phylo_decontam_rar)$Group, 
                                                levels = c("No antibiotics", "Antibiotics"))

plot_richness(phylo_decontam_rar, 
                       x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                       measures = c("Chao1", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                       title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 14),
        panel.border = element_rect(color = "gray", fill = NA, size = 1)) +
  ylab("Alpha Diversity Metric") +
  ggtitle("Diversity Metrics: Cleaned and Rarefied Data")
ggsave("alpha_div_clean.png")

# trying to replicate the paper's figure closer
phylo_decontam_rar_df <- estimate_richness(phylo_decontam_rar, measures = c("Chao1", "Shannon"))
head(phylo_decontam_rar_df)
phylo_decontam_rar_df <- subset(phylo_decontam_rar_df, select = -se.chao1)
phylo_decontam_rar_df$Sample <- row.names(phylo_decontam_rar@sam_data)
phylo_decontam_rar_df$Group <- phylo_decontam_rar@sam_data$Group
phylo_decontam_rar_df <- pivot_longer(phylo_decontam_rar_df, c("Chao1", "Shannon"), names_to = "Metric")
phylo_decontam_rar_df$Group <- str_replace_all(phylo_decontam_rar_df$Group, "No antibiotics", "No antibiotics (n=22)")
phylo_decontam_rar_df$Group <- str_replace_all(phylo_decontam_rar_df$Group, "Antibiotics", "Antibiotics (n=54)")
phylo_decontam_rar_df$Group <- factor(phylo_decontam_rar_df$Group,
                                      levels = c("No antibiotics (n=22)", "Antibiotics (n=54)"))

subset(phylo_decontam_rar_df, Group == "No antibiotics (n=22)") %>%
  subset(Metric == "Chao1") %>%
    nrow() # yes, 22

subset(phylo_decontam_rar_df, Group == "No antibiotics (n=22)") %>%
  subset(Metric == "Chao1") %>%
    view()

phylo_decontam_rar_df %>% ggplot(aes(x = Group, y = value)) +
  # theme_set(theme_minimal(base_size = 14)) +
  # geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(aes(fill = Group), outlier.shape = NA) + 
  scale_fill_manual(values = c("#efe68a", "#679acc")) +
  facet_wrap(~Metric, scales = "free") +
  geom_point(position = position_jitter(width = 0.1)) +
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(plot.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylab("Alpha Diversity Metric") +
  ggtitle("Diversity Metrics: Cleaned and Rarefied Data")
ggsave("paper_dup_plot.png")

### Clean and dirty richness together -----



clean_alpha_plot <- plot_richness(phylo_decontam_rar, 
                                  x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                                  measures = c("Chao1", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                                  title = NULL) + 
                    theme_set(theme_minimal(base_size = 14)) +
                    geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
                    geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
                    theme(legend.position = "none") + #use to get rid of your legend
                    # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
                    theme(axis.title.x = element_blank(),
                          plot.title = element_text(size = 14),
                          panel.border = element_rect(color = "gray", fill = NA, size = 1)) +
                    ylab("Alpha Diversity Measure")
clean_alpha_plot

sample_data(phylo)$Group <- factor(sample_data(phylo)$Group, 
                                                levels = c("IPD_ATB", "IPD", "controlneg"), 
                                                labels = c("Antibiotics", "No antibiotics", "Lab Negative Control")) 
# in the paper they have no-atb on the left

sample_data(phylo)$Group <- factor(sample_data(phylo)$Group, 
                                   levels = c("No antibiotics", "Antibiotics", "Lab Negative Control"))
# save the updated object
saveRDS(phylo, "phylo.rds")

dirty_alpha_plot <- plot_richness(phylo, 
                                  x = "Group", #CHANGE ME, A is whatever factor you want on x-axis
                                  measures = c("Chao1", "Shannon"), #CHANGE ME, whatever richness you want. = c("Observed","Shannon")
                                  title = NULL) + 
  theme_set(theme_minimal(base_size = 14)) +
  geom_violin(trim = TRUE, aes(fill = Group)) + #optional. CHANGE ME, A is whatever factor to color violins
  geom_boxplot(width = 0.1, aes(group = Group)) + #optional. CHANGE ME, A is whatever factor to group boxplots
  theme(legend.position = "none") + #use to get rid of your legend
  # theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 6),
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 14),
        panel.border = element_rect(color = "gray", fill = NA, size = 1)) +
  ylab("Alpha Diversity Measure")
dirty_alpha_plot

ggarrange(plotlist = list(dirty_alpha_plot, clean_alpha_plot),
          labels = c("Before Decontamination", "After Decontamination and Rarefaction"),
          hjust = -.1,
          vjust = 1,
          font.label = list(face = "bold", size = 12),
          nrow = 2) %>%
  annotate_figure(top = text_grob("Alpha Diversity Plots; Before and After Data Cleaning", size = 16, face = "bold"))
ggsave("alpha_diversity_before_after_clean.png")
  


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
decontam_atb_SVs <- row.names(as.data.frame(phylo_decontam_rar_atb@tax_table))
decontam_no_atb_SVs <- row.names(as.data.frame(phylo_decontam_rar_no_atb@tax_table))

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

SVs_at_this_point <- c(clean_atb_SVs, clean_no_atb_SVs, decontam_no_atb_SVs, decontam_atb_SVs)
SVs_at_this_point <- unique(SVs_at_this_point)
length(SVs_at_this_point) # 1625

clean_SVs <- c(clean_atb_SVs, clean_no_atb_SVs)
clean_SVs <- unique(clean_SVs)
length(phylo_clean_SVs) # 1414

decontam_SVs <- c(decontam_no_atb_SVs, decontam_atb_SVs)
decontam_SVs <- unique(decontam_SVs)
length(decontam_SVs)

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
  "Clean, +antibiotics" = clean_atb_SVs,
  "Decontam, +antibiotics" = atb_SVs,
  "Clean, -antibiotics" = clean_no_atb_SVs,
  "Decontam, -antibiotics" = no_atb_SVs
) %>%
  ggVennDiagram() +
    scale_y_continuous(expand = expansion(mult = .1)) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          legend.position = "bottom") +
    ggplot2::scale_fill_gradient(low = "blue",high = "yellow") +
    ggtitle("SVs -- After Decontam or Ishaq Clean and Rarification")
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

### Dirty and clean taxa analysis --------------
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
# this one I'm doing separately because I do want the legend
plot_bar(tax_glom(phylo_decontam_rar, "Phylum"), fill = "Phylum") +
  facet_grid(~Group, space = "free", scales = "free")
ggsave("bar_plot_abundance_by_phylum_glom.png")



classes_to_plot <- c("Class", "Order", "Family", "Genus")
for (i in 1:length(classes_to_plot)) {
  p <- plot_bar(phylo_decontam_rar, fill = classes_to_plot[i]) +
    facet_grid(~Group, space = "free", scales = "free") +
    theme(legend.position = "none") +
    ggtitle(paste0(classes_to_plot[i], ", no glom"))
  print(p)
  ggsave(paste0(classes_to_plot[i], "_abundance_bar_chart.png"), plot = p)
}
# "Agglomerate" means "collect or form into a mass or group"
# agglomerate SVs by taxa level to get rid of all the blacklines on the graph. can be done for any level of taxonomy
# Example:
# phylo_decontam_rar_glom = tax_glom(phylo_decontam_rar, "Genus")

for (i in 1:length(classes_to_plot)) {
  p <- plot_bar(tax_glom(phylo_decontam_rar, classes_to_plot[i]), fill = classes_to_plot[i]) +
    facet_grid(~Group, space = "free", scales = "free") +
    theme(legend.position = "none") +
    ggtitle(classes_to_plot[i])
  print(p)
  ggsave(paste0(classes_to_plot[i], "_glom_abundance_bar_chart.png"), plot = p)
}

## Core SVs ------

### Heat maps ------
phylo.coreW_25 <- readRDS("phylo.coreW_25.rds")

atb_phylo.coreW_25 <- subset_samples(phylo.coreW_25, Group == "Antibiotics")
no_atb_phylo.coreW_25 <- subset_samples(phylo.coreW_25, Group == "No antibiotics")

#aggregate the genera so we don't get a lot of lines separating all the SVs
plot.gen <- aggregate_taxa(phylo.coreW_25, "Genus")

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-4), log10(.2), length = 10), 3)

p1 <- plot_core(plot.gen, 
          plot.type = "heatmap", 
          prevalences = prevalences, 
          detections = detections, min.prevalence = 1/10000) + 
  xlab("Detection Threshold (Relative Abundance (%))") +
  ylab("Bacterial SVs") +
  theme_minimal() + scale_fill_viridis() +
  ggtitle("Core SVs, All Patients")
# There's a warning, ignore it

atb_plot.gen <- aggregate_taxa(atb_phylo.coreW_25, "Genus")

p2 <- plot_core(atb_plot.gen, 
          plot.type = "heatmap", 
          prevalences = prevalences, 
          detections = detections, min.prevalence = 1/10000) + 
  xlab("Detection Threshold (Relative Abundance (%))") +
  ylab("Bacterial SVs") +
  theme_minimal() + scale_fill_viridis() +
  ggtitle("Core SVs, Antibiotics")

no_atb_plot.gen <- aggregate_taxa(no_atb_phylo.coreW_25, "Genus")

p3 <- plot_core(no_atb_plot.gen, 
          plot.type = "heatmap", 
          prevalences = prevalences, 
          detections = detections, min.prevalence = 1/10000) + 
  xlab("Detection Threshold (Relative Abundance (%))") +
  ylab("Bacterial SVs") +
  theme_minimal() + scale_fill_viridis() +
  ggtitle("Core SVs, No Antibiotics")

ggarrange(plotlist = list(p1, p2, p3),
          labels = c("A", "B", "C"),
          common.legend = TRUE,
          nrow = 1,
          ncol = 3,
          legend = "bottom") %>%
  annotate_figure(top = text_grob("Core SVs; Frequency >1/10000 and Prevalence > 0.25", size = 16))
ggsave("core_SV_heatmaps.png")

ggarrange(plotlist = list(p1, p2, p3),
          # labels = c("A", "B", "C"),
          common.legend = TRUE,
          nrow = 3,
          ncol = 1,
          legend = "bottom") # %>%
#annotate_figure(top = text_grob("Before cleaning alpha diversity plots", size = 16))
ggsave("vertical_core_SV_heatmaps.png")

# Again with stricter criteria

phylo.coreW_35 <- readRDS("phylo.coreW_35.rds")
atb_phylo.coreW_35 <- readRDS("atb_phylo.coreW_35.rds")
no_atb_phylo.coreW_35 <- readRDS("no_atb_phylo.coreW_35.rds")

#aggregate the genera so we don't get a lot of lines separating all the SVs
plot.gen <- aggregate_taxa(phylo.coreW_35, "Genus")

prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-4), log10(.2), length = 10), 3)

p1 <- plot_core(plot.gen, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = detections, min.prevalence = 1/10000) + 
  xlab("Detection Threshold (Relative Abundance (%))") +
  ylab("Bacterial SVs") +
  theme_minimal() + scale_fill_viridis() +
  ggtitle("Core SVs, All Patients")

atb_plot.gen <- aggregate_taxa(atb_phylo.coreW_35, "Genus")

p2 <- plot_core(atb_plot.gen, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = detections, min.prevalence = 1/10000) + 
  xlab("Detection Threshold (Relative Abundance (%))") +
  ylab("Bacterial SVs") +
  theme_minimal() + scale_fill_viridis() +
  ggtitle("Core SVs, Antibiotics")

no_atb_plot.gen <- aggregate_taxa(no_atb_phylo.coreW_35, "Genus")

p3 <- plot_core(no_atb_plot.gen, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = detections, min.prevalence = 1/10000) + 
  xlab("Detection Threshold (Relative Abundance (%))") +
  ylab("Bacterial SVs") +
  theme_minimal() + scale_fill_viridis() +
  ggtitle("Core SVs, No Antibiotics")

ggarrange(plotlist = list(p1, p2, p3),
          # labels = c("A", "B", "C"),
          common.legend = TRUE,
          nrow = 1,
          ncol = 3,
          legend = "bottom")# %>%
# annotate_figure(top = text_grob("Before cleaning alpha diversity plots", size = 16))
ggsave("35prev_horizontal_core_SV_heatmaps.png")

ggarrange(plotlist = list(p1, p2, p3),
          # labels = c("A", "B", "C"),
          common.legend = TRUE,
          nrow = 3,
          ncol = 1,
          legend = "bottom")# %>%
# annotate_figure(top = text_grob("Before cleaning alpha diversity plots", size = 16))
ggsave("35prev_vertical_core_SV_heatmaps.png")

 ### Box plots -------
atb_phylo.coreW_35 <- readRDS("atb_phylo.coreW_35.rds")
no_atb_phylo.coreW_35 <- readRDS("no_atb_phylo.coreW_35.rds")

amalgamate_SV_data <- function(input_otu_table, group_name){
  df <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(df) <- c("SV", "Group", "Abundance")
  for (i in 1:ncol(input_otu_table)) {  
    temp <- data.frame(matrix(nrow = nrow(input_otu_table), ncol = 3))
    colnames(temp) <- c("SV", "Group", "Abundance")
    temp$SV <- colnames(input_otu_table)[i]
    temp$Group <- group_name
    temp$Abundance <- as.data.frame(input_otu_table)[,i]
    df <- rbind(df, temp)
  }
  return(df)
}

dat <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(dat) <- c("SV", "Group", "Abundance")
dat <- rbind(dat, amalgamate_SV_data(atb_phylo.coreW_35@otu_table, "Antibiotics"))

dat <- rbind(dat, amalgamate_SV_data(no_atb_phylo.coreW_35@otu_table, "No antibiotics"))
head(dat)
dim(dat)

temp <- as.data.frame(atb_phylo.coreW_35@tax_table)
temp$SV <- row.names(temp)
dat <- left_join(dat, temp, by = "SV")


SVs_sig_diff_on_t_test <- c("GGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAGCTTATGGTTAATACCCATAAGCCCTGACGTTACCCACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTATTTAAGTCAGATGTGAAAGCCCCGGGCTTAACCTGGGAACTGCATCTGATACTGGATAACTAGAGTAGGTGAGAGGGGAGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCTCCCTGGCATCATACTGACACTGAGGTGCGAAAGCGTGGGTAGCAAACAG",
                            "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG",
                            "GGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACAAATGTGTAAGTAACTATGCACGTCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAG",
                            "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGTTTTGTGACGGTACCTGGAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCACGTCGTCTGTGAAATTCCACAGCTTAACTGTGGGCGTGCAGGCGATACGGGCTGACTTGAGTACTGTAGGGGTAACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTTACTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCAAACAG",
                            "GGAATATTGCACAATGGGCGCAAGCCTGATGCAGCGACGCCGCGTGGGGGATGACGGCCTTCGGGTTGTAAACTCCTTTCGCCAGGGACGAAGCGTTTTGTGACGGTACCTGGAGAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTGTCCGGAATTACTGGGCGTAAAGAGCTCGTAGGTGGTTTGTCACGTCGTCTGTGAAATTCCACAGCTTAACTGTGGGCGTGCAGGCGATACGGGCTGACTTGAGTACTGTAGGGGTAACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTTACTGGGCAGTTACTGACGCTGAGGAGCGAAAGCATGGGTAGCAAACAG")

core_SVs_sig_diff_on_wilcox_test <- c("GGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTTTGGTTGTAAAGCACTTTAAGTGGGGAGGAAAAGCTTATGGTTAATACCCATAAGCCCTGACGTTACCCACAGAATAAGCACCGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTATTTAAGTCAGATGTGAAAGCCCCGGGCTTAACCTGGGAACTGCATCTGATACTGGATAACTAGAGTAGGTGAGAGGGGAGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGATGGCGAAGGCAGCTCCCTGGCATCATACTGACACTGAGGTGCGAAAGCGTGGGTAGCAAACAG",
                                      "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTAGGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG",
                                      "GGAATCTTCGGCAATGGACGGAAGTCTGACCGAGCAACGCCGCGTGAGTGAAGAAGGTTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGAGTGTGAGAGTGGAAAGTTCACACTGTGACGGTATCTTACCAGAAAGGGACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTAACTTGAGTGCAAGAGGGGAGAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGATATATGGAGGAACACCGGTGGCGAAAGCGGCTCTCTGGCTTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGAGCAAACAG",
                                      "GGAATCTTCCGCAATGGGCGAAAGCCTGACGGAGCAACGCCGCGTGAGTGATGAAGGTCTTCGGATCGTAAAACTCTGTTATTAGGGAAGAACAAATGTGTAAGTAACTATGCACGTCTTGACGGTACCTAATCAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGCGCGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGGAAAACTTGAGTGCAGAAGAGGAAAGTGGAATTCCATGTGTAGCGGTGAAATGCGCAGAGATATGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGTCTGTAACTGACGCTGATGTGCGAAAGCGTGGGGATCAAACAG")

dat$signif <- NA
for (i in 1:nrow(dat)) {
  if (dat$SV[i] %in% SVs_sig_diff_on_t_test) {
    dat$signif[i] <- 1
  }
}
dat$wilcox <- NA
for (i in 1:nrow(dat)) {
  if (dat$SV[i] %in% core_SVs_sig_diff_on_wilcox_test) {
    dat$wilcox[i] <- 1.5
  }
}
          
p1 <- subset(dat, !is.na(Genus)) %>% 
  subset(Phylum == "Firmicutes") %>%
    ggplot() +
      geom_boxplot(aes(x = SV, y = Abundance, color = Genus)) +
      geom_point(
        aes(x = SV, y = signif),
        shape = "*", 
        size = 8, 
        color = "red",
        show.legend = FALSE) +
      ylim(c(0, 1.05)) +
      scale_color_paletteer_d("LaCroixColoR::Lime") +
      facet_grid(cols = vars(Phylum), rows = vars(Group), scales = "free", space = "free") +
      theme(axis.text.x = element_blank(),
            legend.position = "top",
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            panel.background = element_rect(color = "gray"))

p1
# "beyonce::X51"
# "tvthemes::Alexandrite"

p2 <- subset(dat, !is.na(Genus)) %>% 
  subset(Phylum != "Firmicutes") %>%
    ggplot() +
      geom_boxplot(aes(x = SV, y = Abundance, color = Genus)) +
      geom_point(
        aes(x = SV, y = signif),
        shape = "*", 
        size = 8,
        color = "red",
        show.legend = FALSE) +
      ylim(c(0, 1.05)) +
      scale_color_paletteer_d("lisa::OskarSchlemmer") +
      facet_grid(cols = vars(Phylum), rows = vars(Group), scales = "free", space = "free") +
      theme(axis.text.x = element_blank(),
            legend.position = "top",
            panel.border = element_rect(color = "black", fill = NA, size = 1),
            panel.background = element_rect(color = "gray"))

p2
ggarrange(plotlist = list(p1, p2),
          # labels = c("A", "B"),
          nrow = 2,
          ncol = 1) %>%
  annotate_figure(top = text_grob("Core SVs: >1/10,000 Frequency and >0.35 Prevalence", 
                                  size = 19))
ggsave("top_13_core_SVs_boxplot.png")

# playing with where the stars are
for (i in 1:nrow(dat)) {
  if (!is.na(dat$signif[i])) {
    dat$signif[i] <- 1.2
  }
}

for (i in 1:nrow(dat)) {
  if (!is.na(dat$wilcox[i])) {
    dat$wilcox[i] <- 3
  }
}

dat %>% ggplot() +
  geom_boxplot(aes(x = SV, # reorder(SV, Abundance, decreasing = TRUE),
                   y = Abundance, 
                   # fill = Genus,
                   color = Genus),
               # labels = Genus,
               outlier.shape = NA) +
  geom_point(aes(x = SV,
                 y = Abundance, 
                 color = Genus),
             size = 1,
             position = position_jitter(width = 0.2),
             show.legend = FALSE) +
  geom_point(
    aes(x = SV, y = signif),
    shape = "*", 
    size = 8,
    color = "red",
    show.legend = FALSE) +
  scale_y_continuous(trans = "log10", "Abundance, log10 scale", sec.axis = sec_axis(~ . , name = "Treatement Group")) +
  xlab("SV; Grouped by Phylum") +
  # scale_color_paletteer_d("lisa::OskarSchlemmer") +
  facet_grid(rows = vars(Group), cols = vars(Phylum), space = "free", scales = "free") +
  theme_bw() + # this puts the facet names in nice boxes
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "gray84", color = "black")
        # panel.border = element_rect(color = "black", fill = NA, size = 1)
  )  +
  ggtitle("Core SVs: >1/10,000 Frequency and >0.35 Prevalence", subtitle = "Defined Using core Function of Microbiome Package")
ggsave("core_SVs_log_y_scale.png")

dat %>% ggplot() +
  geom_boxplot(aes(x = SV, # reorder(SV, Abundance, decreasing = TRUE),
                   y = Abundance, 
                   # fill = Genus,
                   color = Genus),
               # labels = Genus,
               outlier.shape = NA) +
  geom_point(aes(x = SV,
                 y = Abundance, 
                 color = Genus),
             size = 1,
             position = position_jitter(width = 0.2),
             show.legend = FALSE) +
  # geom_point(
  #   aes(x = SV, y = signif),
  #   shape = "*", 
  #   size = 8,
  #   color = "red",
  #   show.legend = FALSE) +
  geom_point(
    aes(x = SV, y = wilcox),
    shape = "*", 
    size = 6,
    color = "red",
    show.legend = FALSE) +
  scale_y_continuous(trans = "log10", "Abundance, log10 scale", sec.axis = sec_axis(~ . , name = "Treatement Group")) +
  xlab("SV; Grouped by Phylum") +
  # scale_color_paletteer_d("lisa::OskarSchlemmer") +
  facet_grid(rows = vars(Group), cols = vars(Phylum), space = "free", scales = "free") +
  theme_bw() + # this puts the facet names in nice boxes
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "gray84", color = "black")
        # panel.border = element_rect(color = "black", fill = NA, size = 1)
  )  +
  ggtitle("Core SVs: >1/10,000 Frequency and >0.35 Prevalence")
ggsave("core_SVs_log_y_wilcox.png")

### Simple t-tests -----------

sig_SVs <- readRDS("t_test_sig_SVs.rds")

sig_SVs_adj <- sig_SVs
for (i in 1:nrow(sig_SVs_adj)) {
  if (sig_SVs_adj$Abundance[i] == 0) {
    sig_SVs_adj$Abundance[i] <- 0.0001
  }
}


  # sig_SVs %>% subset(mean_abundance > 0.001) %>% ggplot() +
  #   geom_boxplot(aes(x = SV, y = Abundance, color = Genus), outlier.shape = NA) +
  #   geom_point(aes(x = SV, y = Abundance, color = Genus), alpha = 0.8, size = 1, position = position_jitter(width = 0.1)) +
  #   # facet_grid(cols = vars(Phylum), rows = vars(Group), scales = "free", space = "free") +
  #   facet_grid(rows = vars(Group)) +
  #   theme(axis.text.x = element_blank(),
  #         legend.position = "bottom",
  #         panel.border = element_rect(color = "black", fill = NA, size = 1),
  #         panel.background = element_rect(color = "gray")) +
  #   ggtitle("SVs with p<0.05 Abundance Differences and Mean Abundance > 0.001")
  # ggsave("t-test_top_SVs.png")

sig_SVs_adj$special <- NA
for (i in 1:nrow(sig_SVs_adj)) {
  if (sig_SVs_adj$Genus[i] == "Streptococcus") {
    sig_SVs_adj$Genus[i] <- "Streptococcus +"
    sig_SVs_adj$special[i] <- 1}
}

for (i in 1:nrow(sig_SVs_adj)) {
  if (!is.na(sig_SVs_adj$special[i])) {
    sig_SVs_adj$special[i] <- 1.1}
}
# Solo plot:

sig_SVs_adj %>% ggplot() +
  geom_boxplot(aes(x = SV, # reorder(SV, Abundance, decreasing = TRUE),
                   y = Abundance, 
                   # fill = Genus,
                   color = Genus),
               # labels = Genus,
               outlier.shape = NA) +
  geom_point(aes(x = SV,
                   y = Abundance, 
                   color = Genus),
              size = 1,
              position = position_jitter(width = 0.2),
             show.legend = FALSE) +
  geom_point(
    aes(x = SV, 
        y = special),
    shape = "+", 
    size = 6,
    color = "red",
    show.legend = FALSE) +
  scale_y_continuous(trans = "log10", "Abundance, log10 scale", sec.axis = sec_axis(~ . , name = "Treatement Group")) +
  xlab("SV; Grouped by Phylum") +
  # scale_color_paletteer_d("lisa::OskarSchlemmer") +
  facet_grid(rows = vars(Group), cols = vars(Phylum), space = "free", scales = "free") +
  theme_bw() + # this puts the facet names in nice boxes
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "gray84", color = "black")
        # panel.border = element_rect(color = "black", fill = NA, size = 1)
        )  +
  ggtitle("SVs with Significant Difference by t-test")
ggsave("facets_t-test_signif_SVs.png")

# to be combined with the wilcox-significant plot
t_sig_svs_plot <- sig_SVs_adj %>% ggplot() +
  geom_boxplot(aes(x = SV, # reorder(SV, Abundance, decreasing = TRUE),
                   y = Abundance, 
                   # fill = Genus,
                   color = Genus),
               # labels = Genus,
               outlier.shape = NA) +
  geom_point(aes(x = SV,
                 y = Abundance, 
                 color = Genus),
             size = 1,
             position = position_jitter(width = 0.2),
             show.legend = FALSE) +
  geom_point(
    aes(x = SV, 
        y = special),
    shape = "+", 
    size = 6,
    color = "red",
    show.legend = FALSE) +
  scale_y_continuous(trans = "log10", "Abundance, log10 scale") +
  xlab("SV; Grouped by Phylum") +
  # scale_color_paletteer_d("lisa::OskarSchlemmer") +
  facet_grid(rows = vars(Group), cols = vars(Phylum), space = "free", scales = "free") +
  theme_bw() + # this puts the facet names in nice boxes
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "gray84", color = "black")
        # panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

### Wilcoxon tests --------

wilcox_sig_SVs <- readRDS("wilcox_sig_SVs.rds")

sig_SVs_adj <- wilcox_sig_SVs
for (i in 1:nrow(sig_SVs_adj)) {
  if (sig_SVs_adj$Abundance[i] == 0) {
    sig_SVs_adj$Abundance[i] <- 0.0001
  }
}


sig_SVs_adj$special <- NA
for (i in 1:nrow(sig_SVs_adj)) {
  if (sig_SVs_adj$Genus[i] == "Streptococcus") {
    sig_SVs_adj$Genus[i] <- "Streptococcus +"
    sig_SVs_adj$special[i] <- 1.1}
}

for (i in 1:nrow(sig_SVs_adj)) {
  if (!is.na(sig_SVs_adj$special[i])) {
    sig_SVs_adj$special[i] <- 2}
}

sig_SVs_adj %>% ggplot() +
  geom_boxplot(aes(x = SV, # reorder(SV, Abundance, decreasing = TRUE),
                   y = Abundance, 
                   # fill = Genus,
                   color = Genus),
               # labels = Genus,
               outlier.shape = NA) +
  geom_point(aes(x = SV,
                 y = Abundance, 
                 color = Genus),
             size = 1,
             position = position_jitter(width = 0.2),
             show.legend = FALSE) +
  geom_point(
    aes(x = SV, 
        y = special),
    shape = "+", 
    size = 6,
    color = "red",
    show.legend = FALSE) +
  scale_y_continuous(trans = "log10", "Abundance, log10 scale", sec.axis = sec_axis(~ . , name = "Treatement Group")) +
  xlab("SV; Grouped by Phylum") +
  # scale_color_paletteer_d("lisa::OskarSchlemmer") +
  facet_grid(rows = vars(Group), cols = vars(Phylum), space = "free", scales = "free") +
  theme_bw() + # this puts the facet names in nice boxes
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "gray84", color = "black")
        # panel.border = element_rect(color = "black", fill = NA, size = 1)
  )  +
  ggtitle("SVs with Significant Difference by Wilcoxon Test")
ggsave("facets_wilcox_signif_SVs.png")

wilcox_sig_svs_plot <- sig_SVs_adj %>% ggplot() +
  geom_boxplot(aes(x = SV, # reorder(SV, Abundance, decreasing = TRUE),
                   y = Abundance, 
                   # fill = Genus,
                   color = Genus),
               # labels = Genus,
               outlier.shape = NA) +
  geom_point(aes(x = SV,
                 y = Abundance, 
                 color = Genus),
             size = 1,
             position = position_jitter(width = 0.2),
             show.legend = FALSE) +
  geom_point(
    aes(x = SV, 
        y = special),
    shape = "+", 
    size = 6,
    color = "red",
    show.legend = FALSE) +
  scale_y_continuous(trans = "log10", "Abundance, log10 scale") +
  xlab("SV; Grouped by Phylum") +
  # scale_color_paletteer_d("lisa::OskarSchlemmer") +
  facet_grid(rows = vars(Group), cols = vars(Phylum), space = "free", scales = "free") +
  theme_bw() + # this puts the facet names in nice boxes
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        panel.background = element_rect(fill = "gray84", color = "black")
        # panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) 

ggarrange(plotlist = list(t_sig_svs_plot, wilcox_sig_svs_plot),
          labels = c("t-test", "Wilcox test"),
          ncol = 1)
# re-think; this makes it confusing whether they are the same SVs or not
# and the colors are different for some of the same genuses


### Venn diagrams -----

atb_phylo.coreW
no_atb_phylo.coreW

core_atb_SVs <- row.names(as.data.frame(atb_phylo.coreW@tax_table))
core_atb_genus <- as.data.frame(atb_phylo.coreW@tax_table)$Genus
core_no_atb_SVs <- row.names(as.data.frame(no_atb_phylo.coreW@tax_table))
core_no_atb_genus <- as.data.frame(no_atb_phylo.coreW@tax_table)$Genus

list(
  Antibiotics = core_atb_genus,
  No_Antibiotics = core_no_atb_genus
) %>%
  ggVennDiagram( set_size = 3) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(plot.title = element_text(face = "bold", size = 16)) +
  ggtitle("Core SVs by Genus")

### Strepto focus --------

strep_abun_df <- readRDS("strep_abun_df.rds")
length(unique(strep_abun_df$SV)) # 61 is going to make a crazy plot
summary(strep_abun_df$mean_abundance) # 3rd quartile starts at 0.0026
for (i in 1:nrow(strep_abun_df)) { # makes log transform possible
  if (strep_abun_df$Abundance[i] == 0) {
    strep_abun_df$Abundance[i] <- 0.001
  }
}

for (i in 1:nrow(strep_abun_df)) { # makes log transform possible
  if (is.na(strep_abun_df$Species[i])) {
    strep_abun_df$Species[i] <- "Undefined"
  }
}

strep_abun_df %>%
  subset(mean_abundance > 0.006) %>% 
    ggplot(aes(x = SV_name, y = Abundance)) +
    geom_boxplot(aes(color = Species), outlier.shape = NA) +
    geom_point(aes(color = Species), position = position_jitter(width = 0.1)) +
    scale_y_continuous(transform = "log10", limits = c(0.0001,1.1)) +
    facet_wrap(~Group, ncol = 1, scales = "free") +
    xlab("SV") +
    ylab("Relative Abundance of Streptococcus Reads") +
    theme_bw() +
    theme(axis.text.x = element_blank()) +
    ggtitle("Most Abundant Streptococcus SVs")


## DEseq plot -----
sigtab <- readRDS("sigtab.rds")
ggplot(sigtab, aes(y = Genus, x = log2FoldChange, color = Phylum)) + #play with aesthetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(aes(size = baseMean), position = position_jitter()) + #scale size by mean relative abundance
  scale_size_continuous(range = c(3, 8)) +
  theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 11),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.border = element_rect(color = "gray", fill = NA, size = 1)) + 
  xlab("log2 Fold Change") + 
  labs(size = "Mean Sequence Abundance",
       title = "Fold Change of Read Abundance",
       subtitle = "Antibiotic-treated Compared to Control Group")
  # theme_minimal()
ggsave("deseq_fc_analysis.png")

## RF -------
rf_model <- readRDS("rf_model.rds")

# set color palette  
col.bro <- (rainbow(6))
# add white to that color list
col.bro <- append(col.bro, "#ffffff", after = 6)


# adjusted plot for factorial data, recommend using sample ID as 'factor A' (x value)
ggplot(rf_model, aes(as.factor(Sample), reorder(variable, rescale, mean))) + 
  theme_minimal() + 
  facet_grid(.~Group, space = "free", scales = "free") + #set up graph facets to separate out levels of FactorA
  geom_tile(aes(fill = rescale), color = "gray") + #add the heatmap coloring 
  scale_fill_gradientn(colors = rev(col.bro), na.value = 'white') + #use the preset color palette
  labs(fill = "Log abundance") + #rename legend heading
  theme(legend.position = 'bottom',
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(size = 2),
        axis.text.y = element_text(size = 6),
        panel.background = element_rect(fill = "gray84", color = "black")) +
  ylab('Predictor Taxa') +
  xlab('Sample') +
  ggtitle("Random Forest")
ggsave("rf_heatmap.png")

## Beta ordination ------
### PCOA -----
uJ_pcoa <- readRDS("uJ_pcoa.rds")

plot_ordination(phylo_decontam_rar, uJ_pcoa, type = "samples") + #CHANGE ME
  geom_point(aes(color = as.factor(Group))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  # theme_minimal() +
  theme(panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        legend.position = "top") +
  # scale_color_viridis(discrete = TRUE, option = "A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group = Group, color = Group)) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Treatment Group",
       title = "PCOA Ordination",
       subtitle = "After Cleaning and Rarefaction") 
ggsave("pcoa_ordination_after_clean_and_rar_with_circles.png")

### NMDS ------
uJ_nmds <- readRDS("uJ_nmds.rds")
plot_ordination(phylo_decontam_rar,
                uJ_nmds, type = "samples",
                color = "Group") + #CHANGE ME
  geom_point(aes(color = Group)) +
  theme(legend.position = "top",
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1)) +
  # geom_point(aes(size = as.numeric(Weight_kg), color=as.factor(Week), alpha = as.factor(Week))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  # scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group = Group, color = Group)) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Treatment Group",
       title = "NMDS Ordination",
       subtitle = "After Cleaning and Rarefaction") #rename the headers in the legend
  # ggtitle("NMDS Ordination") #+
# labs(color = "Weeks 0 to 2", size = "Weight in kg", shape = "Diet", linetype="Diet", alpha= "Week") #rename the headers in the legend
ggsave("nmds_ordination_after_clean_and_rar_with_circles.png")

# Unconstrained:
# Does not enforce spacing
# "Letting the data fall where it falls"
# Constrained:
# Shows how important different factors are
# Forces a scaling
# Might want to put in both

uJ_nmds_bray <- ordinate(phylo_decontam_rar, #calculate similarities
                    method = "NMDS", #ordination type
                    "bray", binary = TRUE)

plot_ordination(phylo_decontam_rar,
                uJ_nmds_bray, type = "samples",
                color = "Group") + #CHANGE ME
  geom_point(aes(color = Group)) +
  theme(legend.position = "top",
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1)) +
  # geom_point(aes(size = as.numeric(Weight_kg), color=as.factor(Week), alpha = as.factor(Week))) + # resize the points based on a numeric factor, and make light/dark based on another factor
  # scale_color_viridis(discrete = TRUE, option="A", begin = 0, end = 0.75) + #begin/end indicates where on the color scale to use, A refers to 'magma' color palette
  stat_ellipse(aes(group = Group, color = Group)) + #add circles around a particular grouping, and make the circle lines different
  labs(color = "Treatment Group",
       title = "NMDS Ordination",
       subtitle = "After Cleaning and Rarefaction; Bray-Curtis Distance") #rename the headers in the legend
# ggtitle("NMDS Ordination") #+
# labs(color = "Weeks 0 to 2", size = "Weight in kg", shape = "Diet", linetype="Diet", alpha= "Week") #rename the headers in the legend
ggsave("bray_nmds_ordination_after_clean_and_rar_with_circles.png")



## Order test ---------

phylo_species_then_decontam_rar <- readRDS("phylo_species_then_decontam_rar.rds")
phylo_decontam_rar <- readRDS("phylo_decontam_rar.rds")


list(
  "Decontamination then Species Assignment" = row.names(phylo_decontam_rar@tax_table),
  "Species Assignment then Decontamination" = row.names(phylo_species_then_decontam_rar@tax_table)
) %>%
  ggVennDiagram( set_size = 3) +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.1, vjust = -12),
        legend.position = "none") +
  coord_flip() +
  ggtitle("SVs Found")

# add ggvenn library to packages script if I end up using any of this

# list(
#   row.names(phylo_decontam_rar@tax_table),
#   row.names(phylo_species_then_decontam_rar@tax_table)
# ) %>%
#   ggvenn()
# 
# # 
# # ggvenn(columns = list(row.names(phylo_decontam_rar@tax_table),
# #                         row.names(phylo_species_then_decontam_rar@tax_table)))
#   
# list(
#   row.names(phylo_decontam_rar@tax_table),
#   row.names(phylo_species_then_decontam_rar@tax_table)
# ) %>%
#   ggvenn(columns = c("Decontamination then Species Assignment",
#                      "Species Assignment then Decontamination"))

install.packages("eulerr")
library(eulerr)

dat <- list("Decontamination then Species Assignment" = row.names(phylo_decontam_rar@tax_table),
  "Species Assignment then Decontamination" = row.names(phylo_species_then_decontam_rar@tax_table))

# plot(euler(dat),
#      fills = list(fill = c("Decontamination then Species Assignment" = "#54928D",
#                   "Species Assignment then Decontamination" = "#941C50"),
#                   alpha = 0.9))

# I don't like that either, the output is not a ggplot object


install.packages("VennDiagram")
library(VennDiagram)

# this kinda sucks because it goes to file and not the plot window
venn.diagram(x = dat, 
             filename = "order_of_operations.png",
             imagetype = "png",
             fill = c("yellow", "purple"),
             main = "SVs Found, by Order of Cleaning Steps",
             main.cex = 1.5,
             cat.pos = c(320, 140),
             cat.dist = c(0.1,0.1),
             cat.just = list(c(0.2,-0.1), c(0.8,0)))
             




