

## QA reads ----------

### Qiime+dada

# old, unmerged reads
# qiime_read_stats <- read.table("/Users/emilykibbler/Documents/Classes/BMS_690/cov-project/dada2_filter_stats.tsv", sep = "\t", header = TRUE)

qiime_read_stats <- read.table("filtering_reads_denoising_stats.tsv", sep = "\t", header = TRUE)
head(qiime_read_stats)

qiime_plotData <- qiime_read_stats[,c("sample.id", "input", "filtered", "non.chimeric")]

qiime_plotData <- pivot_longer(qiime_plotData, cols = c("input", "filtered", "non.chimeric"), 
                               names_to = "Step", 
                               values_to = "Reads")
qiime_plotData$Step <- factor(qiime_plotData$Step, 
                        levels = c("input", "filtered", "non.chimeric"), 
                        labels = c("Unfiltered reads", "Filtered and trimmed reads", "Nonchimeric reads"))

qiime_plotData %>%
    ggplot(aes(x = Step, y = Reads)) +
      geom_boxplot() +
      ggtitle("Reads at each QA step: Qiime2 + DADA2")
ggsave("qiime_reads_qa_steps.png")


# This will help me determine what to rarify to
summary(subset(qiime_plotData, Step == "Nonchimeric reads")$Reads)
# Output:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 42623   58264   68542   70556   77730  129642 
# cool, can rarify to 42k without losing any samples!
# LOL that was with unmerged reads
# The real stats are min 27679, median 44k, max 78k

### dada only -------

# Bind columns from filtered output, # of seqs/sample from the no.chim seq table, and the treatment factor, all into a new variable
filtoutput250 <- readRDS("/Users/emilykibbler/Desktop/projects/bms690/trims_250/filtoutput.rds")
seqtab_nochim250 <- readRDS("/Users/emilykibbler/Desktop/projects/bms690/trims_250/seqtab_nochim.rds")

# track <- cbind(filtoutput, rowSums(seqtab_nochim), meta$SarsCov2) 
track <- cbind(filtoutput, rowSums(seqtab_nochim))
track <- as.data.frame(track)
colnames(track) <- c("reads.in","filtered", "nonchimeras") 
head(track)
saveRDS(track, 'track.rds') 

track250 <- cbind(filtoutput250, rowSums(seqtab_nochim250))
track250 <- as.data.frame(track250)
colnames(track250) <- c("reads.in","filtered", "nonchimeras") 
head(track250)
saveRDS(track250, 'track250.rds')


# 8. Plot all reads along the QC workflow
# make a prettier plot by taking the data
# plotData <- as.data.frame(track) %>% gather(type, totals, reads.in, filtered, nonchimeras)

#order the types from raw reads to cleanest
# plotData$type <- factor(plotData$type, levels = c("reads.in", "filtered", "nonchimeras"))


# ggplot(plotData, aes(x = Sample_type, y = as.numeric(totals))) +
ggplot(track, aes(x = Sample_type, y = as.numeric(totals))) +
  geom_boxplot(aes(color = Sample_type)) +
  geom_jitter(aes(color = Sample_type), width = 0.1) + 
  facet_grid(cols = vars(type)) +
  scale_color_hue(name = "Sample type") +
  ylab("Reads") + 
  # xlab("QA stage") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 10),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.background = element_rect(fill = "gray85"),
        plot.title = element_text(size = 16, face = "bold")) +
  # scale_y_continuous(trans = "log") +
  ggtitle("Reads by filtering step")
# 
# plot_dat <- plotData[,2:3]
# colnames(plot_dat) <- c("Step", "Reads")
# plot_dat$method <- "DADA2"
# plot_dat$Step <- str_replace_all(plot_dat$Step, "filtered", "Filtered and trimmed reads")
# plot_dat$Step <- str_replace_all(plot_dat$Step, "reads.in", "Unfiltered reads")
# plot_dat$Step <- str_replace_all(plot_dat$Step, "nonchimeras", "Nonchimeric reads")
# temp <- as.data.frame(qiime_plotData[,2:3])
# temp$method <- "Qiime2"
# plot_dat <- rbind(plot_dat, temp)
# unique(plot_dat$Step)
# plot_dat$Reads <- as.numeric(plot_dat$Reads)
track250$method <- "DADA, 250"
track$method <- "DADA, custom"
track250$sample.id <- row.names(track250)
track250$sample.id <- str_remove_all(track$sample.id, "_1.fastq")
colnames(track250) <- c("input", "filtered", "non.chimeric", "method", "sample.id")
track$sample.id <- row.names(track)
track$sample.id <- str_remove_all(track$sample.id, "_1.fastq")
colnames(track) <- c("input", "filtered", "non.chimeric", "method", "sample.id")
qiime_plotData$method <- "Qiime2"
# reorder columns to match
qiime_plotData <- qiime_plotData[,c("input", "filtered", "non.chimeric", "method", "sample.id")]

plot_dat <- rbind(track, track250, qiime_plotData)
plot_dat <- subset(plot_dat, select = -sample.id)


plot_dat <- pivot_longer(plot_dat, cols = c("input", "filtered", "non.chimeric"), 
                               names_to = "Step", 
                               values_to = "Reads")
plot_dat <- as.data.frame(plot_dat)
plot_dat$Step <- factor(plot_dat$Step, 
                        levels = c("input", "filtered", "non.chimeric"))

plot_dat %>% ggplot(aes(x = Step, y = Reads, fill = method)) +
  geom_boxplot() +
  facet_grid(cols = vars(Step), scales = "free") +
  ggtitle("Read depth and method comparisons") +
  theme_bw() +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        panel.background = element_rect(color = "black", fill = "gray90"))
ggsave("trims_read_depths.png")



## Alpha diversity ------------

### Qiime analysis ------------

# faith <- read.table("faith_group_signif.tsv", sep = "\t", header = TRUE)
# faith <- rename(faith, "Value" = faith_pd)
# faith$Metric <- "Faith's PD"

evenness <- read.table("even_group_signif.tsv", sep = "\t", header = TRUE)
evenness <- rename(evenness, "Value" = pielou_evenness)
evenness$Metric <- "Pielou evenness"

shannon <- read.table("shannon_group_signif.tsv", sep = "\t", header = TRUE)
shannon <- rename(shannon, "Value" = shannon_entropy)
shannon$Metric <- "Shannon"

obs <- read.table("observed_stats.tsv", sep = "\t", header = T)
obs <- rename(obs, "Value" = observed_features)
obs$Metric <- "Observed"

alpha <- rbind(obs, evenness, shannon)

q <- ggplot(alpha, aes(x = SarsCov2, y = Value, fill = SarsCov2)) +
    geom_boxplot() +
    facet_wrap(~Metric, scales = "free") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          plot.title = element_text(size = 12),
          panel.background = element_rect(color = "black", fill = "gray90")) +
  ggtitle("Qiime2")
q
ggsave("qiime_obs_even_shannon_alpha.png")

### R analysis ----------
phylo_rar <- readRDS("phylo_rar.rds") 
# This version of phylo_rar was made using species assignment first, taking out stuff like mitochondria, then rarifaction
df <- estimate_richness(phylo_rar, measures = c("Shannon", "Observed"))
df$sample.id <- row.names(df)
df <- df[,c("Shannon", "Observed", "sample.id")] %>%
  pivot_longer(cols = c("Shannon", "Observed"), names_to = "Metric", values_to = "Value") 
df <- as.data.frame(df)

even <- evenness(phylo_rar, index = "pielou")
colnames(even) <- ("Value")
even$sample.id <- row.names(even)
even$Metric <- "Pielou's evenness"
even <- even[, c("sample.id", "Metric", "Value")]
df <- rbind(df, even)
df <- merge(df, meta, by = "sample.id")
df1 <- df

d1 <- df1 %>% ggplot(aes(x =  SarsCov2, y = Value, fill = SarsCov2)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~Metric, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        plot.title = element_text(size = 12),
        panel.background = element_rect(color = "black", fill = "gray90")) +
  ggtitle("DADA2: Species filtering then rarifaction")
d1



# This one is rarified before any species assignments
phylo_no_spec_rar <- readRDS("phylo_no_spec_rar.rds")
df <- estimate_richness(phylo_no_spec_rar, measures = c("Shannon", "Observed"))

df$sample.id <- row.names(df)
df <- df[,c("Shannon", "Observed", "sample.id")] %>%
  pivot_longer(cols = c("Shannon", "Observed"), names_to = "Metric", values_to = "Value") 
df <- as.data.frame(df)

even <- evenness(phylo_no_spec_rar, index = "pielou")
colnames(even) <- ("Value")
even$sample.id <- row.names(even)
even$Metric <- "Pielou's evenness"
even <- even[, c("sample.id", "Metric", "Value")]
df <- rbind(df, even)
df <- merge(df, meta, by = "sample.id")

d2 <- df %>% ggplot(aes(x =  SarsCov2, y = Value, fill = SarsCov2)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    facet_wrap(~Metric, scales = "free") +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          plot.title = element_text(size = 12),
          panel.background = element_rect(color = "black", fill = "gray90")) +
    ggtitle("DADA2: Rarefied, no species filtering")
d2
# ggsave("obs_shann_even_R_analysis.png")

ggarrange(plotlist = list(q, d1, d2), 
          common.legend = TRUE,
          legend = "right",
          ncol = 1) %>%
    annotate_figure(top = text_grob("Alpha diversity metrics, by method", size = 14))
ggsave("alpha_div_methods_fill.png")


## DESeq ----------------
sigtab <- readRDS("sigtab.rds")

ggplot(sigtab, aes(y = Genus, x = log2FoldChange, color = Phylum)) + #play with aesthetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(aes(size = baseMean), position = position_jitter()) + #scale size by mean relative abundance
  scale_size_continuous(range = c(2, 8)) +
  theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, size = 14), 
            axis.text.y = element_text(size = 12),
           # axis.title.y = element_blank(),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10),
        # legend.direction = "vertical",
        legend.position = "none",
            panel.border = element_rect(color = "gray", fill = NA, linewidth = 1)) + 
  xlab("log2 Fold Change") + 
  labs(size = "Mean Abundance",
            title = "Fold Change of Read Abundance",
           subtitle = "Determined by DESeq")



ggsave("deseq_no_legend.png")

### Trying to make a plot with the size legend on the top and the color legend on the bottom

# this is just the color part
plot1 <- ggplot(sigtab, aes(y = Genus, x = log2FoldChange, color = Phylum)) + #play with aesthetics to make graph informative
  geom_point(position = position_jitter()) +
  theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10),
        # axis.title.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        # legend.direction = "vertical",
        legend.position = "none",
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1)) + 
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  xlab("log2FoldChange")
plot1

# 3.1 create legend1
legend1_plot <- ggplot(sigtab, aes(y = Genus, x = log2FoldChange, color = Phylum)) + #play with aesthetics to make graph informative
  geom_point(position = position_jitter()) +
  theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10),
        # axis.title.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.direction = "horizontal",
        # legend.position = "bottom",
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1)) + 
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  xlab("log2FoldChange")
legend1_plot


# just the size part

plot2 <- ggplot(sigtab, aes(y = Genus, x = log2FoldChange)) + #play with aesthetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(aes(size = baseMean), position = position_jitter()) + #scale size by mean relative abundance
  scale_size_continuous(range = c(2, 8)) +
  theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10),
        # axis.title.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        # legend.direction = "vertical",
        legend.position = "none",
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1))

legend2_plot <- ggplot(sigtab, aes(y = Genus, x = log2FoldChange)) + #play with aesthetics to make graph informative
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(aes(size = baseMean), position = position_jitter()) + #scale size by mean relative abundance
  scale_size_continuous(range = c(2, 8)) +
  theme(axis.text.x = element_text(hjust = 0, vjust = 0.5, size = 10), 
        axis.text.y = element_text(size = 10),
        # axis.title.y = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.direction = "horizontal",
        # legend.position = "top",
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1))  +
  labs(size = "Mean Abundance")
legend2_plot

# 3.3 extract "legends only" from ggplot object
legend1 <- get_legend(legend1_plot)
legend1
legend2 <- get_legend(legend2_plot)
legend2

# 4.1 setup legends grid
legend1_grid <- cowplot::plot_grid(legend1, align = "h", nrow = 1)

# 4.2 add second legend to grid, specifying its location
legends <- legend1_grid +
  ggplot2::annotation_custom(
    grob = legend2,
    xmin = 0.5, xmax = 0.5, ymin = 0.55, ymax = 0.25
  )
legends
# I couldn't get it to all work together so I gave up. Come back to this later
# So I just saved the legends as a separate png and made a frankengraph in ppt
ggsave("legends.png")
plot
ggsave("no_legend_deseq_plot.png")
# 5. plot "plots" + "legends" (with legends in between plots)
cowplot::plot_grid(plot, legends,
                   nrow = 2,
                   rel_widths = c(1, 0.2)
)

## Differential abundance by Qiime --------

da <- read.csv("qiime_da.csv")
head(da)
da <- subset(da, select = -X)
tax <- read.table("silva_tax_assignments.tsv", sep = "\t", header = TRUE)
head(tax)
tax <- rename(tax, "id" = Feature.ID)
da <- merge(da, tax, by = "id")








da %>%
  ggplot(aes(y = reorder(Genus, lfc), x = lfc, fill = Phylum)) +
  geom_bar(stat = "identity") +
  geom_linerange(aes(xmin = error.lower, xmax = error.upper)) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 12)) +
  xlab("log2FoldChange") +
  ylab("Genus") +
  ggtitle("Differentially abundant taxa, determined by Qiime")

ggsave("da_qiime.png")


list(
  DESeq = DEseq_sig_SVs,
  Qiime = da$Sequence
) %>%
  ggVennDiagram() +
  # scale_x_continuous(expand = expansion(mult = .2)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("SVs found to be significant") +
  coord_flip()
ggsave("deseq_v_qiime_venn.png")

## Rarefaction curve -------------

rar_curve <- read.csv("observed_features_rar_curve.csv")

rar_curve <- pivot_longer(rar_curve, cols = 2:(ncol(rar_curve)-1), names_to = "Depth")
rar_curve <- as.data.frame(rar_curve)
rar_curve$Depth <- str_split_i(rar_curve$Depth, "\\.", 2)
rar_curve$Depth <- str_remove_all(rar_curve$Depth, "_iter")
rar_curve$Depth <- as.numeric(rar_curve$Depth)
rar_curve <- rename(rar_curve, "Features" = value)
head(rar_curve)
rar_curve %>% ggplot(aes(y = Features, x = Depth, color = SarsCov2)) +
  geom_point(position = position_jitterdodge(), alpha = 0.6, size = 1) +
  geom_smooth(method = "gam") +
  scale_x_continuous(breaks = (seq(0, max(rar_curve$Depth), by = 5000))) +
  theme(axis.text.x = element_text(size = 12),
        legend.position = "top") +
  ggtitle("Rarefaction curve")
ggsave("rar_curve.png")    


anova_res <- read.table("unweighted_covid_sig_anova_res.tsv", sep = "\t", header = TRUE)
head(anova_res)
anova_res <- subset(anova_res, select = -X)

anova_res %>% 
  # ggplot(aes(x = Group1, y = Distance, fill = stage(Group1, after_scale = alpha(color, 0.5)), color = Group2)) +
  ggplot(aes(x = Group1, y = Distance, fill = Group2)) +
    geom_boxplot() +
    # geom_boxplot(aes(fill = Group2), alpha = 0.5) +
    # geom_boxplot(aes(color = Group1, group_by = Group1)) +
    facet_wrap(~Group1, 
               scales = "free", 
               labeller = as_labeller(c(neg = "Distance to neg", pos = "Distance to pos"))) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    ggtitle("PERMANOVA results")
ggsave("beta_permanova.png")
