# This is the Qiime2 code I used for my analysis
# It is based almost entirely on the Moving Pictures tutorial on the Qiime website

# It was run in batches on the Ohio Supercomputer platform
# Submitted to the Pitzer cluster

# Note that many of these visuals have been re-plotted in ggPlot

#!/bin/bash
#SBATCH --job-name=fresh_run_part1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH -A PUOM0012

cd /users/PUOM0012/emilykibbler/project/fresh_run
module load miniconda3
conda activate qiime
module load sratoolkit/2.11.2

mkdir raw_files
cd raw_files

# Pull from SRA
for (( i = 268; i <= 343; i++ ))
  do
  fasterq-dump SRR16775$i
done

module load fastqc
fastqc *.fastq
source activate multiqc
multiqc *_fastqc.zip

cd ..

# Make the reads into Qiime-compatible format
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data ./paired-end-demux.qza \
  --o-visualization ./demux_seqs.qzv

# Trim
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trim-left-f 10 \
  --p-trim-left-r 10 \
  --p-trunc-len-f 284 \
  --p-trunc-len-r 224 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

# This rep-seqs table has the feature IDs and the corresponding sequences
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

#!/bin/bash
#SBATCH --job-name=alpha_div
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH -A PUOM0012

cd /users/PUOM0012/emilykibbler/project/fresh_run
module load miniconda3
conda activate qiime

# Rarefy to minimum sample depth
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 27000 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ./core-metrics-results/observed_features_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --o-visualization ./core-metrics-results/observed_statistics.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column SarsCov2 \
  --o-visualization core-metrics-results/unweighted-unifrac-covid-significance.qzv \
  --p-pairwise

qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file metadata.tsv \
  --p-custom-axes SarsCov2 \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-covid.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file metadata.tsv \
  --p-custom-axes SarsCov2 \
  --o-visualization core-metrics-results/bray-curtis-emperor-covid.qzv

qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 27000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization core-metrics-results/shannon-pd-group-significance.qzv


#!/bin/bash
#SBATCH --job-name=stats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH -A PUOM0012

cd /users/PUOM0012/emilykibbler/project/fresh_run
module load miniconda3
conda activate qiime

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file metadata.tsv \
  --o-filtered-table ancom-table.qza

qiime composition ancombc \
  --i-table ancom-table.qza \
  --m-metadata-file metadata.tsv \
  --p-formula 'SarsCov2' \
  --o-differentials ancombc-dx.qza

qiime composition da-barplot \
  --i-data ancombc-dx.qza \
  --p-significance-threshold 0.001 \
  --o-visualization da-barplot-dx.qzv


qiime longitudinal anova \
  --m-metadata-file ./core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --p-formula 'faith_pd ~ SarsCov2' \
  --o-visualization ./core-metrics-results/faiths_pd_anova.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column SarsCov2 \
  --o-visualization core-metrics-results/weighted-unifrac-covid-significance.qzv \
  --p-pairwise

qiime longitudinal anova \
  --m-metadata-file ./core-metrics-results/shannon_vector.qza \
  --m-metadata-file ./metadata.tsv \
  --p-formula 'shannon_entropy ~ SarsCov2' \
  --o-visualization ./core-metrics-results/shannon_anova.qzv

#!/bin/bash
#SBATCH --job-name=build_silva_classifier
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH -A PUOM0012

cd /users/PUOM0012/emilykibbler/project/fresh_run
module load miniconda3
conda activate qiime


qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza

qiime rescript cull-seqs \
    --i-sequences silva-138.2-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.2-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.2-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.2-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.2-ssu-nr99-seqs-discard.qza

qiime rescript dereplicate \
    --i-sequences silva-138.2-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.2-ssu-nr99-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138.2-ssu-nr99-tax-derep-uniq.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  silva-138.2-ssu-nr99-seqs-derep-uniq.qza \
  --i-reference-taxonomy silva-138.2-ssu-nr99-tax-derep-uniq.qza \
  --o-classifier silva-138.2-ssu-nr99-classifier.qza


#!/bin/bash
#SBATCH --job-name=part1_C
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH -A PUOM0012

cd /users/PUOM0012/emilykibbler/project/fresh_run
module load miniconda3
conda activate qiime


qiime feature-classifier classify-sklearn \
  --i-classifier ./silva-138.2-ssu-nr99-classifier.qza  \
  --i-reads rep-seqs.qza \
  --o-classification silva_taxonomy.qza

qiime metadata tabulate \
  --m-input-file silva_taxonomy.qza \
  --o-visualization silva_taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy silva_taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization silva_taxa-bar-plots.qzv


#!/bin/bash
#SBATCH --job-name=part1_C
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH -A PUOM0012

cd /users/PUOM0012/emilykibbler/project/fresh_run
module load miniconda3
conda activate qiime

# Green genes classifier did not need to be built like the Silva classifier
# Pull in the pre-formatted classifier from qiime website
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification gg_taxonomy.qza

qiime metadata tabulate \
  --m-input-file gg_taxonomy.qza \
  --o-visualization gg_taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy gg_taxonomy.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization gg_taxa-bar-plots.qzv

#!/bin/bash
#SBATCH --job-name=random_forest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=4:00:00
#SBATCH -A PUOM0012

cd /users/PUOM0012/emilykibbler/project/fresh_run
module load miniconda3
conda activate qiime

# Random forest modeling
qiime sample-classifier classify-samples \
  --i-table ./table.qza \
  --m-metadata-file ./metadata.tsv \
  --m-metadata-column SarsCov2 \
  --p-random-state 666 \
  --p-n-jobs 1 \
  --output-dir ./sample-classifier-results/

qiime sample-classifier heatmap \
  --i-table ./table.qza \
  --i-importance ./sample-classifier-results/feature_importance.qza \
  --m-sample-metadata-file ./metadata.tsv \
  --m-sample-metadata-column SarsCov2 \
  --p-group-samples \
  --p-feature-count 100 \
  --o-heatmap ./sample-classifier-results/heatmap_100-features.qzv \
  --o-filtered-table ./sample-classifier-results/filtered-table_100-features.qza

  qiime metadata tabulate \
  --m-input-file ./sample-classifier-results/feature_importance.qza  \
  --o-visualization ./sample-classifier-results/feature_importance.qzv
