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

for (( i = 268; i <= 343; i++ ))
  do
  fasterq-dump SRR16775$i
done

module load fastqc
fastqc *.fastq
source activate multiqc
multiqc *_fastqc.zip

cd ..

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./manifest.tsv \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data ./paired-end-demux.qza \
  --o-visualization ./demux_seqs.qzv

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
  

  

