# NaMeco
This pipeline was written for processing full-length 16S Nanopore reads

## Citation
If you used our pipelines, please cite this paper: XXX

## Why this pipeline was created?
16S long Nanopore reads have an advantage in terms of length compared with Illumina short reads and have the potential for better taxonomical annotations. However, in practice, due to the relatively high error rate, reads produced by Nanopore are more challenging to handle. The combination of longer length and higher error rate results in a higher number of unique reads. When we tried to use blast on unclustered Nanopore reads, a dataset with 170 samples was not finished on the HPC with 24 threads in a month. However, clustering Nanopore reads at a 97% similarity threshold did not improve the situation, only slightly decreasing the number of unique features (sequences). 

Here, we decided to merge ideas from such pipelines as NanoCLUST, NGSpeciesID and Natrix2 and create one pipeline ~to rule them all~ that will combine the advantages of above mentioned tools. So, it will preprocess the reads, count kmers and then perform clustering with UMAP + HDBscan sample by sample. After that, form each cluster of each sample 50 representatives will be randomly selected for additional clustering between samples to clsuter clusters. New clusters, this time already "shared" between samples, will be polished with combination of SPOA, Racon and Medaka. Taxonomy will be assigned based on either GTDB or NCBI databases.

## Dependencies 
Linux OS with conda installed (anaconda3 or miniconda3). 
- pandas=2.2.2
- scikit-learn=1.4.2
- python=3.10.8
- chopper=0.7.0
- racon=1.5.0
- medaka=1.11.3
- minimap2=2.28
- umap-learn=0.5.5
- biopython=1.83
- matplotlib=3.8.4
- blast=2.15
- spoa=4.1.4

## Installation
During installation, new conda environment NaMeco will be created, which includes all dependencies.

This pipeline can be installed with the following script:


```python
wget ... #edit
conda env create --file ... #edit
```

## Running the pipeline
This pipeline takes Nanopore reads as input in fastq format. It will automatically recognize .fastq, .fq extensions. Reads also can be gziped (.fastq.gz, .fq.gz)

To run the pipeline, please provide the path to raw reads and adjust threads. The rest of the parameters may be changed if needed. 



```python
nameco --help

required arguments:
  --inp_dir INP_DIR     Path to the folder with reads, absolute or relative.
                        Reads should be in the fastq or fq format, gziped or
                        not

optional arguments:
  --out_dir OUT_DIR     Path to the directory to store output files, absolute
                        or relative. If not provided, folder "Nameco_out" will
                        be created in working directory
  --threads THREADS     The number of threads/cpus (default 2)
  --qc                  Run chopper for quality control (default)
  --no-qc               Skip chopper for quality control
  --phred PHRED         Minimum phred average score for chopper (default 8)
  --min_length MIN_LENGTH
                        Minimum read length for chopper (default 1300)
  --max_length MAX_LENGTH
                        Maximum read length for chopper (default 1700)
  --kmer KMER           K-mer length for clustering (default 5)
  --no-low              Don't restrict RAM for UMAP (default)
  --low                 Reduce RAM usage by UMAP
  --cluster_size CLUSTER_SIZE
                        Minimum cluster size for HDBscan (default 50)
  --subsample SUBSAMPLE
                        Subsample bigger than that threshold clusters for
                        consensus creation and polishing by Racon and Medaka
                        (default 1000)
  --random_state RANDOM_STATE
                        Random state for subsampling (default 42)
  --database {GTDB,NCBI}
                        Database for taxonomy assignment (default GTDB). Only
                        GTDB or NCBI are currently supported
  --db_path DB_PATH     Path to store/existing database (default
                        $out_dir/$database). Please use only databases,
                        created by previous NaMeco run to avoid errors
```

If the run was killed, it can be relaunched with the same command without deleting output directory. It should start from the same step as it was aborted. If you want to rerun it from the first step, remove output directory, or change output path in the configuration file.

## Output files
Several folders will be produced:

### Final_output

It is the main output of the pipeline.
- cluster_counts.tsv - tab-separated table, cluster IDs and absolute counts across all samples.
- rep_seqs.fasta - representative sequences for each cluster, corrected by "polishing" (SPOA, two rounds of Racon and one of Medaka)
- DB-taxonomy.tsv - tab-separated table, cluster IDs and taxonomy annotations (ranks by columns), read length and percent identity from blast.
- DB-taxonomy-q2.tsv - same as above, but in Qiime2 format (all ranks pooled, separated by ";" and prefixed with "r__", where r is the first character of the rank). It can be imported to qiime2.

### Quality_control
- Chopper - contains reads after QC with Chopper.

### Clustering
- Folders with sample names:
  - clusters.tsv - a table with sequence IDs and assigned clusters for given sample
  - clusters.png - scatterplot of clusters for given sample
  - kmers.tsv - a table with kmer counts for each sequence for given sample
  - subsampled_ids.tsv - sequences with kmers counts, randomly selected for "between samples" clustering
- Clusters_subsampled - directory with "between samples" subsampled fastq files, consensus sequences and IDs of subsampled representatives
- shared_clusters.tsv - selected features from each cluster and each sample, clustered "between samples"
- clusters.png - scatterplot of clusters "between samples"
- consensus_pooled.fa - fasta file with consensus sequences (SPOA) of each "between samples" cluster
- pooled.fq - pooled fastq file of all samples
  
Counts of clusters are stored in the "Final_output" folder.

### Read_correction
Contains fasta files with "best" read for each cluster, polished by two rounds of racon and one of Medaka.
These reads are merged into one fasta file in the "Final_output" folder.

### Taxonomy_annotation
- GTDB/NCBI - folder with GTDB/NCBI 16S database
- DB-blastn.tsv - output from blastn run with
- DB-blastn_tophit.tsv - tophit annotation for each cluster

For NCBI database, full taxonomies are parsed from the NCBI website based on taxid. Final taxonomy tables are stored in the "Final_output" folder.

### Logs
Log files for steps that produce logs worth reading.

## Export (Import?) to Qiime2
Files from NaMeco "Final_output" folder may be used for all kind of analyses in R or Python, but also may be imported to Qiime2 for downstream analyses:
- cluster_counts.tsv - to table.qza as feature table with frequencies (absolute counts)
- rep_seqs.fasta - to rep-seq.qza as representative sequences 
- DB-taxonomy-q2.tsv - to taxonomy.qza file

Example commands are listed below:


```python
#deactivate NaMeco (if activated)
conda deactivate

#activate qiime2 env (https://docs.qiime2.org/2024.2/install/native/)
conda activate qiime2-amplicon-2024.2

mkdir to_qiime2

#import feature table
biom convert \
    -i NaMeco_out/Final_output/cluster_counts.tsv \
    -o to_qiime2/table.biom \
    --table-type="OTU table" \
    --to-hdf5

qiime tools import \
    --input-path to_qiime2/table.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path to_qiime2/table.qza

#import taxonomy
qiime tools import \
    --type 'FeatureData[Taxonomy]' \
    --input-path NaMeco_out/Final_output/GTDB-taxonomy-q2.tsv \
    --output-path to_qiime2/taxonomy.qza

#import representative sequences
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path NaMeco_out/Final_output/rep_seqs.fasta \
    --output-path to_qiime2/rep-seq.qza
```

## Developer recommendations
All samples that are compared to each other should be run together in one pool, even from different sequencing runs. Do not merge different NaMeco runs since Cluster IDs would not match.

## Unassigned sequences
When I blasted unassigned sequences from different datasets I tested on the NCBI blastn, those sequences were annotated as host DNA. Somehow host DNA was amplified with bacterial primers. So, for downstream analyses, one should either remove unassigned sequences, or blast it on the NCBI to double check. 
