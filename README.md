![](logo/NaMeco_logo.png)

# Processing full-length 16S Nanopore reads

## Why this pipeline was created?
16S long Nanopore reads have an advantage in terms of length compared with Illumina short reads and have the potential for better taxonomical annotations. However, in practice, due to the relatively high error rate, reads produced by Nanopore are more challenging to handle. The combination of longer length and higher error rate results in a higher number of unique reads. Moreover, clustering Nanopore reads at a 97% similarity threshold did not improve the situation, only slightly decreasing the number of unique features (sequences). 

Here, we decided to merge ideas from different exsisting tools and create one pipeline ~to rule them all~ that will be able to provide taxonomy annotations at the species level and will be easy to handle. 

So, Nameco will preprocess the reads, count kmers and then perform clustering with UMAP + HDBscan sample by sample. After that, form each cluster of each sample representatives will be randomly selected for additional clustering between samples to clsuter clusters. New clusters, this time already "shared" between samples, will be polished with combination of SPOA and Racon. Taxonomy will be assigned based on GTDB database.

## Dependencies 
Linux OS with conda installed (anaconda3, miniconda3 or miniforge). 
- qiime2 (partial)
- q2cli
- q2templates
- q2-types
- q2-longitudinal
- q2-feature-classifier
- pandas>=0.25.3
- xmltodict 
- ncbi-datasets-pyclient=18.4.0
- rescript
- python
- chopper=0.7.0
- biopython
- matplotlib
- blast=2.16
- scikit-learn
- umap-learn=0.5.7
- racon=1.5.0
- minimap2=2.28
- spoa=4.1.4
- ipykernel
- pigz
- pip
- pip:
    - nameco

## Installation
During installation, new conda environment NaMeco will be created, which includes all dependencies.

This pipeline can be installed with the following command:


```python
wget https://raw.githubusercontent.com/timyerg/NaMeco/main/NaMeco.yaml
conda env create --file NaMeco.yaml
```

To update (only if needed!) the NaMeco version inside of already created environment, one can use the following command:


```python
pip install nameco --upgrade
```

Attention! The command above will only upgrade NaMeco script. If some dependencies are outdated, just delete the environment and reisnstall it.

Hint: if you are stacked for a long time on "Solving environment" step try:


```python
conda update -n base conda
conda install -n base conda-libmamba-solver
config --set solver libmamba
```

## Running the pipeline
This pipeline takes Nanopore reads as input in fastq format. It will automatically recognize .fastq, .fq extensions. Reads also can be gziped (.fastq.gz, .fq.gz)

To run the pipeline, please provide the path to raw reads and adjust threads. The rest of the parameters may be changed if needed. 



```python
usage: nameco [-h] --inp_dir INP_DIR [--out_dir OUT_DIR] [--threads THREADS]
              [--qc] [--no-qc] [--phred PHRED] [--min_length MIN_LENGTH]
              [--max_length MAX_LENGTH] [--primer_F PRIMER_F]
              [--primer_R PRIMER_R] [--primer_PI PRIMER_PI] [--kmer KMER]
              [--cluster_size CLUSTER_SIZE] [--subsample SUBSAMPLE]
              [--select_epsilon SELECT_EPSILON] [--db_type DB_TYPE]
              [--db_version DB_VERSION] [--gap GAP]
              [--min_fraction MIN_FRACTION] [--mask_taxa] [--no_masking]
              [--random_state RANDOM_STATE] [--n_polish N_POLISH]
              [--db_path DB_PATH] [--version]

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
  --phred PHRED         Minimum phred score for chopper (default 10)
  --min_length MIN_LENGTH
                        Minimum read length for chopper (default 1300)
  --max_length MAX_LENGTH
                        Maximum read length for chopper (default 1700)
  --primer_F PRIMER_F   Forward primer (default AGAGTTTGATCMTGGCTCAG)
  --primer_R PRIMER_R   Reverse primer (default CGGTTACCTTGTTACGACTT)
  --primer_PI PRIMER_PI
                        Percent identity for primers (default 0.6)
  --kmer KMER           K-mer length for clustering (default 5)
  --cluster_size CLUSTER_SIZE
                        Min. unique cluster size (default 10, can't be < 10)
  --subsample SUBSAMPLE
                        Subsample clusters for consensus creation and
                        polishing (default 200)
  --select_epsilon SELECT_EPSILON
                        Selection epsilon for clusters (default 0.1)
  --db_type DB_TYPE     Use all rRNAs from GTDB ("All", higher accuracy,
                        slower) or only representative species ("SpeciesReps",
                        lower accuracy, faster) (default "All")
  --db_version DB_VERSION
                        GTDB version. Choices: "202.0", "207.0", "214.0",
                        "214.1", "220.0", "226.0" (default "226.0")
  --gap GAP             Gap between the bit score of the best hit and others,
                        that are considered with the top hit for taxonomy
                        selection (default 1)
  --min_fraction MIN_FRACTION
                        If numerous hits retained after gap filtering,
                        consensus taxon should have at least this fraction to
                        be selected. Otherwise set as lower level +
                        unclassified (default 0.6)
  --mask_taxa           Mask taxonomy ranks based on percent identity
                        thresholds (default "True"). Thresholds are: d: 65, p:
                        75, c: 78.5,o: 82, f: 86.5, g: 94.5, s: 97
  --no_masking          Skip masking taxonomy step
  --random_state RANDOM_STATE
                        Random state for subsampling (default 888)
  --n_polish N_POLISH   Number of polishing rounds (default 3)
  --db_path DB_PATH     Path to store/existing database (default
                        $out_dir/$database). Please use only databases,
                        created by previous NaMeco run to avoid errors
  --version             Check the version
```


```python
# example 

conda activate NaMeco
nameco --inp_dir Reads --threads 20 

#where fastq files are located in the "Reads" folder, and 20 threads are requested.
```

If the run was killed, it can be relaunched with the same command without deleting output directory. It should start from the same step as it was aborted. If you want to rerun it from the first step, remove output directory, or change output path in the configuration file.

## Output files
Several folders will be produced:

### Final_output

It is the main output of the pipeline.
- cluster_counts.tsv - tab-separated table, cluster IDs and absolute counts across all samples.
- rep_seqs.fasta - representative sequences for each cluster, corrected by "polishing" (SPOA and Racon)
- Taxonomy.tsv - tab-separated table, cluster IDs and taxonomy annotations (ranks by columns), read length and percent identity from blast.
- Taxonomy-q2.tsv - same as above, but in Qiime2 format (all ranks pooled, separated by ";" and prefixed with "r__", where r is the first character of the rank). It can be imported to qiime2.
- rank.tsv - collapsed to corresponding rank taxonomies with counts.

### Quality_control
- Chopper - contains reads after QC with Chopper.
- Fastas - fasta files after primer-specific read extraction with RESCRIPt

### Clustering
- Folders with sample names:
  - clusters.tsv - a table with sequence IDs and assigned clusters for given sample
  - kmers.tsv - a table with kmer counts for each sequence for given sample
  - subsampled_ids.tsv - sequences with kmers counts, randomly selected for "between samples" clustering
- Clusters_subsampled - directory with "between samples" subsampled fastq files, consensus sequences and IDs of subsampled representatives
- shared_clusters.tsv - selected features from each cluster and each sample, clustered "between samples"
- consensus_pooled.fa - fasta file with consensus sequences (SPOA) of each "between samples" cluster
- pooled.fq - pooled fastq file of all samples
  
Counts of clusters are stored in the "Final_output" folder.

### Read_correction
Contains fasta files with "best" read for each cluster, polished by Racon.
These reads are merged into one fasta file in the "Final_output" folder.

### Taxonomy_annotation
- GTDB - folder with GTDB SSU (16S rRNA gene) database
- blastn.tsv - output from blastn run with

### Logs
Log files for steps that produce logs worth reading.

## Export (Import?) to Qiime2
Files from NaMeco "Final_output" folder may be used for all kind of analyses in R or Python, but also may be imported to Qiime2 for downstream analyses:
- cluster_counts.tsv - to table.qza as feature table with frequencies (absolute counts)
- rep_seqs.fasta - to rep-seq.qza as representative sequences 
- Taxonomy-q2.tsv - to taxonomy.qza file

Example commands are listed below:


```python
#can be done in your qiime2 env. or in NaMeco

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
    --input-path NaMeco_out/Final_output/Taxonomy-q2.tsv \
    --output-path to_qiime2/taxonomy.qza

#import representative sequences
qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path NaMeco_out/Final_output/rep_seqs.fasta \
    --output-path to_qiime2/rep-seq.qza
```

If needed, abovementioned commands can be adapted for importing collapsed taxonomy ranks. Or just collapse your feature table with clusters in qiime2 based on the taxonomy to a desired level.

## Developer recommendations
- All samples that are compared to each other should be run together in one pool, even from different sequencing runs. Do not merge different NaMeco runs at Cluster level since Cluster IDs would not match. If needed, we recommend to merge different NaMeco runs at taxonomy level.
- Using multiple threads can significantly speed up the NaMeco run
- If you are facing issues with the drive space on your working drive, export tmp directory before running NaMeco: "export TMPDIR=/big_storage_path/TMP"
- "Creating database" can take a long time, just wait until it is done. Good time to make some coffee!

## Unassigned sequences
- When I blasted unassigned sequences from different datasets I tested on the NCBI blastn, those sequences were annotated as host DNA. Somehow host DNA was amplified with bacterial primers. So, for downstream analyses, one should either remove unassigned sequences, or blast it on the NCBI to double check.

## Errors

- Give the feedback when encountered by creating an issue.

## Citation
If you used NaMeco tool, please cite our paper:


```python
Yergaliyev, T., Rios-Galicia, B. & Camarinha-Silva, A. 
NaMeco - Nanopore full-length 16S rRNA gene reads clustering and annotation. 
BMC Genomics (2025). https://doi.org/10.1186/s12864-025-12415-x
```
