![](logo/NaMeco_logo.png)

# Processing full-length 16S Nanopore reads

## Why this pipeline was created?
16S long Nanopore reads have an advantage in terms of length compared with Illumina short reads and have the potential for better taxonomical annotations. However, in practice, due to the relatively high error rate, reads produced by Nanopore are more challenging to handle. The combination of longer length and higher error rate results in a higher number of unique reads. Moreover, clustering Nanopore reads at a 97% similarity threshold did not improve the situation, only slightly decreasing the number of unique features (sequences). 

Here, we decided to merge ideas from different exsisting tools and create one pipeline ~to rule them all~ that will be able to provide taxonomy annotations at the species level and will be easy to handle. 

So, Nameco will preprocess the reads, count kmers and then perform clustering with UMAP + HDBscan sample by sample. After that, form each cluster of each sample representatives will be randomly selected for additional clustering between samples to clsuter clusters. New clusters, this time already "shared" between samples, will be polished with combination of SPOA and Racon. Taxonomy will be assigned based on either GTDB or NCBI databases.

## Dependencies 
Linux OS with conda installed (anaconda3 or miniconda3). 
- pandas=2.2.2
- scikit-learn=1.4.2
- python=3.10.8
- chopper=0.7.0
- racon=1.5.0
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
wget https://raw.githubusercontent.com/timyerg/NaMeco/main/NaMeco.yaml
conda env create --file NaMeco.yaml
```

To update (only if needed!) the NaMeco version inside of already created environment, one can use the following command:


```python
pip install nameco --upgrade
```

Attention! The command above will only upgrade NaMeco script. If some dependencies are outdated, just delete the environment and reisnstall it.

## Running the pipeline
This pipeline takes Nanopore reads as input in fastq format. It will automatically recognize .fastq, .fq extensions. Reads also can be gziped (.fastq.gz, .fq.gz)

To run the pipeline, please provide the path to raw reads and adjust threads. The rest of the parameters may be changed if needed. 



```python
usage: nameco [-h] --inp_dir INP_DIR [--out_dir OUT_DIR]
                     [--threads THREADS] [--qc] [--no-qc] [--phred PHRED]
                     [--min_length MIN_LENGTH] [--max_length MAX_LENGTH]
                     [--kmer KMER] [--no-low] [--low]
                     [--cluster_size CLUSTER_SIZE] [--subsample SUBSAMPLE]
                     [--select_epsilon SELECT_EPSILON] [--gap GAP]
                     [--min_fraction MIN_FRACTION]
                     [--random_state RANDOM_STATE] [--database {GTDB,NCBI}]
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
  --phred PHRED         Minimum phred average score for chopper (default 8)
  --min_length MIN_LENGTH
                        Minimum read length for chopper (default 1300)
  --max_length MAX_LENGTH
                        Maximum read length for chopper (default 1700)
  --kmer KMER           K-mer length for clustering (default 5)
  --no-low              Don't restrict RAM for UMAP (default)
  --low                 Reduce RAM usage by UMAP
  --cluster_size CLUSTER_SIZE
                        Minimum cluster size for HDBscan (default 500, not <
                        100!)
  --subsample SUBSAMPLE
                        Subsample bigger than that threshold clusters for
                        consensus creation and polishing by Racon (default
                        500)
  --select_epsilon SELECT_EPSILON
                        Selection epsilon for clusters (default 0.5)
  --gap GAP             Gap between the bit score of the best hit and others,
                        that are considered with the top hit for taxonomy
                        selection (default 1)
  --min_fraction MIN_FRACTION
                        If numerous hits retained after gap filtering,
                        consensus taxon should have at least this fraction to
                        be selected. Otherwise it will be set to lower level +
                        unclassified (default 0.6)
  --random_state RANDOM_STATE
                        Random state for subsampling (default 42)
  --database {GTDB,NCBI}
                        Database for taxonomy assignment (default GTDB). Only
                        GTDB or NCBI are currently supported
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
- rep_seqs.fasta - representative sequences for each cluster, corrected by "polishing" (SPOA, two rounds of Racon and one of Medaka)
- DB-taxonomy.tsv - tab-separated table, cluster IDs and taxonomy annotations (ranks by columns), read length and percent identity from blast.
- DB-taxonomy-q2.tsv - same as above, but in Qiime2 format (all ranks pooled, separated by ";" and prefixed with "r__", where r is the first character of the rank). It can be imported to qiime2.
- DB-taxonomy-rank.tsv - collapsed to corresponding rank taxonomies with counts.

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
Contains fasta files with "best" read for each cluster, polished by Racon.
These reads are merged into one fasta file in the "Final_output" folder.

### Taxonomy_annotation
- GTDB/NCBI - folder with GTDB/NCBI 16S database
- DB-blastn.tsv - output from blastn run with

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

#activate qiime2 env (https://docs.qiime2.org/2024.5/install/native/)
conda activate qiime2-amplicon-2024.5

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

If needed, abovementioned commands can be adapted for importing collapsed taxonomy ranks. Or just collapse your feature table with clusters in qiime2 based on the taxonomy to a desired level.

## Developer recommendations
- All samples that are compared to each other should be run together in one pool, even from different sequencing runs. Do not merge different NaMeco runs at Cluster level since Cluster IDs would not match. If needed, we recommend to merge different NaMeco runs at taxonomy level.
- Adjust minimum cluster size according to your reads depth. Default 500 should work for most of the samples, but one don't have a lot of reads in a sample, then it should be decreased
- If both GTDB and NCBI annotations are needed, one can run the same command 2 times, each time indicating different database, but without deleting any files. So NaMeco in second run will reuse files from the first run. At the end, each cluster will have annotations from both databases.
- Using multiple threads can significantly speed up the Nameco run

## Unassigned sequences
When I blasted unassigned sequences from different datasets I tested on the NCBI blastn, those sequences were annotated as host DNA. Somehow host DNA was amplified with bacterial primers. So, for downstream analyses, one should either remove unassigned sequences, or blast it on the NCBI to double check. 

## Errors

- Error with "FullID": that means that you have at least one sample with low amount of reads. Try to decrease minimum cluster size or remove the sample with low amount of reads.
- TypeError("cannot pickle '_io.BufferedReader' object") with NCBI database: some issues with the NCBI website, relaunch later - it will start from the parsing step.

## Citation
If you used NaMeco tool, please cite our paper: (will be added later)
