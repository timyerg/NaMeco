![](logo/NaMeco_logo.png)

# Processing full-length 16S Nanopore reads

## Why this pipeline was created?
16S and ITS long Nanopore reads have an advantage in length compared with Illumina short reads and have the potential for better taxonomic annotations. However, in practice, due to their relatively high error rate, Nanopore reads are more challenging to handle. The combination of longer length and a higher error rate results in more unique reads. Moreover, clustering Nanopore reads at a 97% similarity threshold did not improve the situation, only slightly decreasing the number of unique features (sequences).

Here, we decided to merge ideas from different existing tools and create a single pipeline to rule them all that can provide species-level taxonomy annotations and is easy to handle.

So, Nameco will preprocess the reads, count kmers, and then perform clustering with UMAP + HDBscan, sample by sample. After that, from each cluster of each sample, representatives will be randomly selected for additional clustering between samples to cluster clusters. New clusters, this time already "shared" between samples, will be polished with a combination of SPOA and Racon. Taxonomy will be assigned based on the GTDB or UNITE databases.


## Dependencies 
Linux OS with conda installed (anaconda3, miniconda3 or miniforge). 
- numpy<=2
- pandas>=0.25.3
- xmltodict 
- python
- chopper=0.7.0
- biopython
- matplotlib
- blast=2.16
- scikit-learn=1.8
- umap-learn=0.5.11
- racon=1.5.0
- minimap2=2.28
- spoa=4.1.4
- ipykernel
- pigz
- pip
- pip:
    - nameco

## Installation
During installation, a new conda environment named NaMeco will be created, including all dependencies.

This pipeline can be installed with the following command:


```python
wget https://raw.githubusercontent.com/timyerg/NaMeco/main/NaMeco.yaml
conda env create --file NaMeco.yaml
```

To update (only if needed!) the NaMeco version inside of an already created environment, one can use the following command:


```python
pip install nameco --upgrade
```

Attention! The command above will only upgrade the NaMeco script. If some dependencies are outdated, just delete the environment and reinstall it.

Hint: if you are stuck for a long time on the "Solving environment" step, try:


```python
conda update -n base conda
conda install -n base conda-libmamba-solver
config --set solver libmamba
```

## Data preparation
Nanopore sequencers often output samples as multiple small-sized fastq files for each barcode. ***NaMeco accepts ONE fastq file as ONE sample***, so please ***concatenate*** multiple files that belong to the same sample ***before*** running NaMeco. For example, if you have a folder named "barcode01",  which contains multiple fastq.gz files that belong to Sample1, you can use the following command: 


```python
cat barcode01/*fastq.gz > Sample1.fastq.gz
```

You can iterate through all your samples in the loop. But make sure you are not concatenating files from different samples.

## Running the pipeline
This pipeline takes Nanopore reads as input in fastq format. It will automatically recognize .fastq, .fq extensions. Reads can also be gziped (.fastq.gz, .fq.gz)

To run the pipeline, please provide the path to the raw reads and adjust the number of threads. The remaining parameters may be adjusted as needed.


```python
usage: nameco [-h] --inp_dir INP_DIR [--out_dir OUT_DIR] [--threads THREADS]
              [--qc] [--no-qc] [--phred PHRED] [--min_length MIN_LENGTH]
              [--max_length MAX_LENGTH] [--min_sample_size MIN_SAMPLE_SIZE]
              [--kmer KMER] [--cluster_size CLUSTER_SIZE]
              [--subsample SUBSAMPLE] [--select_epsilon SELECT_EPSILON]
              [--fetch_db FETCH_DB] [--db_version DB_VERSION] [--gap GAP]
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
  --min_sample_size MIN_SAMPLE_SIZE
                        Minimum sample size to be retained (default 500)
  --kmer KMER           K-mer length for clustering (default 5)
  --cluster_size CLUSTER_SIZE
                        Min. unique cluster size (default 10, can't be < 10)
  --subsample SUBSAMPLE
                        Subsample clusters for consensus creation and
                        polishing (default 200)
  --select_epsilon SELECT_EPSILON
                        Selection epsilon for clusters (default 0.1)
  --fetch_db FETCH_DB   Fetch prebuild database. Choices: "GTDB_220",
                        "GTDB_226", "GTDB_232", "UNITE_fungi_V10",
                        "UNITE_fungi-2_V10", "UNITE_eukaryotes_V10"
                        "UNITE_eukaryotes-2_V10" (default False)
  --db_version DB_VERSION
                        GTDB version. Choices: "220.0", "226.0", "232.0",
                        "latest" (default "latest")
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

By default, the NaMeco tool will download and install the latest GTDB database. However, now you can fetch prebuilt GTDB or UNITE databases by using the “--fetch_db” option.


```python
# example 

conda activate NaMeco
nameco --inp_dir Reads --threads 20 

#where fastq files are located in the "Reads" folder, and 20 threads are requested.
```

If the run was killed, it can be relaunched with the same command without deleting the output directory. It should start from the same step as it was aborted. If you want to rerun it from the first step, remove the output directory, or change the output path in the configuration file.

## Working with ITS (UNITE database)
Starting with NaMeco v.1.4.0, we added support for ITS sequences and the UNITE database.
To run NaMeco on ITS sequences:
- Add option “--fetch_db” with corresponding UNITE database (check help message or “Databases” directory in this repository.
- Adjust the minimum and maximum lengths for Chopper based on the expected amplicon size (250-1300?).

You can check available for fetching prebuilt databases in the "Databases" directory or by calling "nameco --help". To get additional information regarding UNITE databases, check their [website.](https://unite.ut.ee/repository.php)

### Attention!
UNITE databases are licensed under [CC BY 4.0.](https://creativecommons.org/licenses/by/4.0/) Please cite databases if you used them in your pipeline (as well as NaMeco and all other tools within it! See corresponding information at the bottom of this page).


```python
# example of running with ITS reads and UNITE database

conda activate NaMeco
!nameco \
    --inp_dir tests/Samples_ITS  \
    --max_length 1500  \
    --min_length 250 \
    --threads 20 \
    --fetch_db "UNITE_fungi-2_V10"

#where fastq files are located in the "tests/Samples_ITS" folder, 
#while requesting 20 threads, setting lengths of sequences to 250-1500,
#and fetching database "UNITE_fungi-2_V10"
```

## Support of custom databases
NaMeco can be used with custom databases, though we can't guarantee that it will work.
To create a custom database, follow the steps:
- Prepare a FASTA file from the database, and put it in the new folder. Rename the fasta file as “db.fa”.
- Create a database with blast within the NaMeco environment (to make sure that the blast version is the same) with the following command: “makeblastdb -in {DBpath}/db.fa -parse_seqids -dbtype ‘nucl’”
- Create a “tab-separated” file named “map.tsv”. It should contain only 2 columns: “SeqID” and “Taxonomy”. “SeqID” column should contain sequence IDs from the fasta file, and “Taxonomy” column - a string with the taxonomy of the corresponding fasta sequence. Important: string should start from “d__” (not “k__”!) and contain ranks from domain to species, with “__” between rank identifier and taxonomy rank. Example: “d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli”
- Run Nameco with “--db_path” option pointing to your custom database. It will override the default database, forcing NaMeco to use your custom database. 

## Output files
Several folders will be produced:

### Final_output

It is the main output of the pipeline.
- cluster_counts.tsv - tab-separated table, cluster IDs, and absolute counts across all samples.
- rep_seqs.fasta - representative sequences for each cluster, corrected by "polishing" (SPOA and Racon)
- Taxonomy.tsv - tab-separated table, cluster IDs and taxonomy annotations (ranks by columns), read length, and percent identity from blast.
- Taxonomy-q2.tsv - same as above, but in Qiime2 format (all ranks pooled, separated by ";" and prefixed with "r__", where r is the first character of the rank). It can be imported to Qiime2.
- rank.tsv - collapsed to corresponding rank taxonomies with counts.

### Quality_control
- Chopper - contains reads after QC with Chopper.

### Clustering
- Folders with sample names:
  - clusters.tsv - a table with sequence IDs and assigned clusters for a given sample
  - kmers.tsv - a table with kmer counts for each sequence for a given sample
  - subsampled_ids.tsv - sequences with kmers counts, randomly selected for "between samples" clustering
- Clusters_subsampled - directory with "between samples" subsampled fastq files, consensus sequences, and IDs of subsampled representatives
- shared_clusters.tsv - selected features from each cluster and each sample, clustered "between samples"
- consensus_pooled.fa - fasta file with consensus sequences (SPOA) of each "between samples" cluster
- pooled.fq - pooled fastq file of all samples
  
Counts of clusters are stored in the "Final_output" folder.

### Read_correction
Contains FASTA files with the "best" read for each cluster, polished with Racon.
These reads are merged into one fasta file in the "Final_output" folder.

### Taxonomy_annotation
- Folder with GTDB, UNITE or other database
- blastn.tsv - output from blastn

### Logs
Log files for steps that produce logs worth reading.

## Export (Import?) to Qiime2
Files from NaMeco "Final_output" folder may be used for all kinds of analyses in R or Python, but also may be imported to Qiime2 for downstream analyses:
- cluster_counts.tsv - to table.qza as feature table with frequencies (absolute counts)
- rep_seqs.fasta - to rep-seq.qza as representative sequences 
- Taxonomy-q2.tsv - to taxonomy.qza file

Example commands are listed below:


```python
#can be done in your qiime2 env

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

If needed, the aforementioned commands can be adapted to import collapsed taxonomy ranks. Or just collapse your feature table with clusters in Qiime2 based on the taxonomy to a desired level.

## Developer recommendations
- All samples that are compared to each other should be run together in one pool, even from different sequencing runs. Do not merge different NaMeco runs at the cluster level since Cluster IDs would not match. If needed, we recommend merging different NaMeco runs at the taxonomy level.
- Using multiple threads can significantly speed up the NaMeco run.
- If you are facing issues with the drive space on your working drive, export the tmp directory before running NaMeco: "export TMPDIR=/big_storage_path/TMP".
- "Counting k-mers" can take a long time, just wait until it is done. Good time to make some coffee!

## Unassigned sequences
- When I blasted unassigned sequences from different datasets, I tested on the NCBI blastn, and those sequences were annotated as host DNA. Somehow, the host DNA was amplified with bacterial primers. So, for downstream analyses, one should either remove unassigned sequences or BLAST them against NCBI to double-check.

## Errors

- "KeyError: 'FullID'". To fix this error, delete samples with very low read counts.
- Give feedback when encountered by creating an issue.

## Citation

### NaMeco
If you used NaMeco tool, please cite our paper:
- Yergaliyev, T., Rios-Galicia, B. & Camarinha-Silva, A. 
NaMeco - Nanopore full-length 16S rRNA gene reads clustering and annotation. 
BMC Genomics (2025). https://doi.org/10.1186/s12864-025-12415-x

### Other tools we strongly recommend to cite when using NaMeco:
#### Quality control
- Chopper: https://doi.org/10.1093/bioinformatics/btad311
#### Clustering
- UMAP: https://doi.org/10.21105/joss.00861
- HDBscan: https://doi.org/10.21105/joss.00205
#### Consensus sequence
- SPOA: https://doi.org/10.1101%2Fgr.214270.116
#### Polishing
- Minimap2: https://doi.org/10.1093/bioinformatics/bty191
- Racon: https://doi.org/10.1101%2Fgr.214270.116
#### Taxonomy annotation tool
- BLAST: https://doi.org/10.1016/s0022-2836(05)80360-2
#### Database:
- If you used GTDB: https://doi.org/10.1038/s41587-020-0501-8
- If you used UNITE database, cite both their paper and database: 
    - Paper: https://doi.org/10.1093/nar/gkad1039
    - Database: Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2025): UNITE general FASTA release for {DATABASE}. Version 19.02.2025. UNITE Community. {DOI}. ##### Replace {DATABASE} and {DOI} by dorresponding data from https://unite.ut.ee/repository.php#panel5a (General FASTA release) #####.


