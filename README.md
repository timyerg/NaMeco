{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "7c215dd6-9f1b-4f63-8c70-18170ed0c262",
   "metadata": {},
   "source": [
    "# NaMeco\n",
    "This pipeline was written for processing full-length 16S Nanopore reads"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "c4e0e850-2023-4545-86da-49e311d5f79f",
   "metadata": {},
   "source": [
    "## Citation\n",
    "If you used our pipelines, please cite this paper: XXX"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "de2e9886-a7c1-424a-9474-25c588027195",
   "metadata": {},
   "source": [
    "## Why this pipeline was created?\n",
    "16S long Nanopore reads have an advantage in terms of length compared with Illumina short reads and have the potential for better taxonomical annotations. However, in practice, due to the relatively high error rate, reads produced by Nanopore are more challenging to handle. The combination of longer length and higher error rate results in a higher number of unique reads. When we tried to use blast on unclustered Nanopore reads, a dataset with 170 samples was not finished on the HPC with 24 threads in a month. However, clustering Nanopore reads at a 97% similarity threshold did not improve the situation, only slightly decreasing the number of unique features (sequences). \n",
    "Here, we decided to merge ideas from such pipelines as NanoCLUST, NGSpeciesID and Natrix2 and create one pipeline ~to rule them all~ that will combine the advantages of above mentioned tools. So, it will preprocess the reads, count kmers and then perform clustering with UMAP + HDBscan sample by sample. After that, form each cluster of each sample 50 representatives will be randomly selected for additional clustering between samples to clsuter clusters. New clusters, this time already \"shared\" between samples, will be polished with combination of SPOA, Racon and Medaka. Taxonomy will be assigned based on either GTDB or NCBI databases."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "5dc7c58d-6fc1-4f01-885a-961acde47203",
   "metadata": {},
   "source": [
    "## Dependencies \n",
    "Linux OS with conda installed (anaconda3 or miniconda3). \n",
    "- pandas=2.2.2\n",
    "- scikit-learn=1.4.2\n",
    "- python=3.10.8\n",
    "- chopper=0.7.0\n",
    "- racon=1.5.0\n",
    "- medaka=1.11.3\n",
    "- minimap2=2.28\n",
    "- umap-learn=0.5.5\n",
    "- biopython=1.83\n",
    "- matplotlib=3.8.4\n",
    "- blast=2.15\n",
    "- spoa=4.1.4"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "9528bde9-d409-4365-a382-6545bdd6e601",
   "metadata": {},
   "source": [
    "## Installation\n",
    "During installation, new conda environment NaMeco will be created, which includes all dependencies.\n",
    "\n",
    "This pipeline can be installed with the following script:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "77b6fdbd-2e9e-4583-8348-c3615af0f23f",
   "metadata": {},
   "outputs": [],
   "source": [
    "wget ... #edit\n",
    "conda env create --file ... #edit"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "d4cd008d-9542-44b2-be34-7284c64097f8",
   "metadata": {},
   "source": [
    "## Running the pipeline\n",
    "This pipeline takes Nanopore reads as input in fastq format. It will automatically recognize .fastq, .fq extensions. Reads also can be gziped (.fastq.gz, .fq.gz)\n",
    "\n",
    "To run the pipeline, please provide the path to raw reads and adjust threads. The rest of the parameters may be changed if needed. \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "15c234e3-d1e3-4cef-a1bc-2df7d916f9c5",
   "metadata": {},
   "outputs": [],
   "source": [
    "nameco --help"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "442b72ea-c2d3-4df0-8ac4-51cc35a5f912",
   "metadata": {},
   "source": [
    "If the run was killed, it can be relaunched with the same command without deleting output directory. It should start from the same step as it was aborted. If you want to rerun it from the first step, remove output directory, or change output path in the configuration file."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "594e8bfd-a1c7-406d-a0a5-d3d014194a58",
   "metadata": {},
   "source": [
    "## Output files\n",
    "Several folders will be produced:"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "07af4082-8e5e-4164-b760-38c7a6e7ce40",
   "metadata": {},
   "source": [
    "### Final_output"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "f89ef571-306a-48bb-aefa-8a25686c533f",
   "metadata": {},
   "source": [
    "It is the main output of the pipeline.\n",
    "- cluster_counts.tsv - tab-separated table, cluster IDs and absolute counts across all samples.\n",
    "- rep_seqs.fasta - representative sequences for each cluster, corrected by \"polishing\" (SPOA, two rounds of Racon and one of Medaka)\n",
    "- DB-taxonomy.tsv - tab-separated table, cluster IDs and taxonomy annotations (ranks by columns), read length and percent identity from blast.\n",
    "- DB-taxonomy-q2.tsv - same as above, but in Qiime2 format (all ranks pooled, separated by \";\" and prefixed with \"r__\", where r is the first character of the rank). It can be imported to qiime2."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "ba8ba41d-2492-4056-a9e7-279a364b68fe",
   "metadata": {},
   "source": [
    "### Quality_control\n",
    "- Chopper - contains reads after QC with Chopper."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "633aff95-944b-4f68-8df4-4562146ed142",
   "metadata": {},
   "source": [
    "### Clustering\n",
    "- Folders with sample names:\n",
    "  - clusters.tsv - a table with sequence IDs and assigned clusters for given sample\n",
    "  - clusters.png - scatterplot of clusters for given sample\n",
    "  - kmers.tsv - a table with kmer counts for each sequence for given sample\n",
    "  - subsampled_ids.tsv - sequences with kmers counts, randomly selected for \"between samples\" clustering\n",
    "- Clusters_subsampled - directory with \"between samples\" subsampled fastq files, consensus sequences and IDs of subsampled representatives\n",
    "- shared_clusters.tsv - selected features from each cluster and each sample, clustered \"between samples\"\n",
    "- clusters.png - scatterplot of clusters \"between samples\"\n",
    "- consensus_pooled.fa - fasta file with consensus sequences (SPOA) of each \"between samples\" cluster\n",
    "- pooled.fq - pooled fastq file of all samples\n",
    "  \n",
    "Counts of clusters are stored in the \"Final_output\" folder."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "4cb12579-2ab4-476a-9033-c7475a55f6fa",
   "metadata": {},
   "source": [
    "### Read_correction\n",
    "Contains fasta files with \"best\" read for each cluster, polished by two rounds of racon and one of Medaka.\n",
    "These reads are merged into one fasta file in the \"Final_output\" folder."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "2b6e1f90-e8c4-42a0-b4b6-7fbcbe4f716e",
   "metadata": {},
   "source": [
    "### Taxonomy_annotation\n",
    "- GTDB/NCBI - folder with GTDB/NCBI 16S database\n",
    "- DB-blastn.tsv - output from blastn run with\n",
    "- DB-blastn_tophit.tsv - tophit annotation for each cluster\n",
    "\n",
    "For NCBI database, full taxonomies are parsed from the NCBI website based on taxid. Final taxonomy tables are stored in the \"Final_output\" folder."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "d16e8bb1-3e89-46ea-b292-b2aa140a4b79",
   "metadata": {},
   "source": [
    "### Logs\n",
    "Log files for steps that produce logs worth reading."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "aac6f51b-7c74-4916-98bb-723a5097841b",
   "metadata": {},
   "source": [
    "## Export (Import?) to Qiime2\n",
    "Files from NaMeco \"Final_output\" folder may be used for all kind of analyses in R or Python, but also may be imported to Qiime2 for downstream analyses:\n",
    "- cluster_counts.tsv - to table.qza as feature table with frequencies (absolute counts)\n",
    "- rep_seqs.fasta - to rep-seq.qza as representative sequences \n",
    "- DB-taxonomy-q2.tsv - to taxonomy.qza file\n",
    "\n",
    "Example commands are listed below:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "5eb15fa1-ff18-4397-8072-b37533f64082",
   "metadata": {},
   "outputs": [],
   "source": [
    "#deactivate NaMeco (if activated)\n",
    "conda deactivate\n",
    "\n",
    "#activate qiime2 env (https://docs.qiime2.org/2024.2/install/native/)\n",
    "conda activate qiime2-amplicon-2024.2\n",
    "\n",
    "mkdir to_qiime2\n",
    "\n",
    "#import feature table\n",
    "biom convert \\\n",
    "    -i NaMeco_out/Final_output/cluster_counts.tsv \\\n",
    "    -o to_qiime2/table.biom \\\n",
    "    --table-type=\"OTU table\" \\\n",
    "    --to-hdf5\n",
    "\n",
    "qiime tools import \\\n",
    "    --input-path to_qiime2/table.biom \\\n",
    "    --type 'FeatureTable[Frequency]' \\\n",
    "    --input-format BIOMV210Format \\\n",
    "    --output-path to_qiime2/table.qza\n",
    "\n",
    "#import taxonomy\n",
    "qiime tools import \\\n",
    "    --type 'FeatureData[Taxonomy]' \\\n",
    "    --input-path NaMeco_out/Final_output/GTDB-taxonomy-q2.tsv \\\n",
    "    --output-path to_qiime2/taxonomy.qza\n",
    "\n",
    "#import representative sequences\n",
    "qiime tools import \\\n",
    "    --type 'FeatureData[Sequence]' \\\n",
    "    --input-path NaMeco_out/Final_output/rep_seqs.fasta \\\n",
    "    --output-path to_qiime2/rep-seq.qza"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "a5cca725-fc8a-469f-b31a-0502c3184ad8",
   "metadata": {},
   "source": [
    "## Developer recommendations\n",
    "All samples that are compared to each other should be run together in one pool, even from different sequencing runs. Do not merge different NaMeco runs since Cluster IDs would not match."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "a2537372-991b-4ce2-bc97-888bc5891169",
   "metadata": {},
   "source": [
    "## Unassigned sequences\n",
    "When I blasted unassigned sequences from different datasets I tested on the NCBI blastn, those sequences were annotated as host DNA. Somehow host DNA was amplified with bacterial primers. So, for downstream analyses, one should either remove unassigned sequences, or blast it on the NCBI to double check. "
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python [conda env:qiime2-amplicon-2024.2]",
   "language": "python",
   "name": "conda-env-qiime2-amplicon-2024.2-py"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.15"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
