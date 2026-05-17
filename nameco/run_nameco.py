#! /usr/bin/env python

import os
import re
import argparse
import subprocess
import numpy as np
import pandas as pd
from importlib.metadata import version
from ._helper import *
from ._qc import *
from ._clustering import *
from ._polishing import *
from ._taxonomy_annotation import *


#Function to run NaMeco 
def main():
    ##### ARGPARSE #####
    inp_dir_help = " ".join(['Path to the folder with reads, absolute or relative.', 
                             'Reads should be in the fastq or fq format, gziped or not'])
    out_dir_help = " ".join(['Path to the directory to store output files, absolute or relative.', 
                             'If not provided, folder "Nameco_out" will be created in working directory'])
    gap_help = " ".join(['Gap between the bit score of the best hit and others,',
                         'that are considered with the top hit for taxonomy selection (default 1)'])
    frac_help = " ".join(['If numerous hits retained after gap filtering, consensus taxon should have at least this',
                          'fraction to be selected. Otherwise set as lower level + unclassified (default 0.6)'])
    fetch_db_help = " ".join(['Fetch prebuild database. Choices: "GTDB_220", "GTDB_226", "GTDB_232",',
                              '"UNITE_fungi_V10", "UNITE_fungi-2_V10", "UNITE_eukaryotes_V10"',
                              '"UNITE_eukaryotes-2_V10" (default False)'])
    db_version_help = " ".join(['GTDB version. Choices: "220.0", "226.0", "232.0", "latest" (default "latest")'])
    mask_taxa_help = " ".join(['Mask taxonomy ranks based on percent identity thresholds (default "True").',
                               'Thresholds are: d: 65, p: 75, c: 78.5,o: 82, f: 86.5, g: 94.5, s: 97'])
    db_path_help = " ".join(['Path to store/existing database (default $out_dir/$database).', 
                             'Please use only databases, created by previous NaMeco run to avoid errors'])

    parser = argparse.ArgumentParser(prog='nameco')
    parser._action_groups.pop()
    req = parser.add_argument_group('required arguments')
    opt = parser.add_argument_group('optional arguments')
    req.add_argument("--inp_dir", help=inp_dir_help, required=True)
    opt.add_argument("--out_dir", help=out_dir_help, default='NaMeco_out')
    opt.add_argument("--threads", help="The number of threads/cpus (default 2)", type=int, default=2)
    opt.add_argument("--qc", help="Run chopper for quality control (default)", action='store_true', default=True)
    opt.add_argument("--no-qc", help="Skip chopper for quality control", dest='qc', action='store_false')
    opt.add_argument("--phred", help="Minimum phred score for chopper (default 10)", type=int, default=10)
    opt.add_argument("--min_length", help="Minimum read length for chopper (default 1300)", type=int, default=1200)
    opt.add_argument("--max_length", help="Maximum read length for chopper (default 1700)", type=int, default=2000)
    opt.add_argument("--min_sample_size", help="Minimum sample size to be retained (default 500)", type=int, default=500)
    opt.add_argument("--kmer", help="K-mer length for clustering (default 5)", type=int, default=5)
    opt.add_argument("--cluster_size", help="Min. unique cluster size (default 10, can't be < 10)", type=int, default=10)
    opt.add_argument("--subsample", help='Subsample clusters for consensus creation and polishing (default 200)', type=int, default=200)
    opt.add_argument("--select_epsilon", help="Selection epsilon for clusters (default 0.1)", type=float, default=0.1)
    opt.add_argument('--fetch_db', help=fetch_db_help, default=False)
    opt.add_argument('--db_version', help=db_version_help, default='latest')
    opt.add_argument("--gap", help=gap_help, type=float, default=1)
    opt.add_argument("--min_fraction", help=frac_help, type=float, default=.6)
    opt.add_argument("--mask_taxa", help=mask_taxa_help, action='store_true', default=True)
    opt.add_argument("--no_masking", help="Skip masking taxonomy step", dest='mask_taxa', action='store_false')
    opt.add_argument("--random_state", help="Random state for subsampling (default 888)", type=int, default=888)
    opt.add_argument("--n_polish", help="Number of polishing rounds (default 3)", type=int, default=3)
    opt.add_argument('--db_path', help=db_path_help, default='{OUT}/{DB}-{DBV}')
    opt.add_argument('--version', help="Check the version", action="version", version=version("nameco"))
    args = parser.parse_args()
    
    INPDIR = args.inp_dir
    LOGS = f'{args.out_dir}/Logs'
    QC = f'{args.out_dir}/Quality_control'
    CL = f'{args.out_dir}/Clustering'
    RC = f'{args.out_dir}/Read_correction'
    TA = f'{args.out_dir}/Taxonomy_annotation'
    FI = f'{args.out_dir}/Final_output'
    exts = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
    greetings(log=f"{LOGS}/NaMeco_version.log", LOGS=LOGS)
    SAMPLES = [f.split('.')[0] for f in os.listdir(INPDIR) if f.endswith(exts)]
    print('Only "*.fastq.gz", "*.fq.gz", "*.fastq" or "*.fq" files will be procsessed')
    if len(SAMPLES) == 0:
        raise ValueError('Input directory does not contain fastq.gz or fq.gz files')
        
    #####  Quality control ######
    module = QC.split('/')[-1]
    hashtags_wrapper(f"{module.replace('_', ' ')} module")
    log = f"{LOGS}/{module}.log"
    #chopping
    if args.qc:
        chopper(INPUT=INPDIR, SAMPLES=SAMPLES, T=args.threads, LOGS=LOGS, log=log, Q=args.phred, 
                MINL=args.min_length, MAXL=args.max_length, OUT=f'{QC}/Chopper')
        INPDIR = f'{QC}/Chopper'
        print('\nPlease, cite chopper: https://doi.org/10.1093/bioinformatics/btad311')
    else:
        print(f"Chopper is disabled. Skipping")
    #ingnoring samples with low count
    LOWS = min_sample_size(INPUT=INPDIR, SAMPLES=SAMPLES, LOGS=LOGS, log=log, MinSS=args.min_sample_size)
    SAMPLES = [s for s in SAMPLES if s not in LOWS]
    
    ###### Clustering #####
    module = CL.split('/')[-1]
    hashtags_wrapper(f"{module.replace('_', ' ')} module")
    log = f"{LOGS}/{module}.log"
    #kmers counting
    print(f"Counting kmers ({args.kmer}-mers) for all samples...")
    kmer_counter(OUT=CL, INPUT=INPDIR, SAMPLES=SAMPLES, T=args.threads, L=args.kmer, log=log)
    #clustering with UMAP + HDBscan
    clustering_UMAP_HDBscan(OUT=CL, T=args.threads, EPS=args.select_epsilon, log=log,
                            CLUST_UQ=args.cluster_size, SAMPLES=SAMPLES, RSTAT=args.random_state)
    #pool clusters from samples to shared clusters and recalculate abundances
    shared_clusters(OUT=CL, FI=FI, SAMPLES=SAMPLES, RSTAT=args.random_state, 
                    SUBS=args.subsample, T=args.threads, log=log)
    #spliting fasta by cluster
    file_by_cluster(INPUT=INPDIR, subs=args.subsample, OUT=CL, T=args.threads, log=log)
    print('\nPlease, cite UMAP: https://doi.org/10.21105/joss.00861')
    print('Please, cite HDBscan: https://doi.org/10.21105/joss.00205')
    print('Please, cite SPOA: https://doi.org/10.1101%2Fgr.214270.116')
    print(f"\nEnd of the {module.replace('_', ' ')} module")
    
    ###### Read correction #####
    module = RC.split('/')[-1]
    hashtags_wrapper(f"{module.replace('_', ' ')} module")
    log = f"{LOGS}/{module}.log"
    read_correction(OUT=RC, N=args.n_polish, FI=FI, T=args.threads, log=log)
    print('\nPlease, cite minimap2: https://doi.org/10.1093/bioinformatics/bty191')
    print('Please, cite racon: https://doi.org/10.1101%2Fgr.214270.116')
    print(f"\nEnd of the {module.replace('_', ' ')} module")

    ###### Taxonomy annotation ######
    module = TA.split('/')[-1]
    hashtags_wrapper(f"{module.replace('_', ' ')} module")
    log = f"{LOGS}/{module}.log"
    DBpath=args.db_path
    FETCH=args.fetch_db
    taxonomy_annotation(DB='GTDB', DBV=args.db_version, gap=args.gap, frac=args.min_fraction, 
                        T=args.threads, OUT=TA, FI=FI, DBpath=DBpath, MASK=args.mask_taxa, 
                        FETCH=FETCH, log=log)
    if "GTDB" in DBpath:
        print('\nPlease, cite GTDB database: https://doi.org/10.1038/s41587-020-0501-8')
    if "UNITE" in DBpath or "UNITE" in FETCH:
        print('\nPlease, cite UNITE database paper: https://doi.org/10.1093/nar/gkad1039')
        print('Also, cite the version of UNITE database you used.')
        print('Check citation here (General FASTA release): https://unite.ut.ee/repository.php#panel5a')
    else:
        print('\nPlease, cite the database you used.')
    print('\nPlease, cite BLAST: https://doi.org/10.1016/s0022-2836(05)80360-2')
    print(f"\nEnd of the {module.replace('_', ' ')} module")
    module = "NaMeco run successfully completed. Enjoy your data!"
    hashtags_wrapper(f"{module.replace('_', ' ')}")

if __name__ == '__main__':
    main()
