import os
import pandas as pd
from ._helper import *

#Function for read correction (Racon)
def read_correction(T, N, OUT, FI, log):
    print(f'Polishing with Racon...')
    bash(f'mkdir -p {OUT}')
    file = f'{OUT}/../Clustering/consensus_pooled.fa'
    corr = f'{FI}/rep_seqs.fasta'
    df = pd.read_csv(f"{OUT}/../Final_output/cluster_counts.tsv", sep='\t', index_col=0)
    clusters = df.index.tolist()
    skip, checks = log_checker(log, ['mock'], f'{OUT}/{{}}/mock')
    if os.path.exists(corr):
        if os.stat(corr).st_size > 0 and os.stat(file).st_size > 0:
            if int(bash(f'grep -c "^>" {file}')) == int(bash(f'grep -c "^>" {corr}')):
                return print(f'Consensuses for all clusters were already corrected. Skipping')
    for cluster in clusters:
        if os.path.exists(f"{OUT}/{cluster}_racon{N-1}.fa") and cluster in skip:
            continue
        fa = f'{OUT}/../Clustering/Clusters_subsampled/{cluster}_consensus.fa'
        sam = f"{OUT}/{cluster}.sam"
        fq = f'{OUT}/../Clustering/Clusters_subsampled/{cluster}.fq.gz'
        bash(f'echo "\n##### Processing {cluster} #####" >> {log}')
        for n in range(N):
            po = f"{OUT}/{cluster}_racon{n}.fa"
            ta = f"{OUT}/{cluster}_racon{n-1}.fa"
            if n == 0:
                ta = fa
            #overlaping with minimap2
            bash(f'echo "\nMapping {cluster} {n}" >> {log}')
            bash(f"minimap2 -ax map-ont -t {T} {ta} {fq} -o {sam} 2>> {log}")
            #polishing with racon
            bash(f'echo "\nPolishing with Racon {cluster}" >> {log}')
            bash(f'racon -m 8 -x -6 -g -8 -t {T} {fq} {sam} {ta} > {po} 2>> {log}')
            bash(f'rm {sam}')
            bash(f'echo "{cluster} done. Enjoy" >> {log}')
    #collect corrected sequences
    with open(corr, "w") as corrected:
        for cluster in clusters:
            with open(f"{OUT}/{cluster}_racon{N-1}.fa", "rt") as rep:
                corrected.write(">{id}\n{seq}\n".format(id=cluster, seq=rep.read().split('\n')[1],))

