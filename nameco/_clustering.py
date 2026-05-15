import os
import re
import glob
import gzip
import random
import numpy as np
import pandas as pd
import multiprocessing as mp
import sklearn.cluster as cluster
from Bio import SeqIO
from itertools import product
from ._helper import *


###Functions to count kmers of given length. 

#working function
def kmer_subcounter(kmers, rec, q):
    count = [str(len(re.findall(f'(?={mer})', str(rec.seq)))) for mer in kmers]
    res = '\t'.join([rec.id] + count)
    q.put(res)


#listening fuction
def kmer_writer(q, out, kmers):
    with open (f'{out}/kmers.tsv', 'wt') as tab:
        tab.write('\t'.join(['ID'] + kmers) + '\n')
        while 1:
            m = q.get()
            if m == 'kill':
                break
            tab.write(m + '\n')
            tab.flush()


#Function to count kmers by sample
def kmer_counter(OUT, INPUT, SAMPLES, L, T, log):    
    skip, checks = log_checker(log, SAMPLES, f'{OUT}/{{}}/kmers.tsv')
    if all(checks):
        return print(f'Kmers were already counted. Skipping')
    #get all kmers possible
    nucleotides = 'ACGT'
    kmers = product(nucleotides, repeat=L)
    kmers = [''.join(c) for c in kmers]
    for sample in SAMPLES:
        file = glob.glob(f"{INPUT}/{sample}.f*q*")[0]
        out = f'{OUT}/{sample}'
        bash(f'mkdir -p {out}')
        if os.path.exists(f'{out}/kmers.tsv') and sample in skip:
            continue
        #must use Manager queue here, or will not work
        q = mp.Manager().Queue()  
        pool = mp.Pool(T)
        #put listener to work first
        watcher = pool.apply_async(kmer_writer, (q, out, kmers))
        #fire off workers
        bash(f'echo "\n##### Processing {sample} #####\n" >> {log}')
        jobs = []
        with gzip.open(file, 'rt') as f:
            for rec in SeqIO.parse(f, 'fastq'):
                job = pool.apply_async(kmer_subcounter, (kmers, rec, q))
                jobs.append(job)
        [job.get() for job in jobs]
        q.put('kill')
        pool.close()
        pool.join()
        bash(f'echo "\n{sample} done. Enjoy" >> {log}')
    bash(f'echo "\nK-mers counted." >> {log}')
    
    
###Functions for clustering
    
#Function to cluster with UMAP + HDBscan
def clustering_UMAP_HDBscan(OUT, SAMPLES, T, EPS, CLUST_UQ, RSTAT, log,):
    import umap
    print('\nClustering sequences with UMAP and HDBscan...')
    print('"Noisy" (not assigned to any cluster) reads will be removed')
    skip, checks = log_checker(log, SAMPLES, f'{OUT}/{{}}/clusters.tsv')
    if all(checks) and os.path.exists(f'{OUT}/shared_clusters.tsv'):
        return print(f'Clusters for all samples were already created. Skipping')
    os.environ["MKL_NUM_THREADS"] = "1" 
    os.environ["NUMEXPR_NUM_THREADS"] = "1" 
    os.environ["OMP_NUM_THREADS"] = "1" 
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
    bash(f'rm -rf {OUT}/../Final_output {OUT}/../Read_correction {OUT}/../Taxonomy_annotation/*.tsv \
           {OUT}/../Logs/Read_correction.log {OUT}/../Logs/Taxonomy_annotation.log')
    for sample in SAMPLES:
        if os.path.exists(f'{OUT}/{sample}/clusters.tsv') and sample in skip:
            continue
        #get clusters
        data = pd.read_csv(f'{OUT}/{sample}/kmers.tsv', sep='\t', index_col=0)
        size = max([CLUST_UQ, 10])
        clust_emb = umap.UMAP(n_jobs=T, metric='braycurtis', min_dist=0, n_components=10).fit_transform(data.values)
        labels = cluster.HDBSCAN(min_cluster_size=size, n_jobs=T, cluster_selection_epsilon=EPS,copy=False).fit_predict(clust_emb)
        clusters = pd.DataFrame({'Feature': data.index, 'Cluster': labels})
        clusters = clusters.loc[clusters.Cluster >= 0]
        clusters.Cluster = 'Cluster_' + clusters.Cluster.astype(str)
        for cid in clusters.Cluster.unique():
            sub = clusters.loc[clusters.Cluster == cid].copy()
            if len(sub) > 100:
                sub = sub.sample(n=100, random_state=RSTAT)
            data = data.copy()
            data.loc[sub.Feature.tolist(),'FullID'] = sample+'___'+cid+'___'
        data = data[data['FullID'].notna()]
        data.FullID = data.FullID + data.index.astype(str)
        data.set_index('FullID', inplace=True)
        data.to_csv(f'{OUT}/{sample}/subsampled_ids.tsv', sep='\t')
        clusters.to_csv(f'{OUT}/{sample}/clusters.tsv', sep='\t', index=False)
        bash(f'echo "{len(clusters)} features were clustered into {len(set(labels))} clusters" >> {log}')   
    #cluster clusters
    dfs = []
    for sample in SAMPLES:
        df = pd.read_csv(f'{OUT}/{sample}/subsampled_ids.tsv', sep='\t', index_col=0)
        dfs.append(df)
    dfs = [df.apply(pd.to_numeric, downcast='integer') for df in dfs]
    data = pd.concat(dfs)
    clust_emb = umap.UMAP(n_jobs=T, metric='braycurtis',  min_dist=0, n_components=10).fit_transform(data.values)
    labels = cluster.HDBSCAN(min_cluster_size=8, n_jobs=T, cluster_selection_epsilon=EPS,copy=False).fit_predict(clust_emb)
    clusters = pd.DataFrame({'Feature': data.index, 'Cluster': labels})
    clusters = clusters.loc[clusters.Cluster >= 0]
    clusters.Cluster = 'Cluster_' + clusters.Cluster.astype(str)
    clusters.to_csv(f'{OUT}/shared_clusters.tsv', sep='\t', index=False)
    bash(f'echo "Subsampled by cluster features were clustered between samples" >> {log}')
    print("Subsampled by cluster features were clustered between samples")

    
#Function to pool clusters from samples to shared clusters and recalculate abundances
def shared_clusters(OUT, FI, SAMPLES, RSTAT, SUBS, T, log):
    import umap
    print('\nPooling shared clusters and recalculating abundances...')
    if os.path.exists(f'{FI}/cluster_counts.tsv'):
        return print(f'Clusters already pooled and recalculated. Skipping')
    shared = pd.read_csv(f'{OUT}/shared_clusters.tsv', sep='\t', index_col=0)
    counts = pd.DataFrame(columns=SAMPLES, index=shared.Cluster.unique())
    counts = counts.astype(float).fillna(0)
    clust_dict = {c:[] for c in shared.Cluster.unique()}
    i = len(clust_dict)-1
    for sample in SAMPLES:
        unique = pd.read_csv(f'{OUT}/{sample}/clusters.tsv', sep='\t', index_col=0)
        subunique = pd.read_csv(f'{OUT}/{sample}/subsampled_ids.tsv', sep='\t', index_col=0)
        for uclust in unique.Cluster.unique():
            uniq = unique.loc[unique.Cluster == uclust]
            subsize = len(subunique.loc[subunique.index == uclust]) 
            shar = shared.loc[shared.index.str.contains(f"{sample}___{uclust}___")]
            if len(shar) > 0:
                shar = shar.groupby('Cluster').size().reset_index(name='counts')
                shar = shar.sort_values('counts', ascending=False).reset_index()
                if shar.loc[0, 'counts'] >= subsize/2:
                    clust_dict[shar.loc[0, 'Cluster']] += uniq.index.tolist()
                    counts.loc[shar.loc[0, 'Cluster'], sample] += len(uniq)
                    continue
            i += 1
            counts.loc[f'Cluster_{i}', sample] = len(uniq)
            clust_dict.update({f'Cluster_{i}': uniq.index.tolist()})     
    counts = counts.astype(float).fillna(0)
    counts = counts.loc[counts.sum(axis=1) >= 10]
    counts.index.name = 'Cluster'
    
    #write features by clusters
    print(f'Big clusters will be subsampled to {SUBS} reads for read correction!')
    random.seed(RSTAT)
    bash(f"mkdir -p {OUT}/Clusters_subsampled {FI}")
    with open(f"{OUT}/Clusters_subsampled/Pooled.txt", "w") as pooled:
        for k, v in clust_dict.items():
            pooled.write("{}\n{}\n".format(k, '\n'.join(v)))
            if len(v) > SUBS:
                v = random.sample(v, SUBS)
            with open(f"{OUT}/Clusters_subsampled/{k}.txt", "w") as clust:
                clust.write("{}".format('\n'.join(v)))
    counts.sort_index(key=lambda x: (x.to_series().str[8:].astype(int)), inplace=True)
    counts.to_csv(f'{FI}/cluster_counts.tsv', sep='\t')
    bash(f'echo "\nFeatures were stored by each shared cluster" >> {log}')
    print("\nFeatures were stored by each shared cluster")
    
    
###Functions to split sample files by clusters and finding consensus for each cluster

#working function
def file_splitter(out, cluster, log):
    fq = f'{out}/{cluster}.fq'
    fastq = f"{out}/../pooled.fq"
    bash(f"zgrep -f {out}/{cluster}.txt -F -A 3 {fastq} | grep -v '^--$' > {fq}")
    consensus = bash(f"spoa {fq}")
    with open(f'{out}/{cluster}_consensus.fa', 'wt') as cons:
        cons.write(">{}\n{}\n".format(cluster, str(consensus).split('\\n')[1]))
    bash(f'gzip -f {fq}')
    bash(f'echo "{cluster} done. Enjoy" >> {log}')
   
       
def file_by_cluster(INPUT, subs, OUT, T, log):
    out = f'{OUT}/Clusters_subsampled'
    df = pd.read_csv(f"{OUT}/../Final_output/cluster_counts.tsv", sep='\t', index_col=0)
    clusters = df.index.tolist()
    file = f'{OUT}/consensus_pooled.fa'
    skip, checks = log_checker(log, clusters, f'{out}/{{}}_consensus.fa')
    if all(checks) and os.path.exists(file):
        if os.stat(file).st_size > 0:
            return print(f'Fastq files and consensuses for all clusters exists. Skipping')
    #pool
    print('\nCreating files for subsampled clusters...')
    print(f'Big clusters will be subsampled to {subs} reads!')
    big = f"{OUT}/pooled.fq"
    if not os.path.exists(big):
        inp = glob.glob(f"{INPUT}/*.f*q*")
        if inp[0].endswith('.gz'):
            big = f"{OUT}/pooled.fq.gz"
        bash(f'cat {" ".join(inp)} > {big}')
        if big.endswith('.gz'):
            bash(f'gunzip {big}')
    #must use Manager queue here, or will not work
    pool = mp.Pool(T)
    #fire off workers
    jobs = []
    for cluster in clusters:
        if os.path.exists(f'{out}/{cluster}.fq') and os.path.exists(f'{out}/{cluster}_consensus.fa'):
            continue
        job = pool.apply_async(file_splitter, (out, cluster, log))
        jobs.append(job)
    [job.get() for job in jobs]
    #pool consensuses
    with open(file, 'w') as pooled:
        for cluster in clusters:
            with open(f'{out}/{cluster}_consensus.fa', 'r') as cons:
                pooled.write(cons.read())
                
