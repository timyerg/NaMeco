import os
import pandas as pd
from ._helper import *


###Functions for taxonomy annotation
#def fetch_db(fetch, DBpath, log):
    



#Percent identity thresholds to mask false positive annotations
def taxonomy_thresholds(bclust, thresholds):
    for ind in bclust.index:
        taxon = bclust.loc[ind, 'Taxon']
        last = ''
        for rank, perc in thresholds.items():
            prefix = f"{rank[0].lower()}__"
            pat = taxon.split(prefix)[-1].split(';')[0]
            if bclust.loc[ind, 'pind'] >= perc:
                last = pat
            if bclust.loc[ind, 'pind'] < perc:
                taxon = taxon.replace(prefix+pat, f"{prefix}{last} unclassified")
                bclust.loc[ind, 'Taxon'] = taxon
    return(bclust)

#select top-hit based on consensus taxonomy
def top_hit(bclust, taxa, frac):
    taxa_counts = bclust["Taxon"].value_counts()
    bclust["Taxa_counts"] = bclust["Taxon"].map(taxa_counts)
    bclust.sort_values(["Taxa_counts", 'bitscore', 'pind'], ascending=[False, False, False], inplace=True)
    if len(bclust) == 0:
        taxon, pind = ';'.join([f'{r}__Unclassified' for r in 'dpcofgs']), 0
    else:
        taxon, pind = bclust.Taxon.iloc[0], bclust.pind.iloc[0]
        if len(bclust.loc[bclust.Taxa_counts==bclust.Taxa_counts.max()])/len(bclust) < frac:
            taxon = taxon.rsplit(';',1)[0] +';'+ taxon.rsplit(';',1)[-1].split(' ')[0] + ' unclassified'
    return taxon, pind

def taxonomy_annotation(DB, DBV, MASK, gap, frac, T, OUT, FI, DBpath, log):
    print(f'Starting taxonomy annotations with blastn against {DB}...')
    Q=f'{FI}/rep_seqs.fasta'
    DBpath = DBpath.format(OUT=OUT, DBV=DBV, DB=DB)
    queries = [l[1:].split(' ')[0].split('\n')[0] for l in open(Q, 'rt') if l.startswith('>')]
    thresholds = {'Domain': 65, 'Phylum': 75, 'Class': 78.5,
                  'Order': 82, 'Family': 86.5, 'Genus': 94.5, 'Species': 97}
    taxa = pd.DataFrame(columns=['Taxon', 'Perc. id.'])
    bash(f'mkdir -p {OUT} {FI}')
    
    #create DB
    if not os.path.exists(f"{DBpath}/ssu_all.fna.ndb"):
        print(f'Creating database...')
        substr = DBV if DBV == 'latest' else f"release{DBV.split('.')[0]}/{DBV}"
        suffix = "ssu_all" if DBV == 'latest' else f"ssu_all_r{DBV.split('.')[0]}"
        seq = f"https://data.ace.uq.edu.au/public/gtdb/data/releases/{substr}/genomic_files_all/{suffix}.fna.gz"
        bash(f'mkdir -p {DBpath}')
        bash(f'wget -P {DBpath} {seq} 2>> {log}')
        if not os.path.exists(f"{DBpath}/db.fa.gz"):
            bash(f'mv {DBpath}/ssu_all*.fna.gz {DBpath}/db.fa.gz')
        bash(f'gunzip {DBpath}/db.fa.gz 2>> {log}')
        bash(f'makeblastdb -in {DBpath}/db.fa -parse_seqids -dbtype "nucl"')
        
        #mapping file
        with open(f'{DBpath}/db.fa', 'rt') as fa:
            ls = [l[1:].replace(' d_','\td_').split(' [')[0] for l in fa.readlines() if l.startswith('>')]
            with open(f'{DBpath}/map.tsv', 'wt') as ref:
                ref.write('SeqID\tTaxonomy\n')
                ref.write('\n'.join(ls))
        bash(f'rm {DBpath}/db.fa')
        
    else:
        print('Database exists. Skipping')
        bash(f'echo "{DB} database exists. Skipping." >> {log}')
        
    #annotate
    if not os.path.exists(f"{OUT}/blastn.tsv"):
        print(f'\nAssigning taxonomy...')
        bash(f'blastn -query {Q} -db {DBpath}/db.fa -task blastn -qcov_hsp_perc 80 \
               -num_threads {T} -out {OUT}/blastn.tsv -max_target_seqs 50 -max_hsps 50 \
               -outfmt "6 qseqid sseqid evalue length pident nident bitscore score gaps" 2>> {log}')
    else:
        print('\nBlastn output exists. Skipping')
        bash(f'echo "Blastn output exists. Skipping." >> {log}')

    #select tophit taxonomy
    blast = pd.read_csv(f"{OUT}/blastn.tsv", sep='\t', header=None, 
            names=['Cluster', 'SeqID', 'eval', 'length', 'pind', 'nind', 'bitscore', 'score', 'gaps'])
    blast = blast.sort_values(['bitscore', 'eval'], ascending=[False, False])

    #get full taxonomies
    if not os.path.exists(f"{FI}/Taxonomy.tsv"):
        print('\nMapping GTDB to get full taxonomies...')
        mapp = pd.read_csv(f'{DBpath}/map.tsv', sep='\t')
        mapp.Taxonomy = mapp.Taxonomy.apply(lambda x: x.rsplit(';', 1)[0] +';'+ 
                     ' '.join(x.rsplit(';', 1)[-1].split(' ')[:2]).replace('_', ' ').replace('  ', '__'))
        mapping = dict(mapp[['SeqID', 'Taxonomy']].values)
        for cluster in queries:
            bclust = blast.loc[blast.Cluster == cluster].copy()
            #apply "Gap" filtering
            bclust = bclust.loc[bclust.bitscore > bclust.bitscore.max() - gap]
            #add taxonomies with proper percent identity thresholds
            bclust['Taxon'] = bclust['SeqID'].map(mapping)
            if MASK:
                bclust = taxonomy_thresholds(bclust, thresholds)
            #select top hit based on frequency
            taxa.loc[cluster, ['Taxon', 'Perc. id.']] = top_hit(bclust, taxa, frac)  
    else:
        print('\nTaxonomy exists. Skipping')
        bash(f'echo "Taxonomy exists. Skipping." >> {log}')  
    if len(taxa) != 0:
        for cluster in queries:
            if cluster not in blast.Cluster.tolist():
                taxa.loc[cluster, 'Taxon'] = 'Unclassified'
        taxa.index.rename('Feature ID', inplace=True)
        taxa.to_csv(f'{FI}/Taxonomy-q2.tsv', sep='\t')
        for rank in thresholds:
            taxa[rank] = taxa.Taxon.apply(lambda x: x.split(f"{rank[0].lower()}__")[-1].split(';')[0])
        taxa.drop('Taxon', axis=1, inplace=True)
        taxa.index.rename('Cluster', inplace=True)
        taxa.to_csv(f'{FI}/Taxonomy.tsv', sep='\t')
        
    #collapse taxonomies
    print('\nChecking if collapsed taxonomies exist...')
    taxa = pd.read_csv(f'{FI}/Taxonomy.tsv', sep='\t', index_col=0)
    counts = pd.read_csv(f'{FI}/cluster_counts.tsv', sep='\t', index_col=0)
    for rank in thresholds:
        if os.path.exists(f"{rank}_counts.tsv"):
            continue
        print(f'Collapsing to {rank}')
        coll = counts.copy()
        coll[rank] = taxa[rank]
        coll = coll.groupby(rank).sum()
        coll.to_csv(f"{FI}/{rank}_counts.tsv", sep='\t')
        
