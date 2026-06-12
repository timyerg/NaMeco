import os
import glob
from ._helper import *



#Function to run Chopper
def chopper(INPUT, SAMPLES, T, Q, MINL, MAXL, OUT, LOGS, log):
    print('Running chopper...')
    skip, checks = log_checker(log, SAMPLES, f'{OUT}/{{}}.fq.gz')
    if all(checks):
        return print(f'All samples were already chopped. Skipping')
    bash(f'mkdir -p {OUT} {LOGS}')
    for sample in SAMPLES:
        fq_out = f'{OUT}/{sample}.fq.gz'
        if os.path.exists(fq_out) and sample in skip:
            continue
        file = glob.glob(f"{INPUT}/{sample}.f*q*")[0]
        bash(f'echo "\n##### Processing {sample} #####\n" >> {log}')
        bash(f'echo "Chopping {sample}" >> {log}')
        pre = f'gunzip -c {file} -q |' if file.endswith('gz') else f'cat {file} |'
        bash(f'{pre} chopper -q {Q} -l {MINL} --maxlength {MAXL} -t {T} 2>> {log} | gzip > {fq_out}')
        bash(f'echo "{sample} done. Enjoy" >> {log}')
        

#Function to find samples with too low counts
def min_sample_size(INPUT, SAMPLES, OUT, LOGS, MinSS, log):
    print(f'\nSamples with sample size after QC less than {MinSS} reads will be ignored')
    skip, checks = log_checker(log, SAMPLES, f'{OUT}/{{}}.fq.gz', 'has enough reads\n')
    if all(checks):
        return print(f'All samples were already checked for sequencing depth. Skipping')
    bash(f'echo "\n##### Checking sample sizes after QC #####\n" >> {log}')
    LOWS = []
    bash(f'mkdir -p {LOGS}')
    for sample in SAMPLES:
        fq = glob.glob(f"{INPUT}/{sample}.f*q*")[0]
        ss = int(bash(f'zcat {fq}|wc -l'))/4 if fq.endswith('gz') else int(bash(f'cat {fq}|wc -l'))/4
        if ss < MinSS:
            bash(f'echo "{sample} with {int(ss)} reads will be ignored" >> {log}')
            LOWS.append(sample)
            print(f'{sample} with {int(ss)} reads will be ignored. If needed, adjust "--min_sample_size"')
        bash(f'echo "{sample} has enough reads" >> {log}')
    return LOWS
