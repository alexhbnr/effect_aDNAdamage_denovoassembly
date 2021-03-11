####################################################################################################
# Simulation of short-read sequenicng data with ancient DNA damage for three microbial taxa,
# assembly with MEGAHIT, and comparison of the assembled contigs to the respective reference genome
#
# Alex Huebner, 15/02/21
####################################################################################################

from glob import glob
import json
import os
import re
from tqdm import tqdm

from scipy.stats import lognorm
import numpy as np
import pandas as pd
import pysam
import pyfastx

#### Microbial species #############################################################################
URLS = {'Msmithii': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/525/GCF_000016525.1_ASM1652v1/GCF_000016525.1_ASM1652v1_genomic.fna.gz',
        'Tforsythia': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/215/GCF_000238215.1_ASM23821v1/GCF_000238215.1_ASM23821v1_genomic.fna.gz',
        'Adentalis': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/429/225/GCF_000429225.1_ASM42922v1/GCF_000429225.1_ASM42922v1_genomic.fna.gz',
       }
####################################################################################################

wildcard_constraints:
    genome = "[A-Za-z0-9_]+",
    rl = "[a-z]+",
    damage = "[0-9]",
    covbin = "[0-9]",
    i = "[1-2]"

rule all:
    input:
        "reference_genomes.done",
        "evaluation.done"

#### Prepare reference genomes #####################################################################

rule reference_genomes:
    input:
        expand("genomes/{genome}.fna.gz.fai", genome=URLS.keys()),
        "results/gccontent.txt",
        "results/genomelength.txt",
        "results/chromlength.txt"
    output:
        touch("reference_genomes.done")

rule download_genomes:
    output:
        "genomes/{genome}.fna.gz"
    message: "Download genome of species {wildcards.genome}"
    params: 
        url = lambda wildcards: URLS[wildcards.genome]
    shell:
        """
        curl {params.url} | zcat - | bgzip > {output}
        """

rule index_fastas:
    input:
        "genomes/{genome}.fna.gz"
    output:
        "genomes/{genome}.fna.gz.fai"
    message: "Index downloaded FastA file of {wildcards.genome}"
    conda: "simulation.yaml"
    shell:
        "samtools faidx {input}"

rule determine_gc_content:
    input:
        expand('genomes/{genome}.fna.gz', genome=URLS.keys())
    output:
        'results/gccontent.txt'
    message: "Calculate GC content of the downloaded genomes"
    conda: "simulation.yaml"
    shell:
        """
        echo -e 'sample\tGC' > {output}
        for g in {input}; do
            gc=$(bioawk -c fastx '{{
                    gcs = gcs + gc($seq)
                 }} END {{
                    print gcs / NR
                 }}' ${{g}})
            echo -e "$(basename ${{g}} .fna.gz)\t${{gc}}"
        done >> {output}
        """

checkpoint determine_genome_length:
    input:
        expand('genomes/{genome}.fna.gz', genome=URLS.keys())
    output:
        'results/genomelength.txt'
    message: "Calculate genome lengths of reference genomes"
    params:
        dir = 'genomes'
    run:
        with open(output[0], "wt") as outfile:
            outfile.write("sample\tgenomelength\n")
            for ref in list(URLS.keys()):
                length = 0
                for name, seq in pyfastx.Fasta(f"{params.dir}/{ref}.fna.gz", build_index=False):
                    length += len(seq)
                outfile.write(f"{ref}\t{length}\n")

rule determine_chrom_length:
    input:
        expand('genomes/{genome}.fna.gz', genome=URLS.keys())
    output:
        'results/chromlength.txt'
    message: "Calculate the individual chromosome lengths of reference genomes"
    params:
        dir = 'genomes'
    run:
        chrom_lengths = pd.concat([pd.DataFrame([(name, len(seq))
                                                 for name, seq in pyfastx.Fasta(f"{params.dir}/{genome}.fna.gz",
                                                 build_index=False)], columns=['chr', 'length'])
                                     .assign(genome=genome)
                                   for genome in list(URLS.keys())])
        chrom_lengths[['genome', 'chr', 'length']] \
            .to_csv(output[0], sep="\t")

####################################################################################################

#### Simulate short-read sequencing data using gargammel ###########################################

# Parameters for simulations with gargammel
## read length profiles
READLENGTHS = {'short': (0.25, 75),  # sigma, scale for lognorm.rvs
               'medium': (0.35, 100),
               'long': (0.45, 150)}
## number of reads per read length profile 
NREADS = {'short': 30000000,
          'medium': 21000000,
          'long': 15000000}
## damage levels
DAMAGE = {level: amount
          for level, amount in zip(list(range(1, 10)), [0.01, 0.015, 0.02, 0.025, 0.03, 0.05, 0.1, 0.15, 0.2])}
# Coverage bins
COVBINS = [(0, 1), (1, 2), (2, 5), (5, 10), (10, 20), (20, 50), (50, 100), (100, 100)]
# Briggs et al. aDNA model parameters for each profile
BRIGGSPARAMS = pd.DataFrame.from_dict({'readlength': ['short', 'medium', 'long'],
                                       'l': [0.5, 0.5, 0.2],
                                       'intercept': [0.00285351934051996, 0.006, 0.000806560431244065],
                                       's': [0.510462904248573, 0.46, 0.793827313143463]}).set_index('readlength')

#### Auxilliary functions for simulating short-read data #######################

def infer_s(rl, ctfreq):
    """ Infer the parameter s for the amount of C-to-T substitutions for the
        specified readlength profile.

        Parameters:
        rl     : read length profile name 
        ctfreq : final frequency of C-to-T substitutions at the 5' end
    """
    p = BRIGGSPARAMS.loc[rl]
    return round((ctfreq - p.intercept) / p.s, 3)


def subsampling_fraction(wildcards, covbin, seed=0):
    """ Determines the fraction required for subsampling to a specific coverage.
    
        Parameters:
        covbin       : tuple with interval for final coverage along reference genome
        seed         : seed used for random state initialisation
    """
    np.random.seed(seed) 
    coverage = np.random.uniform(covbin[0], covbin[1], size=1)[0]
    rlp = pd.read_csv(checkpoints.generate_readlength_profile.get(**wildcards).output[0],
                      sep = "\t", header=None, names=['length', 'frac'])
    rlp = pd.read_csv("/mnt/scratch/alexh/wibowo_assembly/results/rldistribution_short_sizefreq.size.gz", sep="\t",
                      header=None, names=['length', 'frac'])
    totalbases = (rlp['length'] * rlp['frac']).sum() * NREADS[wildcards.rl]
    genomelength = pd.read_csv(checkpoints.determine_genome_length.get(**wildcards).output[0], sep="\t",
                               index_col=['sample']).at[wildcards.genome, 'genomelength']
    estbases = coverage * genomelength
    return round(seed + (estbases / totalbases), 6)

rule simulate:
    input:
        expand("fastqs/{genome}-{rl}.{damage}.{covbin}_1.trimmed_fastp.fastq.gz", genome=URLS.keys(), rl=READLENGTHS.keys(), damage=DAMAGE.keys(), covbin=range(len(COVBINS)))
    output:
        touch("simulate.done")

checkpoint generate_readlength_profile:
    output:
        "results/rldistribution_{rl}_sizefreq.size.gz"
    message: "Generate fragment size distribution for read length profile {wildcards.rl}"
    run:
        sigma, scale = READLENGTHS[wildcards.rl]
        lengths_df = pd.DataFrame.from_dict({'length': lognorm.rvs(sigma, loc=-10,
                                                                   scale=scale,
                                                                   size=100000,
                                                                   random_state=1)})
        lengths_df['length'] = np.floor(lengths_df['length']).astype(int)
        lengths_df = lengths_df.query('length >= 35 & length < 300')
        frequency = lengths_df.groupby(['length'])[['length']].count() \
            .rename({'length': 'n'}, axis=1)
        frequency['freq'] = frequency['n'] / frequency['n'].sum()
        frequency[['freq']].to_csv(output[0], sep="\t", header=False,
                                   compression="gzip", float_format="%.6f")

rule fragSim:
    input:
        rl = expand("results/rldistribution_{rl}_sizefreq.size.gz", rl=READLENGTHS.keys()),
        refgenome = "reference_genomes.done"
    output:
        temp("gargammel/{genome}-{rl}.bam")
    message: "Simulate DNA molecules with read length distribution {wildcards.rl} for genome {wildcards.genome}"
    conda: "simulation.yaml"
    params:
        rl_profile = "results/rldistribution_{rl}_sizefreq.size.gz",
        nreads = lambda wildcards: NREADS[wildcards.rl],
        genome = "genomes/{genome}.fna.gz"
    shell:
        """
        fragSim \
            -n {params.nreads} \
            -b {output} \
            -u \
            -m 35 \
            -M 300 \
            -f {params.rl_profile} \
            {params.genome} 
        """

rule deamSim:
    input:
        "gargammel/{genome}-{rl}.bam"
    output:
        "gargammel/{genome}-{rl}.{damage}.bam"
    message: "Simulate reads of read length profile {wildcards.rl} with damage level {wildcards.damage} for genome {wildcards.genome}"
    conda: "simulation.yaml"
    params:
        s = lambda wildcards: infer_s(wildcards.rl, DAMAGE[int(wildcards.damage)]),
        l = lambda wildcards: BRIGGSPARAMS.loc[wildcards.rl, 'l']
    shell:
        """
        deamSim \
            -b {output} \
            -damage 0.03,{params.l},0.01,{params.s} \
            {input}
        """

rule adptSim:
    input:
        "gargammel/{genome}-{rl}.{damage}.bam"
    output:
        temp("gargammel/{genome}-{rl}.{damage}.adapters.bam")
    message: "Add Illumina adapters to genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage}"
    conda: "simulation.yaml"
    params:
        foradpt = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
        revadpt = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGCGCCGTATCATT",
        length = 150
    shell:
        """
        /home/alexander_huebner/opt/gargammel/src/adptSim \
            -bp {output} \
            -f {params.foradpt} \
            -s {params.revadpt} \
            -l {params.length} \
            {input}
        """

rule subsample_bam2fq:
    input:
        "gargammel/{genome}-{rl}.{damage}.adapters.bam"
    output:
        pe1 = temp("fastqs/{genome}-{rl}.{damage}.{covbin}_1.nobq.fastq.gz"),
        pe2 = temp("fastqs/{genome}-{rl}.{damage}.{covbin}_2.nobq.fastq.gz")
    message: "Subsample the sequencing data of genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin}"
    conda: "simulation.yaml"
    params: 
        # frac = lambda wildcards: subsampling_fraction(COVBINS[int(wildcards.covbin)], wildcards.rl, NREADS[wildcards.rl], wildcards.genome, seed=1000)
        frac = lambda wildcards: subsampling_fraction(wildcards, COVBINS[int(wildcards.covbin)], seed=1000)
    shell:
        """
        samtools view -bh -s {params.frac} {input} | samtools fastq -1 {output.pe1} -2 {output.pe2} -
        """

rule fix_bq:
    input:
        "fastqs/{genome}-{rl}.{damage}.{covbin}_{i}.nobq.fastq.gz"
    output:
        temp("fastqs/{genome}-{rl}.{damage}.{covbin}_{i}.fastq.gz")
    message: "Set base qualities of genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin} to 37"
    conda: "simulation.yaml"
    shell:
        """
        bioawk -c fastx '{{
            gsub(/!/, "F", $qual);
            print "@" $name " " $comment;
            print $seq;
            print "+";
            print $qual
        }}' {input} | gzip > {output}
        """

rule adapterremoval:
    input:
        lambda wildcards: [f"fastqs/{wildcards.genome}-{wildcards.rl}.{wildcards.damage}.{wildcards.covbin}_{i}.fastq.gz" for i in range(1, 3)]
    output:
        pe1 = temp("fastqs/{genome}-{rl}.{damage}.{covbin}_1.trimmed.fastq.gz"),
        pe2 = temp("fastqs/{genome}-{rl}.{damage}.{covbin}_2.trimmed.fastq.gz")
    message: "Remove adapters from genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin}"
    conda: "simulation.yaml"
    params:
        prefix = "tmp/{genome}-{rl}.{damage}.{covbin}"
    log: "fastqs/{genome}-{rl}.{damage}.{covbin}_adapterremoval.log"
    threads: 8
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        AdapterRemoval --basename {params.prefix} \
                           --file1 {input[0]} \
                           --file2 {input[1]} \
                           --trimns \
                           --trimqualities \
                           --minquality 20 \
                           --minlength 30 \
                           --output1 {output.pe1} \
                           --output2 {output.pe2} \
                           --threads {threads} \
                           --qualitybase 33 \
                           --settings {log}
        rm {params.prefix}*
        """

rule fastp:
    input:
        pe1 = "fastqs/{genome}-{rl}.{damage}.{covbin}_1.trimmed.fastq.gz",
        pe2 = "fastqs/{genome}-{rl}.{damage}.{covbin}_2.trimmed.fastq.gz"
    output:
        pe1 = "fastqs/{genome}-{rl}.{damage}.{covbin}_1.trimmed_fastp.fastq.gz",
        pe2 = "fastqs/{genome}-{rl}.{damage}.{covbin}_2.trimmed_fastp.fastq.gz"
    message: "Clean reads for of genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin} using fastp"
    conda: "simulation.yaml"
    log: "fastqs/{genome}-{rl}.{damage}.{covbin}_fastp.json"
    threads: 8
    shell:
        """
        fastp --in1 {input.pe1} \
              --in2 {input.pe2} \
              --out1 {output.pe1} \
              --out2 {output.pe2} \
              -A -g -Q -L \
              -w {threads} --json {log} --html /dev/null
        """

#################################################################################

#### De-novo assemble genomes using MEGAHIT ####################################

rule assembly:
    input:
        expand("megahit/{genome}-{rl}.{damage}.{covbin}.megahit.fa", genome=URLS.keys(), rl=READLENGTHS.keys(), damage=DAMAGE.keys(), covbin=range(len(COVBINS)))
    output:
        touch("assembly.done")

rule megahit:
    input:
        pe1 = "fastqs/{genome}-{rl}.{damage}.{covbin}_1.trimmed_fastp.fastq.gz",
        pe2 = "fastqs/{genome}-{rl}.{damage}.{covbin}_2.trimmed_fastp.fastq.gz"
    output:
        "megahit/{genome}-{rl}.{damage}.{covbin}.megahit.fa"
    message: "Assemble reads of genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin} using MEGAHIT"
    conda: "simulation.yaml"
    params:
        memory = 256000000000,
        prefix = "tmp/{genome}-{rl}.{damage}.{covbin}",
        tmpdir = "tmp"
    log: "megahit/{genome}-{rl}.{damage}.{covbin}.megahit.log"
    resources:
        megahit = 1
    threads: 8
    shell:
        """
        megahit -1 {input.pe1} \
                -2 {input.pe2} \
                -t {threads} \
                -m {params.memory} \
                --tmp-dir {params.tmpdir} \
                --min-contig-len 1000 \
                --out-dir {params.prefix}
        mv {params.prefix}/final.contigs.fa {output}
        mv {params.prefix}/log {log}
        rm -r {params.prefix}
        """

################################################################################

#### Evaluation of assemblies ##################################################

rule evaluation:
    input:
        alnstats = "results/alignmentstats_blastn.csv",
        mismatches = "results/mismatches_blastn.csv",
        mismatchfrac = "results/mismatchfrac_blastn.csv"
    output:
        touch("evaluation.done")

rule makeblastdb:
    input:
        "genomes/{genome}.fna.gz"
    output:
        "blastdb/{genome}.fa.nhr"
    message: "Make blastn database from genome {wildcards.genome}"
    conda: "simulation.yaml"
    params:
        prefix = "blastdb/{genome}.fa"
    shell:
        """
        zcat {input} > {params.prefix}
        makeblastdb -in {params.prefix} -parse_seqids -dbtype nucl
        """

rule filter_contigs:
    input:
        "megahit/{genome}-{rl}.{damage}.{covbin}.megahit.fa"
    output:
        "evaluation/{genome}-{rl}.{damage}.{covbin}.megahit_filtered.fa"
    message: "Filter contigs for minimum length of 2.5 kb for genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin}"
    conda: "simulation.yaml"
    shell:
        """
        bioawk -c fastx '{{if (length($seq) >= 2500){{print ">" $name; print $seq}}}}' {input} > {output}
        """

rule blastn:
    input:
        reffa = "blastdb/{genome}.fa.nhr",
        fa = "evaluation/{genome}-{rl}.{damage}.{covbin}.megahit_filtered.fa"
    output:
        tab = "evaluation/{genome}-{rl}.{damage}.{covbin}.tab.gz",
        sam = "evaluation/{genome}-{rl}.{damage}.{covbin}.sam.gz"
    message: "Run BlastN on experiment for genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin}"
    conda: "simulation.yaml"
    params:
        reffa = "blastdb/{genome}.fa"
    shell:
        """
        blastn -query {input.fa} -db {params.reffa} -outfmt 6 | gzip > {output.tab}
        blastn -query {input.fa} -db {params.reffa} -outfmt 17 | gzip > {output.sam}
        """

rule evaluate_blastn:
    input:
        tab = "evaluation/{genome}-{rl}.{damage}.{covbin}.tab.gz",
        sam = "evaluation/{genome}-{rl}.{damage}.{covbin}.sam.gz",
        fa = "evaluation/{genome}-{rl}.{damage}.{covbin}.megahit_filtered.fa",
        reffa = "blastdb/{genome}.fa.nhr"
    output:
        "evaluation/{genome}-{rl}.{damage}.{covbin}.blastn.json"
    message: "Evaluate BLASTn results for genome {wildcards.genome} with read length profile {wildcards.rl} and damage level {wildcards.damage} to coverage bin {wildcards.covbin}"
    params:
        reffa = "blastdb/{genome}.fa"
    run:
        tab = pd.read_csv(input.tab, sep="\t", header=None,
                          names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                                 'gapopen', 'qstart', 'qend', 'sstart', 'send',
                                 'evalue', 'bitscore'])
        if tab.shape[0] > 0:
            tab['rank'] = tab.groupby(['qseqid'])['evalue'].rank(method="first", axis=1)
            tab['rank'] = tab['rank'].astype(int)

            samfile = pysam.AlignmentFile(input.sam, "rb")
            contigs = {name: seq
                       for name, seq in pyfastx.Fasta(input.fa, build_index=False)}
            refseq = {re.sub(r'\.[0-9]+$', '', name): seq
                      for name, seq in pyfastx.Fasta(input.reffa, build_index=False)}

        stats = {'noContigs': 0,
                 'totallength': 0,
                 'alignmentlength': 0,
                 'noDel': 0,
                 'noIns': 0,
                 'mismatches': {'AC': 0,
                                'AG': 0,
                                'AT': 0,
                                'CA': 0,
                                'CG': 0,
                                'CT': 0,
                                'GA': 0,
                                'GC': 0,
                                'GT': 0,
                                'TA': 0,
                                'TC': 0,
                                'TG': 0},
                 'mismatchesPerContig': [],
                 'pident': [],
                 'contiglength': []}
        comp = {'A': 'T',
                'C': 'G',
                'G': 'C',
                'T': 'A',
                '-': '-'}

        tab['sseqid'] = tab['sseqid'].str.replace(r'\.[0-9]+$', '', regex=True)
        if tab.shape[0] > 0:
            for r, aln in zip(tab.iterrows(), samfile):
                _, row = r
                if row['rank'] == 1:
                    stats['noContigs'] += 1
                    stats['totallength'] += len(contigs[row['qseqid']])
                    stats['alignmentlength'] += abs(row['qend'] - row['qstart'] + 1)
                    stats['mismatchesPerContig'].append(row['mismatch'])
                    stats['pident'].append(row['pident'])
                    stats['contiglength'].append(len(contigs[row['qseqid']]))
                    if row['mismatch'] > 0:
                        for op, oplength in aln.cigartuples:
                            if op == 1:
                                stats['noIns'] += 1
                            elif op == 2:
                                stats['noDel'] += 1

                        if row['sstart'] < row['send']:
                            refaln = list(refseq[row['sseqid']][(row['sstart'] - 1):row['send']])
                        else:
                            refaln = list([comp[a] for a in refseq[row['sseqid']][(row['send'] - 1):row['sstart']][::-1]])
                        qaln = list(contigs[row['qseqid']][(row['qstart'] - 1):row['qend']])
                        optypes = {0: "S", 1: "D", 2: "I"}
                        substtype = [optypes[op]
                                    for op, oplength in aln.cigartuples
                                    for i in range(oplength)
                                    if op in [0, 1, 2]]
                        for i, s in enumerate(substtype):
                            if s == "I":
                                refaln.insert(i, "-")
                            elif s == "D":
                                qaln.insert(i, "-")

                        for rallele, qallele in zip(refaln, qaln):
                            if rallele != qallele:
                                if row['sstart'] < row['send']:
                                    if f'{rallele}{qallele}' in stats['mismatches']:
                                        stats['mismatches'][f'{rallele}{qallele}'] += 1
                                else:
                                    if f'{comp[rallele]}{comp[qallele]}' in stats['mismatches']:
                                        stats['mismatches'][f'{comp[rallele]}{comp[qallele]}'] += 1

        with open(output[0], "wt") as outfile:
            json.dump(stats, outfile)

rule summary_blastn:
    input:
        expand("evaluation/{genome}-{rl}.{damage}.{covbin}.blastn.json", genome=URLS.keys(), rl=READLENGTHS.keys(), damage=DAMAGE.keys(), covbin=range(len(COVBINS)))
    output:
        alnstats = "results/alignmentstats_blastn.csv",
        mismatches = "results/mismatches_blastn.csv",
        mismatchfrac = "results/mismatchfrac_blastn.csv"
    message: "Summarise alignment and mismatch stats into a tabular format"
    params:
        dir = "evaluation"
    run:
        alnstats = []
        mismatches = []
        mismatchrate = []

        extract_fn = re.compile(r'([A-Z][a-z]+)-([a-z]+)\.([0-9])-0.([0-9]).blastn.json')

        for fn in glob(f"{params.dir}/*.blastn.json"):
            stat_dict = json.load(open(fn, "rt"))
            genome, rl, damagelevel, covbin = extract_fn.search(os.path.basename(fn)).groups()
            # Alignment stats
            alnstats.append((genome, rl, int(damagelevel), int(covbin),
                             stat_dict['noContigs'], stat_dict['totallength'],
                             stat_dict['alignmentlength']))
            # Mismatch stats
            stat_dict['mismatches']['Del'] = stat_dict['noDel']
            stat_dict['mismatches']['Ins'] = stat_dict['noIns']
            mismatches.append(pd.DataFrame.from_dict(stat_dict['mismatches'], orient="index").transpose()
                .assign(genome=genome)
                .assign(rl=rl)
                .assign(damagelevel=int(damagelevel))
                .assign(covbin=int(covbin)))
            # Mismatch rates
            mismatchrate.append(pd.DataFrame({'rate': stat_dict['mismatchesPerContig'],
                                              'pident': stat_dict['pident'],
                                              'contiglength': stat_dict['contiglength'],
                                              'contigID': list(range(len(stat_dict['mismatchesPerContig'])))})
                .assign(genome=genome)
                .assign(rl=rl)
                .assign(damagelevel=int(damagelevel))
                .assign(covbin=int(covbin)))

        # Write to file
        pd.DataFrame(alnstats, columns=['genome', 'rl', 'damagelevel', 'covbin',
                                        'noContigs', 'totallength', 'alignmentlength']) \
            .to_csv(output.alnstats, sep="\t", index=False)
        pd.concat(mismatches)[['genome', 'rl', 'damagelevel', 'covbin'] + list(stat_dict['mismatches'].keys())] \
            .to_csv(output.mismatches, sep="\t", index=False)
        pd.concat(mismatchrate)[['genome', 'rl', 'damagelevel', 'covbin', 'contigID', 'contiglength', 'rate', 'pident']] \
            .to_csv(output.mismatchfrac, sep="\t", index=False)
