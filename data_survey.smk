configfile: 'data_survey.config.yaml'
workdir: config['workdir']

localrules: all, aggregate_results, plot_results


def len_cutoff(wildcards, trim_percent_cutoff=.8):
    """Compute minimal length based on average read length."""
    import numpy as np
    from Bio import SeqIO

    if wildcards.accession in config['samples_se']:
        fname = f'data/SE/{wildcards.accession}.fastq'
    elif wildcards.accession in config['samples_pe']:
        fname = f'data/PE/{wildcards.accession}_1.fastq'
    else:
        raise RuntimeError(f'Invalid accession: {wildcards.accession}')

    read_len = np.mean(
        [len(r.seq) for r in SeqIO.parse(fname, 'fastq')]
    ).astype(int)

    len_cutoff = int(trim_percent_cutoff * read_len)
    return len_cutoff


def gather_trimmed_input_files(wildcards):
    if wildcards.accession in config['samples_se']:
        return [
            f'trimmed/SE/{wildcards.accession}.fastq'
        ]
    elif wildcards.accession in config['samples_pe']:
        return [
            f'trimmed/PE/{wildcards.accession}_1.fastq',
            f'trimmed/PE/{wildcards.accession}_2.fastq'
        ]
    else:
        raise RuntimeError(f'Invalid accession: {wildcards.accession}')


rule all:
    input:
        'plots/'


rule get_fastq_se:
    output:
        'data/SE/{accession}.fastq'
    resources:
        mem_mb = 5_000
    threads: 6
    wrapper:
        '0.51.3/bio/sra-tools/fasterq-dump'


rule get_fastq_pe:
    output:
        'data/PE/{accession}_1.fastq',
        'data/PE/{accession}_2.fastq'
    resources:
        mem_mb = 5_000
    threads: 6
    wrapper:
        '0.51.2/bio/sra-tools/fasterq-dump'


rule vpipe_trim_se:
    input:
        fname = 'data/SE/{accession}.fastq'
    output:
        'trimmed/SE/{accession}.fastq',
    params:
        extra = '-ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10',
        len_cutoff = len_cutoff
    log:
        outfile = 'logs/prinseq.{accession}.out.log',
        errfile = 'logs/prinseq.{accession}.err.log'
    conda:
        'envs/preprocessing.yaml'
    shell:
        """
        echo "The length cutoff is: {params.len_cutoff}" > {log.outfile}

        prinseq-lite.pl \
            -fastq {input.fname} \
            {params.extra} \
            -out_format 3 \
            -out_good trimmed/SE/{wildcards.accession} \
            -out_bad null \
            -min_len {params.len_cutoff} \
            -log {log.outfile} 2> >(tee {log.errfile} >&2)
        """


rule vpipe_trim_pe:
    input:
        fname1 = 'data/PE/{accession}_1.fastq',
        fname2 = 'data/PE/{accession}_2.fastq'
    output:
        'trimmed/PE/{accession}_1.fastq',
        'trimmed/PE/{accession}_2.fastq'
    params:
        extra = '-ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10',
        len_cutoff = len_cutoff
    log:
        outfile = 'logs/prinseq.{accession}.out.log',
        errfile = 'logs/prinseq.{accession}.err.log'
    conda:
        'envs/preprocessing.yaml'
    shell:
        """
        echo "The length cutoff is: {params.len_cutoff}" > {log.outfile}

        prinseq-lite.pl \
            -fastq {input.fname1} \
            -fastq2 {input.fname2} \
            {params.extra} \
            -out_format 3 \
            -out_good trimmed/PE/{wildcards.accession} \
            -out_bad null \
            -min_len {params.len_cutoff} \
            -log {log.outfile} 2> >(tee {log.errfile} >&2)
        """


rule bwa_index:
    input:
        srcdir(config['reference'])
    output:
        'references/reference.amb',
        'references/reference.ann',
        'references/reference.bwt',
        'references/reference.pac',
        'references/reference.sa'
    log:
        'logs/bwa_index.log'
    params:
        prefix = 'references/reference'
    wrapper:
        '0.51.2/bio/bwa/index'


rule bwa_mem:
    input:
        reads = gather_trimmed_input_files,
        index = 'references/reference.amb'
    output:
        'alignment/{accession}.bam'
    log:
        'logs/bwa_mem_{accession}.log'
    params:
        index = 'references/reference',
        sort = 'samtools',
        sort_order = 'coordinate',
    resources:
        mem_mb = 16_000
    threads: 8
    wrapper:
        '0.51.2/bio/bwa/mem'


rule samtools_index:
    input:
        'alignment/{accession}.bam'
    output:
        'alignment/{accession}.bam.bai'
    wrapper:
        '0.51.2/bio/samtools/index'


rule compute_coverage:
    input:
        fname = 'alignment/{accession}.bam',
        index = 'alignment/{accession}.bam.bai'
    output:
        fname = 'coverage/coverage.{accession}.csv'
    run:
        import pysam

        import numpy as np
        import pandas as pd

        bam = pysam.AlignmentFile(input.fname, 'rb')
        assert bam.has_index()
        assert len(bam.references) == 1
        ref = bam.references[0]

        coverage = np.sum(bam.count_coverage(ref), axis=0)

        pd.DataFrame({
            'accession': [wildcards.accession],
            'coverage': [np.mean(coverage)]
        }).to_csv(output.fname, index=False)


rule aggregate_results:
    input:
        fname_list = expand(
            'coverage/coverage.{accession}.csv',
            accession=config['samples_se'] + config['samples_pe'])
    output:
        fname = 'results/results.csv'
    run:
        import pandas as pd

        df_list = []
        for fname in input.fname_list:
            tmp = pd.read_csv(fname)
            df_list.append(tmp)
        df = pd.concat(df_list)

        df.to_csv(output.fname, index=False)


rule plot_results:
    input:
        fname = 'results/results.csv'
    output:
        dname = directory('plots/')
    run:
        from pathlib import Path

        import pandas as pd

        import seaborn as sns
        import matplotlib.pyplot as plt

        # fix big boom on macOS
        import matplotlib
        matplotlib.use('Agg')

        # read data
        dname = Path(output.dname)
        df = pd.read_csv(input.fname)

        # plot data
        plt.figure(figsize=(8, 6))

        sns.distplot(df['coverage'], kde=False)

        plt.xlabel('Coverage per Accession')
        plt.ylabel('Count')

        plt.tight_layout()
        plt.savefig(dname / 'coverage_histogram.pdf')
