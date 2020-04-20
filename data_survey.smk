configfile: 'data_survey.config.yaml'
workdir: config['workdir']


def len_cutoff(wildcards, trim_percent_cutoff=.8):
    """Compute minimal length based on average read length."""
    import numpy as np
    from Bio import SeqIO

    fname = f'data/{wildcards.accession}_1.fastq'
    read_len = np.mean(
        [len(r.seq) for r in SeqIO.parse(fname, 'fastq')]
    ).astype(int)

    len_cutoff = int(trim_percent_cutoff * read_len)
    return len_cutoff


rule all:
    input:
        'plots/'


rule get_fastq_pe:
    output:
        'data/{accession}_1.fastq',
        'data/{accession}_2.fastq'
    wrapper:
        '0.51.2/bio/sra-tools/fasterq-dump'


rule vpipe_trim:
    input:
        fname1 = 'data/{accession}_1.fastq',
        fname2 = 'data/{accession}_2.fastq'
    output:
        'trimmed/{accession}_1.fastq',
        'trimmed/{accession}_2.fastq'
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
            -out_good trimmed/{wildcards.accession} \
            -out_bad null \
            -min_len {params.len_cutoff} \
            -log {log.outfile} 2> >(tee {log.errfile} >&2)
        """


# rule kallisto_index:
#     input:
#         fasta = srcdir(config['reference'])
#     output:
#         index = 'references/reference.idx'
#     log:
#         'logs/kallisto_index.log'
#     wrapper:
#         '0.51.2/bio/kallisto/index'
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


# rule kallisto_quant:
#     input:
#         fastq = [
#             'trimmed/{accession}_1.fastq',
#             'trimmed/{accession}_2.fastq'
#         ],
#         index = 'references/reference.idx'
#     output:
#         directory('pseudo_alignment/{accession}/')
#     params:
#         extra = ''
#     log:
#         'logs/kallisto_quant_{accession}.log'
#     wrapper:
#         '0.51.2/bio/kallisto/quant'
rule bwa_mem:
    input:
        reads = ['trimmed/{accession}_1.fastq', 'trimmed/{accession}_2.fastq'],
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


# rule aggregate_results:
#     input:
#         dname_list = expand('pseudo_alignment/{accession}/', accession=config['samples_pe'])
#     output:
#         fname = 'results/results.csv'
#     run:
#         import os
#         import pandas as pd
#
#         df_list = []
#         for dname in input.dname_list:
#             fname = os.path.join(dname, 'abundance.tsv')
#             tmp = pd.read_csv(fname, sep='\t')
#             tmp['accession'] = dname.split('/')[1]
#             df_list.append(tmp)
#         df = pd.concat(df_list)
#
#         df.to_csv(output.fname, index=False)
rule aggregate_results:
    input:
        fname_list = expand(
            'coverage/coverage.{accession}.csv', accession=config['samples_pe'])
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
