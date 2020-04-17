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
        'results/results.csv'


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


rule kallisto_index:
    input:
        fasta = srcdir(config['reference'])
    output:
        index = 'references/reference.idx'
    log:
        'logs/kallisto_index.log'
    wrapper:
        '0.51.2/bio/kallisto/index'


rule kallisto_quant:
    input:
        fastq = [
            'trimmed/{accession}_1.fastq',
            'trimmed/{accession}_2.fastq'
        ],
        index = 'references/reference.idx'
    output:
        directory('pseudo_alignment/{accession}/')
    params:
        extra = ''
    log:
        'logs/kallisto_quant_{accession}.log'
    wrapper:
        '0.51.2/bio/kallisto/quant'


rule aggregate_results:
    input:
        dname_list = expand('pseudo_alignment/{accession}/', accession=config['samples_pe'])
    output:
        fname = 'results/results.csv'
    run:
        import os
        import pandas as pd

        df_list = []
        for dname in input.dname_list:
            fname = os.path.join(dname, 'abundance.tsv')
            tmp = pd.read_csv(fname, sep='\t')
            tmp['accession'] = dname.split('/')[1]
            df_list.append(tmp)
        df = pd.concat(df_list)

        df.to_csv(output.fname, index=False)
