configfile: 'data_survey.config.yaml'
workdir: config['workdir']

localrules: all, aggregate_results, plot_results, select_samples


def len_cutoff(wildcards, trim_percent_cutoff=.8):
    """Compute minimal length based on average read length."""
    import glob

    import numpy as np
    from Bio import SeqIO

    fname = sorted(glob.glob(f'data/{wildcards.accession}*.fastq'))[0]

    read_len = np.mean(
        [len(r.seq) for r in SeqIO.parse(fname, 'fastq')]
    ).astype(int)

    len_cutoff = int(trim_percent_cutoff * read_len)
    return len_cutoff


rule all:
    input:
        'plots/',
        'results/selected_samples.csv'


rule download_fastq:
    output:
        fname_marker = touch('data/{accession}.marker')
    log:
        outfile = 'logs/download.{accession}.out.log',
        errfile = 'logs/download.{accession}.err.log'
    resources:
        mem_mb = 5_000
    threads: 6
    group: 'data_processing'
    run:
        import os

        outdir = os.path.dirname(output.fname_marker)
        tmpdir = os.path.join(outdir, f'tmp.{wildcards.accession}')

        shell(
            'fasterq-dump --threads {threads} --outdir {outdir} --temp {tmpdir} {wildcards.accession} > {log.outfile} 2> {log.errfile}'
        )


rule vpipe_trim:
    input:
        fname_marker = 'data/{accession}.marker'
    output:
        touch('trimmed/{accession}.marker')
    params:
        extra = '-ns_max_n 4 -min_qual_mean 30 -trim_qual_left 30 -trim_qual_right 30 -trim_qual_window 10',
        len_cutoff = len_cutoff
    log:
        outfile = 'logs/trimming.{accession}.out.log',
        errfile = 'logs/trimming.{accession}.err.log'
    group: 'data_processing'
    conda:
        'envs/preprocessing.yaml'
    shell:
        """
        echo "The length cutoff is: {params.len_cutoff}" > {log.outfile}

        # detect SE/PE read type
        filecount=$(ls data/{wildcards.accession}*.fastq | wc -l)
        case $filecount in
            1)
                # SE reads
                echo "Read type: SE" > {log.outfile}
                input_spec="-fastq data/{wildcards.accession}.fastq"
                ;;
            2)
                # PE reads
                echo "Read type: PE" > {log.outfile}
                input_spec="-fastq data/{wildcards.accession}_1.fastq -fastq2 data/{wildcards.accession}_2.fastq"
                ;;
            *)
                # oh no
                exit 1
                ;;
        esac

        # do trimming
        prinseq-lite.pl \
            $input_spec \
            {params.extra} \
            -out_format 3 \
            -out_good trimmed/{wildcards.accession} \
            -out_bad null \
            -min_len {params.len_cutoff} \
            -log {log.outfile} 2> >(tee {log.errfile} >&2)

        # remove singletons
        rm -f trimmed/{wildcards.accession}*_singletons.fastq

        # if no reads survive, create empty fastq file
        # (such that mapping does not fail)
        trimmedfilecount=$(shopt -s nullglob; files=(trimmed/{wildcards.accession}*.fastq); echo ${{#files[@]}})
        if [ "$trimmedfilecount" -eq "0" ]; then
            echo "No non-singletons survived trimming, creating empty FastQ file" > {log.outfile}
            touch trimmed/{wildcards.accession}.fastq
        fi
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
        fname_marker = 'trimmed/{accession}.marker',
        index = 'references/reference.amb'
    output:
        fname_bam = 'alignment/{accession}.bam'
    log:
        outfile = 'logs/alignment.{accession}.out.log',
        errfile = 'logs/alignment.{accession}.err.log'
    params:
        index = 'references/reference',
        sort = 'samtools',
        sort_order = 'coordinate',
    resources:
        mem_mb = 16_000
    threads: 8
    group: 'data_processing'
    shell:
        """
        (bwa mem \
            -t {threads} \
            {params.index} \
            trimmed/{wildcards.accession}*.fastq \
            | samtools sort -o {output.fname_bam}) \
        > {log.outfile} 2> {log.errfile}
        """


rule samtools_index:
    input:
        'alignment/{accession}.bam'
    output:
        'alignment/{accession}.bam.bai'
    group: 'data_processing'
    wrapper:
        '0.51.2/bio/samtools/index'


rule compute_coverage:
    input:
        fname = 'alignment/{accession}.bam',
        index = 'alignment/{accession}.bam.bai'
    output:
        fname = 'coverage/coverage.{accession}.csv'
    group: 'data_processing'
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
            wildcards.accession: coverage
        }).to_csv(output.fname, index=False)


rule aggregate_results:
    input:
        fname_list = expand(
            'coverage/coverage.{accession}.csv',
            accession=config['samples_se'] + config['samples_pe'])
    output:
        fname = 'results/coverage.csv',
        fname_stats = 'results/statistics.csv',
        fname_lowquar = 'results/coverage_lowerquartile.csv',
        fname_median = 'results/coverage_median.csv',
        fname_upperquar = 'results/coverage_upperquartile.csv'
    run:
        import pandas as pd

        df_list = []
        for fname in input.fname_list:
            tmp = pd.read_csv(fname)
            df_list.append(tmp)
        df = pd.concat(df_list, axis=1)

        # save data and basic statistics
        df.to_csv(output.fname, index=False)
        df.describe().to_csv(output.fname_stats)

        # save useful metrics
        (pd.melt(df).groupby('variable')
                    .quantile(q=.25)
                    .sort_values('value')
                    .to_csv(output.fname_lowquar))
        (pd.melt(df).groupby('variable')
                    .quantile(q=.5)
                    .sort_values('value')
                    .to_csv(output.fname_median))
        (pd.melt(df).groupby('variable')
                    .quantile(q=.75)
                    .sort_values('value')
                    .to_csv(output.fname_upperquar))


rule plot_results:
    input:
        fname_lowquar = 'results/coverage_lowerquartile.csv',
        fname_median = 'results/coverage_median.csv',
        fname_upperquar = 'results/coverage_upperquartile.csv'
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

        df_lq = pd.read_csv(input.fname_lowquar, index_col=0)
        df_mq = pd.read_csv(input.fname_median, index_col=0)
        df_uq = pd.read_csv(input.fname_upperquar, index_col=0)

        # plot data
        accession_order = df_mq.sort_values('value').index.tolist()

        plt.figure(figsize=(8, 6))

        plt.plot(
            accession_order,
            df_mq.loc[accession_order, 'value'],
            label='median')
        plt.plot(
            accession_order,
            df_lq.loc[accession_order, 'value'],
            alpha=.5,
            label='lower quartile')
        plt.plot(
            accession_order,
            df_uq.loc[accession_order, 'value'],
            alpha=.5,
            label='upper quartile')

        plt.axhline(
            config['thresholds']['median_minium'],
            color='black', ls='dashed', alpha=0.2)
        plt.axhline(
            config['thresholds']['quartile_range'][0],
            color='black', ls='dashed', alpha=0.2)
        plt.axhline(
            config['thresholds']['quartile_range'][1],
            color='black', ls='dashed', alpha=0.2)

        plt.xlabel('SRA Accession')
        plt.ylabel('Per base read count')
        plt.yscale('log')
        plt.tick_params(
            axis='x',
            which='both',
            labelbottom=False)

        plt.legend(loc='best')

        plt.tight_layout()
        plt.savefig(dname / 'coverage.pdf')


rule select_samples:
    input:
        fname_lowquar = 'results/coverage_lowerquartile.csv',
        fname_median = 'results/coverage_median.csv',
        fname_upperquar = 'results/coverage_upperquartile.csv'
    output:
        fname = 'results/selected_samples.csv'
    run:
        import pandas as pd

        # read data
        df_lq = pd.read_csv(input.fname_lowquar, index_col=0)
        df_mq = pd.read_csv(input.fname_median, index_col=0)
        df_uq = pd.read_csv(input.fname_upperquar, index_col=0)

        # apply thresholds
        selection_mq = set(df_mq[
            df_mq['value'] >= config['thresholds']['median_minium']
        ].index)
        selection_lq = set(df_lq[
            df_lq['value'] >= config['thresholds']['quartile_range'][0]
        ].index)
        selection_uq = set(df_uq[
            df_uq['value'] <= config['thresholds']['quartile_range'][1]
        ].index)

        selection = selection_mq & selection_lq & selection_uq

        # save results
        with open(output.fname, 'w') as fd:
            fd.write('accession\n')
            fd.write('\n'.join(selection))
