configfile: 'data_survey.config.yaml'


rule all:
    input:
        expand('samples/{accession}/pseudo_alignment/{accession}/', accession=config['samples'])


rule download:
    output:
        temp('samples/{accession}/19700101/raw_data/info.tsv')
    log:
        'logs/download.{accession}.log'
    shell:
        """
        # TODO: make all of this more robust
        outdir="$(dirname {output[0]})"

        # download FastQ
        fasterq-dump \
            --outdir "$outdir" \
            {wildcards.accession} \
            &> "{log}"

        # make V-pipe recognize output files
        rename 's/_1/_R1/' "$outdir"/*.fastq
        rename 's/_2/_R2/' "$outdir"/*.fastq

        # get fastq filename
        fname=$(ls "$outdir"/*.fastq | head -1)
        echo "Selected filename: $fname" &>> "{log}"

        # get read length
        read_length=$(bioawk -c fastx '{{ bases += length($seq); count++ }} END{{print int(bases/count)}}' $fname)
        echo "Computed read length: $read_length" &>> "{log}"

        # store in info file
        echo -e "{wildcards.accession}\t19700101\t$read_length" > {output}
        """


rule create_samples_tsv:
    input:
        expand(
            'samples/{accession}/19700101/raw_data/info.tsv',
            accession=config['samples']
        )
    output:
        'samples.tsv'
    shell:
        """
        cat {input} > {output}
        """


rule vpipe_trim:
    input:
        'samples.tsv'
    output:
        'samples/{accession}/19700101/preprocessed_data/R1.fastq.gz',
        'samples/{accession}/19700101/preprocessed_data/R2.fastq.gz'
    shell:
        """
        snakemake -s vpipe.snake --cores 1 alltrimmed
        """


rule kallisto_index:
    input:
        fasta = config['reference']
    output:
        index = 'references/reference.idx'
    log:
        'logs/kallisto_index.log'
    wrapper:
        '0.51.2/bio/kallisto/index'


rule kallisto_quant:
    input:
        fastq = [
            'samples/{accession}/19700101/preprocessed_data/R1.fastq.gz',
            'samples/{accession}/19700101/preprocessed_data/R2.fastq.gz'
        ],
        index = 'references/reference.idx'
    output:
        directory('samples/{accession}/pseudo_alignment/{accession}/')
    params:
        extra = ''
    log:
        'logs/kallisto_quant_{accession}.log'
    threads: 1
    wrapper:
        '0.51.2/bio/kallisto/quant'
