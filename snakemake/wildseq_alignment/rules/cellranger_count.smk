# rules/cellranger_count.smk


rule cellranger_count:
    input:
        fastqs=f'{config["FASTQ_DIR"]}/{SLX_ID}'
    params:
        transcriptome='/mnt/scratcha/jblab/middle01/singlecell_seq2/Homo_sapiens_GRCh38_WILDseq_tdTom_WILDseq',
        cores=20,
        memory=100
    output:
        matrix= f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/filtered_feature_bc_matrix.h5',
        bam = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/possorted_genome_bam.bam'
    shell:
        """
        cellranger count --id={wildcards.sample_index} \
                         --transcriptome={params.transcriptome} \
                         --fastqs={input.fastqs} \
                         --sample={wildcards.sample_index} \
                         --localcores={params.cores} \
                         --localmem={params.memory}
        """
