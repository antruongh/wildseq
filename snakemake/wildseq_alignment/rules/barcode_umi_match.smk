# rules/barcode_umi_match.smk

rule barcode_umi_match:
    input:
        bowtie_output = rules.barcode_extraction_10X_tdTom.output.bowtie_output,
        #f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/bamtofastq_R2_theoretical.bowtie',
        cb_ub_txt = rules.barcode_extraction_10X_tdTom.output.cb_ub_txt, 
        # f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/{{sample_index}}_WILDseq.CB.UB.txt'
    output:
        cb_table = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/{{sample_index}}_CBC_table_full.txt'
    shell:
        """
        Rscript /mnt/scratcha/jblab/truong01/wildseq/snakemake/wildseq_alignment/scripts/barcode_10X.R \
                --sample {{sample_index}} \
                --theoretical_bowtie {input.bowtie_output} \
                --cb_ub_txt {input.cb_ub_txt}\
                --output_file {output.cb_table}
        """
