# rules/barcode_assignment.smk

#sample_metadata_file = config["sample_metadata"]
#sample_metadata = pd.read_csv(sample_metadata_file).set_index("sample_index")
#sample_map = sample_metadata["enrich_index"].to_dict()  # Creates a dictionary for lookup

rule barcode_assignment:
    input:
        #cb_ub_txt = f"{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/{{sample_index}}_WILDseq.CB.UB.txt",
        cb_table = rules.barcode_umi_match.output.cb_table,
        table_full_filtered = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/EnrichPCR_table_full_filtered.txt'
    output:
        barcode_assigned_table = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/barcode_assignment_{SLX_ID}_{{sample_index}}.txt'
    script:
        """
        Rscript /mnt/scratcha/jblab/truong01/wildseq/snakemake/wildseq_alignment/scripts/barcode_assignment.R \
                --sample {{sample_index}} \
    			--cb_table {input.cb_table} \
    			--table_full_filtered {input.table_full_filtered} \
    			--output_file {output.barcode_assigned_table}
        """
