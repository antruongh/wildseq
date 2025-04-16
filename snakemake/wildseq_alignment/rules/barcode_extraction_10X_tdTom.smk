# rules/barcode_extraction_10X_tdTom.smk

# cellranger_dir = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/mnt/scratcha/jblab/middle01/singlecell_seq2/results/{SLX_ID}.align/{{{sample_index}}}/outs'

rule barcode_extraction_10X_tdTom:
    input:
        bam = rules.cellranger_count.output.bam,
        
    output:
        wildseq_bam = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/WILDseq.bam',
        
        filtered_fa = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/bamtofastq_R2_filtered.fa',
        bamtobamtofastq_dir = directory(f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/bamtofastq'),
        #cut_fa = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/bamtofastq_R2_filtered_cut.fa',
        cut_fa = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/bamtofastq_R2_filtered_cut.fa',
        bowtie_output = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/bamtofastq_R2_theoretical.bowtie',

        unmapped_txt = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/unmapped_theoretical.txt',
        unmapped_count = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/unmapped_theoretical_count.txt',
        cb_ub_txt = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/{{sample_index}}_WILDseq.CB.UB.txt',
    params:
        bowtie_index = "/scratcha/imaxt/sawick01/theoretical_barcodes2/WILDseq_index2",
        samtools = "/home/bioinformatics/software/samtools/samtools-1.9/bin/samtools",
        bamtofastq = "/home/bioinformatics/software/bamtofastq-1.2.0",
        bowtie = "/home/bioinformatics/software/bowtie/bowtie-1.3.0/bowtie",
    shell:
        """
        # Step 1: Extract reads mapping to the specific region
        {params.samtools} view -b {input.bam} EF1a_tdTomato_WILDseq:5398-5441 > {output.wildseq_bam}

        # Step 2: Convert BAM to FASTQ
        {params.bamtofastq} {output.wildseq_bam} {output.bamtofastq_dir}

        # Step 3: Filter reads containing the full barcode sequence
        FILE=$(find {output.bamtofastq_dir} -type f -name "*R2*.fastq.gz" | head -n 1)
        zcat $FILE | egrep -A 2 -B 1 [ATGC]{12}GCTACCTGGTCGTTCTATCG[ATGC]{12} | sed '/^--$/d' | sed -n '1~4s/^@/>/p;2~4p' > {output.filtered_fa}

        # Step 4: Trim flanking regions
        cat {output.filtered_fa} | sed 's/\(.*\)\([ATCG]\{12\}GCTACCTGGTCGTTCTATCG[ATCG]\{12\}\)\(.*\)/\2/g' > {output.cut_fa}

        # Step 5: Align to theoretical barcode index
        {params.bowtie} -f -v 2 -x {params.bowtie_index} --un {output.unmapped_txt} {output.cut_fa} > {output.bowtie_output}

        # Step 6: Count unmapped barcodes
        sort {output.unmapped_txt} | uniq -c | sort -nrk 1 > {output.unmapped_count}

        # Step 7: Extract cell barcode and UMI info

        {params.samtools} view {output.wildseq_bam} | sed -r -n 's/.*(LH00417[:A-Z0-9]*).*(CB:Z:[ACTG]*).*(UB:Z:[ACTG]*).*/\1 \2 \3/p' > {output.cb_ub_txt}
        """
