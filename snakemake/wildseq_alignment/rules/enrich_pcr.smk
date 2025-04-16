# rules/enrich_pcr.smk
rule enrich_pcr:
    input:
        fastq_dir = f'{config["FASTQ_DIR"]}/{SLX_ID}',
        out_dir = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs',
    params:
        bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNNNN",
        bowtie_index = "/scratcha/imaxt/sawick01/theoretical_barcodes2/WILDseq_index2",
        bowtie = "/home/bioinformatics/software/bowtie/bowtie-1.3.0/bowtie",
        enrich_index = lambda wildcards: f'{config["sample_map"][wildcards.sample_index]}',    
    output:
        table_full = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/{{sample_index}}_EnrichPCR_table_full.txt',
        table_full_filtered = f'{config["OUT_DIR"]}/{SLX_ID}/{{sample_index}}/outs/{{sample_index}}_EnrichPCR_table_full_filtered.txt'
    threads: 8
    shell:
        """
        # whitelist of all 10X cell barcodes in original sample
        umi_tools whitelist --stdin={input.fastq_dir}/{wildcards.sample_index}_*_R1*.fastq.gz \
                            --bc-pattern={params.bc_pattern} \
                            --log2stderr > {input.fastq_dir}/whitelist_{wildcards.sample_index}_mRNA_All.txt

        # extract 10x cell barcodes and UMIs from Enriched library based on the whitelist (i.e. if you enriched a barcode that's not in the sample, it will be ignored)

        umi_tools extract --bc-pattern={params.bc_pattern} \
                          --stdin={input.fastq_dir}/{params.enrich_index}*_R1*.fastq.gz \
                          --read2-stdout --read2-in={input.fastq_dir}/{params.enrich_index}*_R2*.fastq.gz \
                          --stdout={input.fastq_dir}/{wildcards.sample_index}_PCR_All_R2_extracted.fastq.gz \
                          --whitelist={input.fastq_dir}/whitelist_{wildcards.sample_index}_mRNA_All.txt

        # extract the part of the read sequence matching WILDseq barcode pattern
        zcat {input.fastq_dir}/{params.enrich_index}_PCR_All_R2_extracted.fastq.gz | \
        egrep -A 2 -B 1 '[ATGC]{{12}}GCTACCTGGTCGTTCTATCG[ATGC]{{12}}' | \
        sed '/^--$/d' | sed -n '1~4s/^@/>/p;2~4p' | \
        sed 's/\\(.*\\)\\([ATCG]\\{{12\\}}GCTACCTGGTCGTTCTATCG[ATCG]\\{{12\\}}\\)\\(.*\\)/\\2/g' > \
        {input.fastq_dir}/{params.enrich_index}_PCR_All_R2_extracted_barcode_sequences_only.fasta

        # ${ENRICH}_PCR_All_R2_extracted_barcode_sequences_only.fasta

        # barcode sequences aligned to library of all barcodes
        {params.bowtie} -f -v 2 --un {input.out_dir}/{params.enrich_index}_unmapped.txt \
            -x {params.bowtie_index} {input.fastq_dir}/{params.enrich_index}_PCR_All_R2_extracted_barcode_sequences_only.fasta > \
            {input.out_dir}/{params.enrich_index}_PCR.bowtie
        
        # /${ORIG}/outs/${ENRICH}_unmapped.txt
        # ${ENRICH}_PCR_All_R2_extracted_barcode_sequences_only.fasta
        # ./${ORIG}/outs/${ENRICH}_PCR.bowtie


        # Count unmapped barcodes
        sort {input.out_dir}/{params.enrich_index}_unmapped.txt | \
        uniq -c | sort -nrk 1  > {input.out_dir}/{params.enrich_index}_unmapped_count.txt

        #./${ORIG}/outs/${ENRICH}_unmapped.txt 
        #./${ORIG}/outs/${ENRICH}_unmapped_count.txt


        # Generate table from bowtie output
        awk -F'[\\t"_ ]' '{print $2 "\t" $3 "\t" $6}' {input.out_dir}/{params.enrich_index}_PCR.bowtie | sort | uniq > {output.table_full}

        # ./${ORIG}/outs/${ENRICH}_PCR.bowtie
        # ./${ORIG}/outs/${ENRICH}_EnrichPCR_table_full.txt

        # Filter table by barcode count

        awk -F'[\\t"_ ]' '{print $2 "\t" $3 "\t" $6}' {input.out_dir}/{params.enrich_index}_PCR.bowtie | sort | uniq -c | sed 's/^[ \t]*//' | awk -F'[\\t" ]' '{print $1 "\t" $2 "\t" $3 "\t" $4}' | \
        awk '{k=$2FS$3; if (a[k]<$1) {a[k]=$1; data[k]=$0} else if (a[k]==$1) {a[k]=$1; data[k]=data[k]","$4}} END {for (i in a) print  data[i]}' > {output.table_full_filtered}
        #./${ORIG}/outs/${ENRICH}_EnrichPCR_table_full_filtered.txt

        """
