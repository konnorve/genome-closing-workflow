rule pb_barcode_dataset:
    input:
        raw_barcode_path,
    output:
        pb_barcode_dir / "barcode.dataset.xml",
    resources:
        mem_mb=125000,
    threads: 5
    shell:
        "dataset create --generateIndices --name my_barcodes --type BarcodeSet {output} {input}"

rule pb_demultiplexing:
    input:
        barcode_dataset=pb_barcode_dir / "barcode.dataset.xml",
        subreadset=raw_reads_path,
    output:
        demul_dir = directory(pb_demultiplexing_dir),
        demul_samples = expand(pb_demultiplexing_dir / "outputs" / "demultiplex.{barcode_sample}--{barcode_sample}.subreadset.xml", barcode_sample=BARCODE_SAMPLES),
    resources:
        mem_mb=125000,
    threads: 5
    shell:
        "pbcromwell run pb_demux_subreads -e {input.subreadset} -e {input.barcode_dataset} -n {threads} --output-dir {output.demul_dir} --overwrite"

rule pb_assembly:
    input:
        demul_sample = pb_demultiplexing_dir / "outputs" / "demultiplex.{barcode_sample}--{barcode_sample}.subreadset.xml",
    output:
        out_dir = directory(pb_assembly_dir / "{barcode_sample}"),
        reference_genome = pb_assembly_dir / "{barcode_sample}" / "outputs" / "assembly.rotated.polished.renamed.fsa",
        assembly_report = pb_assembly_dir / "{barcode_sample}" / "outputs" / "polished_assembly.report.json",
        done_flag = touch(done_file_dir / "pb_assembly.{barcode_sample}.done"),
    resources:
        mem_mb=250000,
    threads: 10
    shell:
        """
        pbcromwell run pb_assembly_microbial \
            -e {input.demul_sample} \
            -n {threads} \
            --output-dir {output.out_dir} \
            --overwrite \
            --task-option microasm_genome_size=2000000 \
            --task-option microasm_coverage=30
        """

rule pb_genome_dataset:
    input:
        pb_assembly_dir / "{barcode_sample}" / "outputs" / "assembly.rotated.polished.renamed.fsa",
    output:
        pb_assembly_dir / "{barcode_sample}" / "outputs" / "reference.dataset.xml",
    resources:
        mem_mb=125000,
    threads: 5
    shell:
        "dataset create --generateIndices --name {wildcards.barcode_sample}_ref --type ReferenceSet {output} {input}"


rule pb_methylation:
    input:
        demul_sample = pb_demultiplexing_dir / "outputs" / "demultiplex.{barcode_sample}--{barcode_sample}.subreadset.xml",
        reference_xml = pb_assembly_dir / "{barcode_sample}" / "outputs" / "reference.dataset.xml",
    output:
        out_dir = directory(pb_methylation_dir / "{barcode_sample}"),
        done_flag = touch(done_file_dir / "pb_methylation.{barcode_sample}.done"),
    resources:
        mem_mb=125000,
    threads: 5
    shell:
        "pbcromwell run pb_basemods -e {input.demul_sample} -e {input.reference_xml} -n {threads} --output-dir {output.out_dir} --overwrite"


# rule lima_demultiplexing:
#     input:
#         barcodes=raw_barcode_path,
#         subreadset=raw_reads_path,
#         out_dir=directory(lima_dir)
#     output:
#         barcoded_samples = expand(lima_dir / "lima_demultiplexed.{barcode_sample}--{barcode_sample}.bam", barcode_sample=BARCODE_SAMPLES),
#     resources:
#         mem_mb=125000,
#     threads: 5
#     log: lima_log_dir / "lima.log"
#     shell:
#         "lima --same --peek-guess --split-named --log-file {log} --log-level DEBUG --num-threads {threads} {input.subreadset} {input.barcodes} ${input.out_dir}/lima_demultiplexed.bam"

# rule bam2fastq:
#     input:
#         bam = lima_dir / "lima_demultiplexed.{barcode_sample}--{barcode_sample}.bam"
#     output:
#         fastq = lima_dir / "lima_demultiplexed.{barcode_sample}.fastq"
#     resources:
#         mem_mb=125000,
#     threads: 5
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         "samtools fastq --threads {threads} {input.bam} > {output.fastq}"
        

# rule flye_assembly:
#     input:
#         subreadset = lima_dir / "lima_demultiplexed.{barcode_sample}.fastq"
#     output:
#         out_dir = directory(flye_dir / "{barcode_sample}"),
#         done_flag = touch(done_file_dir / "flye_assembly.{barcode_sample}.done"),
#     resources:
#         mem_mb=250000,
#     threads: 10
#     conda:
#         "../envs/flye.yaml"
#     shell:
#         "flye --threads {threads} -pacbio-raw {input.subreadset} --out-dir {output.out_dir}"

rule make_post_assembly_report:
    input:
        sample_reports = expand(pb_assembly_dir / "{barcode_sample}" / "outputs" / "polished_assembly.report.json", barcode_sample=BARCODE_SAMPLES),
        s_ids_tsv = strain_ids_tsv
    output:
        dataframe_path = results_dir / "genome_df.tsv",
    resources:
        mem_mb=5000,
    conda:
        "../envs/post_assembly_analysis.yaml"
    script:
        "../scripts/post_assembly.py"

rule pb_gorg_classifier:
    input:
        pb_assembly_dir / "{barcode_sample}" / "outputs" / "assembly.rotated.polished.renamed.fsa"
    output:
        out_dir = directory(pb_gorg_dir / "{barcode_sample}"),
        done_flag = touch(done_file_dir / "pb_classification.{barcode_sample}.done"),
    resources:
        mem_mb=250000,
    threads: 10
    shell:
        "nextflow run BigelowLab/gorg-classifier -profile conda --seqs {input} --outdir {output.out_dir} --cpus {threads}"

