# To run: snakemake --cores <number_of_cores>

import os

configfile: "config.yaml"

TISSUES = list(config["tissues"].keys())
SAMPLES = [sample for tissue_samples in config["tissues"].values() for sample in tissue_samples]

rule all:
    input:
        expand("results/01_isoquant_raw/{sample}", sample=SAMPLES),
        expand("results/02_merged_readinfo/{tissue}.done", tissue=TISSUES),
        expand("results/03_polya_matrix/{tissue}/{tissue}_cpm_polya_gene.csv", tissue=TISSUES),
        expand("results/04_tagged_bams/{sample}/{sample}.elongating.bam.bai", sample=SAMPLES),
        expand("results/05_extracted_polya/{tissue}/{tissue}_dataarray.pkl", tissue=TISSUES),
        "results/06_merged_data/merged_info_df.pkl",
        "results/06_merged_data/merged_array.npy",
        "results/07_distance_matrix/merge_dis_complete.h5",

rule run_isoquant:
    input:
        bam=os.path.join(config["paths"]["bam_dir"], "{sample}.rmRNA.bam")
    output:
        directory("results/01_isoquant_raw/{sample}")
    threads: 4
    shell:
        """
        isoquant.py \\
            -d nanopore \\
            -o {output} \\
            -p {threads} \\
            --bam {input.bam} \\
            --reference {config[paths][ref_genome]} \\
            --genedb {config[paths][gene_db]} \\
            --prefix {wildcards.sample} \\
            --labels {wildcards.sample} \\
            --force \\
            --complete_genedb \\
            --report_novel_unspliced true \\
            --check_canonical \\
            --sqanti_output
        """

rule merge_isoquant2polya:
    input:
        isoquant_dirs=lambda wildcards: expand("results/01_isoquant_raw/{sample}", sample=config["tissues"][wildcards.tissue])
    output:
        flag_file="results/02_merged_readinfo/{tissue}.done"
    threads: config["params"]["merge_threads"]
    run:
        replicates = " ".join(config["tissues"][wildcards.tissue])
        shell(f"""
            python3 {config[paths][utils_dir]}/combine_isoform_v3.py \\
                --replicates {replicates} \\
                --isoquant_dir "results/01_isoquant_raw" \\
                --gtf_dir "results/02_discovery_gtf/{wildcards.tissue}" \\
                --output_dir "results/02_merged_readinfo/{wildcards.tissue}" \\
                --results_dir "results/02_replicates_results/{wildcards.tissue}" \\
                --num_processes {threads}

            touch {output.flag_file}
        """)

rule calculate_polya_matrix:
    input:
        merge_done=rules.merge_isoquant2polya.output.flag_file,
        isoquant_dirs=lambda wildcards: expand("results/01_isoquant_raw/{sample}", sample=config["tissues"][wildcards.tissue])
    output:
        cpm_csv="results/03_polya_matrix/{tissue}/{tissue}_cpm_polya_gene.csv"
    threads: 1
    run:
        merged_data = expand("results/02_merged_readinfo/{wildcards.tissue}/{sample}.parquet", sample=config["tissues"][wildcards.tissue])
        isoquant_paths = expand("results/01_isoquant_raw/{sample}", sample=config["tissues"][wildcards.tissue])
        output_dir = os.path.dirname(output.cpm_csv)

        shell(f"""
            mkdir -p {output_dir}
            python3 {config[paths][utils_dir]}/gene_polya_matrix_v2.py \\
                --merged_data {" ".join(merged_data)} \\
                --isoquant_path {" ".join(isoquant_paths)} \\
                --label {wildcards.tissue} \\
                --output_dir {output_dir}
        """)

rule combine_cpm_matrices:
    input:
        expand("results/03_polya_matrix/{tissue}/{tissue}_cpm_polya_gene.csv", tissue=TISSUES)
    output:
        csv="results/03_polya_matrix/combined_cpm_polya_matrix.csv"
    params:
        input_dir="results/03_polya_matrix"
    shell:
        """
        python3 scripts/utils/concat_cpm.py \\
            --input_dir {params.input_dir} \\
            --output_file {output.csv}
        """

rule add_tag:
    input:
        bam=os.path.join(config["paths"]["bam_dir"], "{sample}.rmRNA.bam"),
        parquet=lambda wildcards: f"results/02_merged_readinfo/{[t for t, s in config['tissues'].items() if wildcards.sample in s][0]}/{wildcards.sample}.parquet",
        merge_done=lambda wildcards: f"results/02_merged_readinfo/{[t for t, s in config['tissues'].items() if wildcards.sample in s][0]}.done",
    output:
        tagged_bam="results/04_tagged_bams/{sample}/{sample}.tagged.bam",
        polya_bam="results/04_tagged_bams/{sample}/{sample}.polyadenylated.bam",
        elong_bam="results/04_tagged_bams/{sample}/{sample}.elongating.bam",
        elong_bai="results/04_tagged_bams/{sample}/{sample}.elongating.bam.bai"
    threads: 1
    shell:
        """
        OUT_DIR=$(dirname {output.tagged_bam})
        mkdir -p $OUT_DIR

        python3 {config[paths][utils_dir]}/add_tag_to_bam_iso.py \\
            -i {input.bam} \\
            --read_info {input.parquet} \\
            --adapter_info "{config[paths][aligned_data_dir]}/{wildcards.sample}.adapter.result.txt" \\
            --polya_info "{config[paths][aligned_data_dir]}/{wildcards.sample}.polyA_tail.result.txt" \\
            -o {output.tagged_bam}
        samtools index {output.tagged_bam}

        python3 {config[paths][utils_dir]}/get_polyadenylated_reads.py -i {output.tagged_bam} -o {output.polya_bam}
        samtools index {output.polya_bam}

        python3 {config[paths][utils_dir]}/get_elongating_reads.py -i {output.tagged_bam} -o {output.elong_bam}
        samtools index {output.elong_bam}
        """

rule extract_polya:
    input:
        merge_done=rules.merge_isoquant2polya.output.flag_file
    output:
        data_array="results/05_extracted_polya/{tissue}/{tissue}_dataarray.pkl"
    threads: config["params"]["polya_threads"]
    run:
        merged_data = expand("results/02_merged_readinfo/{wildcards.tissue}/{sample}.parquet", sample=config["tissues"][wildcards.tissue])
        gtf_paths = expand("results/01_isoquant_raw/{sample}/{sample}.extended_annotation.gtf", sample=config["tissues"][wildcards.tissue])
        output_dir = os.path.dirname(output.data_array)
        
        shell(f"""
            mkdir -p {output_dir}
            python3 {config[paths][utils_dir]}/process_dataarray_v2.py \\
                --merged_data {" ".join(merged_data)} \\
                --gtf_path {" ".join(gtf_paths)} \\
                --label {wildcards.tissue} \\
                --output_dir {output_dir} \\
                --bin_width {config[params][polya_bin_width]} \\
                --thread {threads}
        """)

rule merge_pickled_outputs:
    input:
        pkls=expand("results/05_extracted_polya/{tissue}/{tissue}_dataarray.pkl", tissue=TISSUES)
    output:
        info_df="results/06_merged_data/merged_info_df.pkl",
        array="results/06_merged_data/merged_array.npy"
    threads: 1
    shell:
        """
        python3 scripts/preprocessing_pipeline/06_merge_pickles/merge_pickles.py \\
            --input_dir "results/05_extracted_polya" \\
            --output_info {output.info_df} \\
            --output_array {output.array}
        """

rule calculate_polya_cophenetic_pearson:
    input:
        info_df="results/06_merged_data/merged_info_df.pkl",
        array="results/06_merged_data/merged_array.npy"
    output:
        h5_matrix="results/07_distance_matrix/polya_cophenetic_pearson.h5"
    params:
        bin_width=config["params"]["distance_matrix"]["polya_bin_width"],
        block_size=config["params"]["distance_matrix"]["block_size"]
    threads: config["params"]["distance_matrix"]["threads"]
    shell:
        """
        python3 scripts/preprocessing_pipeline/07_distance_matrix/polya_cophenet_distance.py \\
            --info_df {input.info_df} \\
            --array {input.array} \\
            --output_h5 {output.h5_matrix} \\
            --bin_width {params.bin_width} \\
            --block_size {params.block_size} \\
            --threads {threads}
        """

rule calculate_polya_direct_pearson:
    input:
        info_df="results/06_merged_data/merged_info_df.pkl",
        array="results/06_merged_data/merged_array.npy"
    output:
        h5_matrix="results/07_distance_matrix/polya_direct_pearson.h5"
    params:
        block_size=config["params"]["distance_matrix"]["block_size"]
    threads: config["params"]["distance_matrix"]["threads"]
    shell:
        """
        python3 scripts/preprocessing_pipeline/07_distance_matrix/polya_pearson_distance.py \\
            --info_df {input.info_df} \\
            --array {input.array} \\
            --output_h5 {output.h5_matrix} \\
            --block_size {params.block_size} \\
            --threads {threads}
        """

rule merge_polya_distances:
    input:
        coph_h5=rules.calculate_polya_cophenetic_pearson.output.h5_matrix,
        pearson_h5=rules.calculate_polya_direct_pearson.output.h5_matrix
    output:
        h5_matrix="results/07_distance_matrix/polya_distance.h5"
    params:
        block_size=config["params"]["distance_matrix"]["block_size"]
    shell:
        """
        python3 scripts/preprocessing_pipeline/07_distance_matrix/polya_cophenet_distance.py \\
            --polya_h5 {input.coph_h5} \\
            --expression_h5 {input.pearson_h5} \\
            --output_h5 {output.h5_matrix} \\
            --block_size {params.block_size}
        """

rule calculate_expression_distance:
    input:
        info_df="results/06_merged_data/merged_info_df.pkl",
        cpm_csv=rules.combine_cpm_matrices.output.csv
    output:
        h5_matrix="results/07_distance_matrix/expression_distance.h5"
    params:
        block_size=config["params"]["distance_matrix"]["block_size"]
    threads: config["params"]["distance_matrix"]["threads"]
    shell:
        """
        python3 scripts/preprocessing_pipeline/07_distance_matrix/distance_utils.py \\
            --info_df {input.info_df} \\
            --cpm_matrix {input.cpm_csv} \\
            --output_h5 {output.h5_matrix} \\
            --block_size {params.block_size} \\
            --threads {threads}
        """

