
'''
Snakefile for generating simulated RNA-seq reads by Polyester

Usage:
    snakemake -s SnakePolyester.smk --configfile config.yaml --cores <int> --use-singularity --singularity-args "--bind $HOME:$HOME" --rerun-incomplete
'''

import os
import sys
from pathlib import Path
import pandas as pd
workdir: config["workdir"]
base_dir = os.path.dirname(workflow.snakefile)
# Make list of transcript of interest
transcript_of_interest_list = pd.read_csv(config["transcript_table"], header = None, sep = "\t")[0].tolist()
# Make list of scenarios
scenarios = config["scenarios"]
scenario_list = list(scenarios.keys())
coverage_dict = {scenario: scenarios[scenario]["coverage"] for scenario in scenario_list}
replicates_dict = {scenario: scenarios[scenario]["replicates"] for scenario in scenario_list}
size_factor_dict = {scenario: scenarios[scenario]["size_factor"] for scenario in scenario_list}

rule all:
    input:
        expand("data/{scenario}/flag.txt", scenario = scenario_list)

rule makeTranscriptFasta:
    container:
        "docker://pegi3s/seqkit:2.4.0"
    output:
        transcript_cdna_fasta = "data/transcript_cdna.fa"
    params:
        transcript_of_interest = " ".join(transcript_of_interest_list)
    threads:
        workflow.cores
    benchmark:
        "benchmark/makeTranscriptFasta.txt"
    log:
        "logs/makeTranscriptFasta.log"
    shell:
        """
        seqkit faidx \
        {config[cdna_fasta]} \
        -r {params.transcript_of_interest} 2> {log} | \
        seqkit sort --quiet \
        > {output.transcript_cdna_fasta} \
        2>> {log}
        """

for scenario in scenario_list:

    rule:
        name:
            f"all_{scenario}"
        input:
            bam = expand(f"data/{scenario}/star/sample_{{replicate}}/sample_{{replicate}}_Aligned.out.bam", replicate = [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]),
            bam_table = f"data/{scenario}/bam_table.tsv",
            fastq_R1 = expand(f"data/{scenario}/polyester/sample_{{replicate}}_1.fastq.gz", replicate = [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]),
            fastq_R2 = expand(f"data/{scenario}/polyester/sample_{{replicate}}_2.fastq.gz", replicate = [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]),
            fastq_table = f"data/{scenario}/fastq_table.tsv"
        output:
            flag = f"data/{scenario}/flag.txt"
        shell:
            "touch {output.flag}"

    rule:
        name:
            f"polyester_{scenario}"
        container:
            "docker://naotokubota/snakepolyester:0.1"
        input:
            transcript_cdna_fasta = "data/transcript_cdna.fa"
        output:
            fasta_R1 = temp(expand(f"data/{scenario}/polyester/sample_{{replicate}}_1.fasta", replicate = [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)])),
            fasta_R2 = temp(expand(f"data/{scenario}/polyester/sample_{{replicate}}_2.fasta", replicate = [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]))
        params:
            outdir = f"data/{scenario}/polyester",
            coverage = coverage_dict[scenario],
            replicates = replicates_dict[scenario],
            size_factor = size_factor_dict[scenario],
            base_dir = base_dir
        threads:
            workflow.cores / 4
        benchmark:
            f"benchmark/polyester_{scenario}.txt"
        log:
            f"logs/polyester_{scenario}.log"
        shell:
            """
            Rscript {params.base_dir}/R/Polyester.R \
            --fasta {input.transcript_cdna_fasta} \
            --transcript-table {config[transcript_table]} \
            --coverage {params.coverage} \
            --replicates {params.replicates} \
            --size-factor {params.size_factor} \
            --output {params.outdir} \
            >& {log}
            """

    rule:
        name:
            f"fasta2fastq_{scenario}"
        wildcard_constraints:
            replicate = "|".join([str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]),
        container:
            "docker://staphb/seqtk:1.4"
        input:
            fasta_R1 = f"data/{scenario}/polyester/sample_{{replicate}}_1.fasta",
            fasta_R2 = f"data/{scenario}/polyester/sample_{{replicate}}_2.fasta"
        output:
            fastq_R1 = f"data/{scenario}/polyester/sample_{{replicate}}_1.fastq.gz",
            fastq_R2 = f"data/{scenario}/polyester/sample_{{replicate}}_2.fastq.gz"
        threads:
            1
        benchmark:
            f"benchmark/fasta2fastq_{scenario}_{{replicate}}.txt"
        log:
            f"logs/fasta2fastq_{scenario}_{{replicate}}.log"
        shell:
            """
            seqtk seq -F 'I' {input.fasta_R1} | gzip -c > {output.fastq_R1} 2> {log} && \
            seqtk seq -F 'I' {input.fasta_R2} | gzip -c > {output.fastq_R2} 2>> {log}
            """

    rule:
        name:
            f"star_{scenario}"
        wildcard_constraints:
            replicate = "|".join([str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]),
        container:
            "docker://quay.io/biocontainers/star:2.7.11a--h0033a41_0"
        input:
            fastq_R1 = f"data/{scenario}/polyester/sample_{{replicate}}_1.fastq.gz",
            fastq_R2 = f"data/{scenario}/polyester/sample_{{replicate}}_2.fastq.gz"
        output:
            sam = temp(f"data/{scenario}/star/sample_{{replicate}}/sample_{{replicate}}_Aligned.out.sam")
        params:
            outdir = f"data/{scenario}/star/sample_{{replicate}}/sample_{{replicate}}_",
        threads:
            workflow.cores / 2
        benchmark:
            f"benchmark/star_{scenario}_{{replicate}}.txt"
        log:
            f"logs/star_{scenario}_{{replicate}}.log"
        shell:
            """
            STAR \
            --runThreadN {threads} \
            --genomeDir {config[star_index]} \
            --outFilterMultimapNmax 1 \
            --readFilesIn {input.fastq_R1} {input.fastq_R2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.outdir} \
            >& {log}
            """

    rule:
        name:
            f"sort_{scenario}"
        wildcard_constraints:
            replicate = "|".join([str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]),
        container:
            "docker://quay.io/biocontainers/samtools:1.18--h50ea8bc_1"
        input:
            sam = f"data/{scenario}/star/sample_{{replicate}}/sample_{{replicate}}_Aligned.out.sam"
        output:
            bam = f"data/{scenario}/star/sample_{{replicate}}/sample_{{replicate}}_Aligned.out.bam"
        threads:
            workflow.cores / 4
        benchmark:
            f"benchmark/sort_{scenario}_{{replicate}}.txt"
        log:
            f"logs/sort_{scenario}_{{replicate}}.log"
        shell:
            """
            samtools sort -@ {threads} -O bam -o {output.bam} {input.sam} &> {log} && \
            samtools index {output.bam} &>> {log}
            """

    rule:
        name:
            f"makeBamTable_{scenario}"
        input:
            bam = expand(f"data/{scenario}/star/sample_{{replicate}}/sample_{{replicate}}_Aligned.out.bam", replicate = [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]),
        output:
            bam_table = f"data/{scenario}/bam_table.tsv"
        params:
            replicate = lambda wildcards: [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)],
            scenario = scenario
        threads:
            1
        benchmark:
            f"benchmark/makeBamTable_{scenario}.txt"
        log:
            f"logs/makeBamTable_{scenario}.log"
        run:
            Ref_group = params.replicate[:replicates_dict[params.scenario]]
            Alt_group = params.replicate[replicates_dict[params.scenario]:]
            with open(output.bam_table, "w") as f:
                f.write("sample\tbam\tgroup\n")
                for sample, bam in zip(params.replicate, input.bam):
                    # Absolute path to bam
                    bam = Path(bam).resolve()
                    if sample in Ref_group:
                        f.write(f"sample_{sample}\t{bam}\tRef\n")
                    elif sample in Alt_group:
                        f.write(f"sample_{sample}\t{bam}\tAlt\n")

    '''bam_table
    sample	bam	group
    sample_01	/path/to/sample_01_Aligned.out.bam	Ref
    sample_02	/path/to/sample_02_Aligned.out.bam	Ref
    sample_03	/path/to/sample_03_Aligned.out.bam	Ref
    sample_04	/path/to/sample_04_Aligned.out.bam	Alt
    sample_05	/path/to/sample_05_Aligned.out.bam	Alt
    sample_06	/path/to/sample_06_Aligned.out.bam	Alt
    '''

    rule:
        name:
            f"makeFastqTable_{scenario}"
        input:
            fastq_R1 = expand(f"data/{scenario}/polyester/sample_{{replicate}}_1.fastq.gz", replicate = [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)]),
            fastq_R2 = expand(f"data/{scenario}/polyester/sample_{{replicate}}_2.fastq.gz", replicate = [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)])
        output:
            fastq_table = f"data/{scenario}/fastq_table.tsv"
        params:
            replicate = lambda wildcards: [str(x).zfill(2) for x in range(1, replicates_dict[scenario]*2 + 1)],
            scenario = scenario
        threads:
            1
        benchmark:
            f"benchmark/makeFastqTable_{scenario}.txt"
        log:
            f"logs/makeFastqTable_{scenario}.log"
        run:
            Ref_group = params.replicate[:replicates_dict[params.scenario]]
            Alt_group = params.replicate[replicates_dict[params.scenario]:]
            with open(output.fastq_table, "w") as f:
                f.write("sample\tfastq\tgroup\n")
                for sample, fastq_R1, fastq_R2 in zip(params.replicate, input.fastq_R1, input.fastq_R2):
                    # Absolute path to fastq
                    fastq_R1 = Path(fastq_R1).resolve()
                    fastq_R2 = Path(fastq_R2).resolve()
                    if sample in Ref_group:
                        f.write(f"sample_{sample}\t{fastq_R1},{fastq_R2}\tRef\n")
                    elif sample in Alt_group:
                        f.write(f"sample_{sample}\t{fastq_R1},{fastq_R2}\tAlt\n")

    '''fastq_table
    sample	fastq	group
    sample_01	/path/to/sample_01_1.fastq.gz,/path/to/sample_01_2.fastq.gz	Ref
    sample_02	/path/to/sample_02_1.fastq.gz,/path/to/sample_02_2.fastq.gz	Ref
    sample_03	/path/to/sample_03_1.fastq.gz,/path/to/sample_03_2.fastq.gz	Ref
    sample_04	/path/to/sample_04_1.fastq.gz,/path/to/sample_04_2.fastq.gz	Alt
    sample_05	/path/to/sample_05_1.fastq.gz,/path/to/sample_05_2.fastq.gz	Alt
    sample_06   /path/to/sample_06_1.fastq.gz,/path/to/sample_06_2.fastq.gz Alt
    '''
