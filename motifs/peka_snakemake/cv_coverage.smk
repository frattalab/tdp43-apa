import pandas as pd
import os

# Extract values from sample table based on the wildcard
def get_param(wildcards, param_name):
    sample = wildcards.sample
    if sample in sample_table['comparison_name'].values:
        index = sample_table[sample_table['comparison_name'] == sample].index[0]
        return sample_table.at[index, param_name]
    else:
        raise ValueError(f"Sample {sample} not found in sample_table")


# Read the sample table using pandas
sample_table = pd.read_csv(config["sample_table"])

# Define the list of sample names based on the 'comparison_name' column
sample_names = sample_table['comparison_name'].tolist()

# base output directory
base_outdir = config["base_output_dir"]


include: "peka_beds.smk"

rule all_cv_coverage:
    input:
        expand(os.path.join(base_outdir, "cv_coverage", "{sample}", "cv_coverage_{group}", "cv_coverage.done"), sample=sample_names, group=["foreground", "background"]),
        expand(os.path.join(base_outdir, "cv_coverage", "{sample}", "{kmer}", "cv_coverage_{group}", "cv_coverage.done"), sample=sample_names, kmer = config["cv_coverage_motifs"], group=["foreground", "background"]) if config["run_per_motif"] else [] # rules.pas_to_peka_beds.output

# plot & get % coverage for specific motifs in given window around foreground regions
# final TSV output to ./results/no_of_txn_{sample_name}' 
# where sample_name is input BED file with .bed suffix removed and any text with '.xl.' removed 
# again too lazy to code out
# outputs in current directory, so have to make and switch inside script

rule per_motif_cv_coverage_foreground:
    input:
        bed=rules.pas_to_peka_beds.output.reg_bed,
        genome=os.path.abspath(config["genome"]),
        genome_idx=os.path.abspath(config["genome_idx"])

    output:
         os.path.join(base_outdir, "cv_coverage", "{sample}", "{kmer}", "cv_coverage_foreground", "cv_coverage.done")
    
    params:
        script=os.path.abspath("scripts/cv_coverage.py"),
        bed=os.path.abspath(rules.pas_to_peka_beds.output.reg_bed),
        motifs="{kmer}",
        kmer=config["kmer"],
        percentile=config["percentile"],
        regions=os.path.abspath(config["regions_gtf"]),
        specific_region=config["specific_region"],
        window=config["distal_window"],
        smoothing=config["cv_coverage_smoothing"],
        use_scores=config["cv_coverage_usescores"],
        cap=config["cv_coverage_cap"],
        chunk_size=10000,
        outdir=os.path.join(base_outdir, "cv_coverage", "{sample}", "{kmer}", "cv_coverage_foreground"),
        tmpdir=os.path.abspath(os.path.join(base_outdir, "cv_coverage", "{sample}", "{kmer}", "cv_coverage_foreground", "tmp")),
        out="cv_coverage.done"

    threads:
        config["cv_coverage_n_cores"]


    log:
        stdout=os.path.abspath(os.path.join(base_outdir, "logs", "cv_coverage", "cv_coverage_foreground.{kmer}.{sample}.stdout.txt")),
        stderr=os.path.abspath(os.path.join(base_outdir, "logs", "cv_coverage", "cv_coverage_foreground.{kmer}.{sample}.stderr.txt"))

    conda:
        "envs/peka_fork.yaml"


    shell:
        """
        mkdir -p {params.tmpdir} && export TMPDIR={params.tmpdir} && cd {params.outdir}
        python {params.script} \
        {params.bed} \
        {params.motifs} \
        {params.specific_region} \
        {params.kmer} \
        {input.genome} \
        {input.genome_idx} \
        {params.regions} \
        {params.smoothing} \
        {params.percentile} \
        {params.window} \
        {params.use_scores} \
        {threads} \
        {params.chunk_size} \
        {params.cap} \
        1> {log.stdout} \
        2> {log.stderr} && \
        touch {params.out}
        """


rule per_motif_cv_coverage_background:
    input:
        bed=rules.pas_to_peka_beds.output.bg_bed,
        genome=os.path.abspath(config["genome"]),
        genome_idx=os.path.abspath(config["genome_idx"])
    
    output:
         os.path.join(base_outdir, "cv_coverage", "{sample}", "{kmer}", "cv_coverage_background", "cv_coverage.done")
    
    params:
        script=os.path.abspath("scripts/cv_coverage.py"),
        bed=os.path.abspath(rules.pas_to_peka_beds.output.bg_bed),
        motifs="{kmer}",
        kmer=config["kmer"],
        percentile=config["percentile"],
        regions=os.path.abspath(config["regions_gtf"]),
        specific_region=config["specific_region"],
        window=config["distal_window"],
        smoothing=config["cv_coverage_smoothing"],
        use_scores=config["cv_coverage_usescores"],
        cap=config["cv_coverage_cap"],
        chunk_size=10000,
        outdir=os.path.join(base_outdir, "cv_coverage", "{sample}", "{kmer}", "cv_coverage_background"),
        tmpdir=os.path.abspath(os.path.join(base_outdir, "cv_coverage", "{sample}", "{kmer}", "cv_coverage_background", "tmp")),
        out="cv_coverage.done"

    threads:
        config["cv_coverage_n_cores"]


    log:
        stdout=os.path.abspath(os.path.join(base_outdir, "logs", "cv_coverage", "cv_coverage_background.{kmer}.{sample}.stdout.txt")),
        stderr=os.path.abspath(os.path.join(base_outdir, "logs", "cv_coverage", "cv_coverage_background.{kmer}.{sample}.stderr.txt"))

    conda:
        "envs/peka_fork.yaml"


    shell:
        """
        mkdir -p {params.tmpdir} && export TMPDIR={params.tmpdir} && cd {params.outdir}
        python {params.script} \
        {params.bed} \
        {params.motifs} \
        {params.specific_region} \
        {params.kmer} \
        {input.genome} \
        {input.genome_idx} \
        {params.regions} \
        {params.smoothing} \
        {params.percentile} \
        {params.window} \
        {params.use_scores} \
        {threads} \
        {params.chunk_size} \
        {params.cap} \
        1> {log.stdout} \
        2> {log.stderr} && \
        touch {params.out}
        """


rule cv_coverage_foreground:
    input:
        bed=rules.pas_to_peka_beds.output.reg_bed,
        genome=os.path.abspath(config["genome"]),
        genome_idx=os.path.abspath(config["genome_idx"])

    output:
         os.path.join(base_outdir, "cv_coverage", "{sample}", "cv_coverage_foreground", "cv_coverage.done")
    
    params:
        script=os.path.abspath("scripts/cv_coverage.py"),
        bed=os.path.abspath(rules.pas_to_peka_beds.output.reg_bed),
        motifs=",".join(config["cv_coverage_motifs"]),
        kmer=config["kmer"],
        percentile=config["percentile"],
        regions=os.path.abspath(config["regions_gtf"]),
        specific_region=config["specific_region"],
        window=config["distal_window"],
        smoothing=config["cv_coverage_smoothing"],
        use_scores=config["cv_coverage_usescores"],
        cap=config["cv_coverage_cap"],
        chunk_size=10000,
        outdir=os.path.join(base_outdir, "cv_coverage", "{sample}", "cv_coverage_foreground"),
        tmpdir=os.path.abspath(os.path.join(base_outdir, "cv_coverage", "{sample}", "cv_coverage_foreground", "tmp")),
        out="cv_coverage.done"

    threads:
        config["cv_coverage_n_cores"]


    log:
        stdout=os.path.abspath(os.path.join(base_outdir, "logs", "cv_coverage", "cv_coverage.{sample}.stdout.txt")),
        stderr=os.path.abspath(os.path.join(base_outdir, "logs", "cv_coverage", "cv_coverage.{sample}.stderr.txt"))

    conda:
        "envs/peka_fork.yaml"


    shell:
        """
        mkdir -p {params.tmpdir} && export TMPDIR={params.tmpdir} && cd {params.outdir}
        python {params.script} \
        {params.bed} \
        {params.motifs} \
        {params.specific_region} \
        {params.kmer} \
        {input.genome} \
        {input.genome_idx} \
        {params.regions} \
        {params.smoothing} \
        {params.percentile} \
        {params.window} \
        {params.use_scores} \
        {threads} \
        {params.chunk_size} \
        {params.cap} \
        1> {log.stdout} \
        2> {log.stderr} && \
        touch {params.out}
        """


rule cv_coverage_background:
    input:
        bed=rules.pas_to_peka_beds.output.bg_bed,
        genome=os.path.abspath(config["genome"]),
        genome_idx=os.path.abspath(config["genome_idx"])
    
    output:
         os.path.join(base_outdir, "cv_coverage", "{sample}", "cv_coverage_background", "cv_coverage.done")
    
    params:
        script=os.path.abspath("scripts/cv_coverage.py"),
        bed=os.path.abspath(rules.pas_to_peka_beds.output.bg_bed),
        motifs=",".join(config["cv_coverage_motifs"]),
        kmer=config["kmer"],
        percentile=config["percentile"],
        regions=os.path.abspath(config["regions_gtf"]),
        specific_region=config["specific_region"],
        window=config["distal_window"],
        smoothing=config["cv_coverage_smoothing"],
        use_scores=config["cv_coverage_usescores"],
        cap=config["cv_coverage_cap"],
        chunk_size=10000,
        outdir=os.path.join(base_outdir, "cv_coverage", "{sample}", "cv_coverage_background"),
        tmpdir=os.path.abspath(os.path.join(base_outdir, "cv_coverage", "{sample}", "cv_coverage_background", "tmp")),
        out="cv_coverage.done"

    threads:
        config["cv_coverage_n_cores"]


    log:
        stdout=os.path.abspath(os.path.join(base_outdir, "logs", "cv_coverage", "cv_coverage.{sample}.stdout.txt")),
        stderr=os.path.abspath(os.path.join(base_outdir, "logs", "cv_coverage", "cv_coverage.{sample}.stderr.txt"))

    conda:
        "envs/peka_fork.yaml"


    shell:
        """
        mkdir -p {params.tmpdir} && export TMPDIR={params.tmpdir} && cd {params.outdir}
        python {params.script} \
        {params.bed} \
        {params.motifs} \
        {params.specific_region} \
        {params.kmer} \
        {input.genome} \
        {input.genome_idx} \
        {params.regions} \
        {params.smoothing} \
        {params.percentile} \
        {params.window} \
        {params.use_scores} \
        {threads} \
        {params.chunk_size} \
        {params.cap} \
        1> {log.stdout} \
        2> {log.stderr} && \
        touch {params.out}
        """