import pandas as pd
import os

window_keys = ["upstream", "downstream"]
strand_keys = {"plus": "+", "minus": "-"}

in_bam_dir = config["input_dir"]
bam_suffix = config["bam_suffix"]
output_dir = config["output_dir"]

sample_names = [f.replace(bam_suffix, "") for f in os.listdir(in_bam_dir) if f.endswith(bam_suffix)]

# Check that each input sample has a BAM index (.bai)
for s in sample_names:
    assert os.path.exists(os.path.join(in_bam_dir, s + bam_suffix + ".bai")), f".bai index file does not exist at same location as input BAM file for sample {s}"

wildcard_constraints:
    window = "|".join(window_keys),
    strand = "|".join(strand_keys.keys())


rule all:
    input:
         expand(os.path.join(output_dir, "coverage", "{sample}.pas_windows.{window}.summarised_coverage.bed"), sample=sample_names,window=window_keys)


rule bam_to_idxstats:
    input:
        bam=os.path.join(in_bam_dir, "{sample}" + bam_suffix)
    output:
        temp(os.path.join(output_dir, "chromsizes", "{sample}.idxstats.txt"))
    log:
        stderr=os.path.join(output_dir, "logs", "bam_to_idxstats.{sample}.stderr.txt")

    container:
        "docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0"

    shell:
        """
        samtools idxstats {input.bam} 1> {output} 2> {log.stderr}
        """

rule idxstats_to_chromsizes:
    input:
        rules.bam_to_idxstats.output
    
    output:
        os.path.join(output_dir, "chromsizes", "{sample}.chromsizes.txt")
    
    log:
        stderr=os.path.join(output_dir, "logs", "idxstats_to_chromsizes.{sample}.stderr.txt")

    shell:
        """
        awk -F"\t" -v OFS="\t" '{{print $1,$2}}' {input} 1> {output} 2> {log.stderr}
        """ 


rule create_pas_windows:
    input: 
        config["pas_bed"]
    output:
        window_bed = expand(os.path.join(output_dir, "regions", "pas_windows.{window}.bed"), window=window_keys),
        merged_bed = os.path.join(output_dir, "regions", "pas_windows.merged.bed")

    params:
        script = "scripts/create_pas_windows.py",
        output_prefix = os.path.join(output_dir, "regions", "pas_windows"),
        window_size = config["window_size"],
        exclude_pas = "--exclude_pas" if config["exclude_pas"] else "",
        return_merged = "--return_merged"

    log:
        stdout = os.path.join(output_dir, "logs", "create_pas_windows.stdout.txt"),
        stderr = os.path.join(output_dir, "logs", "create_pas_windows.stderr.txt")
    
    container:
        "docker://quay.io/biocontainers/pyranges:0.0.120--pyh7cba7a3_0"
    
    shell:
        """
        python {params.script} \
        --window_size {params.window_size} \
        {params.exclude_pas} {params.return_merged} \
        {input} \
        {params.output_prefix} \
        1> {log.stdout} \
        2> {log.stderr}
        """


# Make stranded bw files of per-base coverage over intervals that respects RNA-seq strandedness
rule bam_to_stranded_bigwig: 
    input:
        bam=os.path.join(in_bam_dir, "{sample}" + config["bam_suffix"]),
        bed=rules.create_pas_windows.output.merged_bed,
        chromsizes=rules.idxstats_to_chromsizes.output
    
    output:
        expand(os.path.join(output_dir, "coverage", "{{sample}}.regions.{strand}.bw"), strand=strand_keys.keys()) #double brace to mask + retain sample wildcard

    params:
        script="scripts/bam_to_bw.py",
        output_prefix=os.path.join(output_dir, "coverage", "{sample}.regions"),
        strandedness=config["strandedness"],
        keep_orphans="--keep-orphans" if True else "" # placeholder - hardcoding for this analysis

    log:
        stdout=os.path.join(output_dir, "logs", "bam_to_stranded_bigwig.{sample}.stdout.txt"),
        stderr=os.path.join(output_dir, "logs", "bam_to_stranded_bigwig.{sample}.stderr.txt")

    container:
        "docker://docker.io/sambrycesmith/py_ranges_sam_bigwig:0.0.120_0.22.1_0.3.22"

    shell:
        """
        python {params.script} \
        -b {input.bam} \
        -r {input.bed} \
        -c {input.chromsizes} \
        -s {params.strandedness} \
        -o {params.output_prefix} \
        {params.keep_orphans} \
        1> {log.stdout} \
        2> {log.stderr}
        """


rule split_regions_by_strand:
    '''
    Split BED file by strand for generating coverage with external tools
    '''
    input:
        window_bed=os.path.join(output_dir, "regions", "pas_windows.{window}.bed")
    output:
        os.path.join(output_dir, "regions", "pas_windows.{window}.{strand}.bed")

    params:
        strand_key=lambda wildcards: strand_keys[wildcards.strand]

    log:
        stderr=os.path.join(output_dir, "logs", "split_regions_by_strands.{window}.{strand}.stderr.txt")

    shell:
        """
        awk -F"\\t" '{{if ($6=="{params.strand_key}") {{print $0}} }}' {input.window_bed} > {output} 2> {log.stderr}
        """

rule megadepth:
    '''
    '''
    input:
        bw=os.path.join(output_dir, "coverage", "{sample}.regions.{strand}.bw"),
        regions=rules.split_regions_by_strand.output

    output:
        os.path.join(output_dir, "coverage", "{sample}.megadepth.{window}.{strand}.annotation.tsv")
    
    params:
        output_prefix=os.path.join(output_dir, "coverage", "{sample}.megadepth.{window}.{strand}"),
        operation=config["megadepth_operation"]

    log:
        stderr=os.path.join(output_dir, "logs", "megadepth.{sample}.{window}.{strand}.stderr.txt")

    container:
        "docker://quay.io/biocontainers/megadepth:1.2.0--h43eeafb_6"

    shell:
        """
        megadepth {input.bw} \
        --annotation {input.regions} \
        --op {params.operation} \
        --prefix {params.output_prefix} \
        --no-annotation-stdout \
        2> {log.stderr}
        """


def aggregate_strands_coverage(wildcards):
    '''
    Return list of files produced by megadepth for each strand for given sample 
    '''

    # this get us the wo
    # checkpoint_output = checkpoints.split_regions_by_strand.get(strand=wildcards.strand).output[0]#

    return expand(os.path.join(output_dir, "coverage", "{sample}.megadepth.{window}.{strand}.annotation.tsv"),
    sample=wildcards.sample,
    window=wildcards.window,
    strand=strand_keys.keys()
    )


rule add_coverage_to_bed:
    '''
    for each sample 
    '''
    input:
        cov=aggregate_strands_coverage,
        window_bed=os.path.join(output_dir, "regions", "pas_windows.{window}.bed")

    output:
        os.path.join(output_dir, "coverage", "{sample}.pas_windows.{window}.summarised_coverage.bed")
    
    params:
        script="scripts/add_coverage_to_bed.py"

    log:
        stdout = os.path.join(output_dir, "logs", "add_coverage_to_bed.{sample}.{window}.stdout.txt"),
        stderr = os.path.join(output_dir, "logs", "add_coverage_to_bed.{sample}.{window}.tderr.txt")
    
    container:
        "docker://quay.io/biocontainers/pyranges:0.0.120--pyh7cba7a3_0"

    shell:
        """
        python {params.script} \
        -r {input.window_bed} \
        -o {output} \
        {input.cov} \
        1> {log.stdout} \
        2> {log.stderr}
        """
        



