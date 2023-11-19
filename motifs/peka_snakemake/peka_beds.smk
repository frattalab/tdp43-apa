rule pas_to_peka_beds:
    input:
        pas=lambda wildcards: get_param(wildcards, 'regions_bed')
    output:
        comb_bed=os.path.join(base_outdir, "peka_beds", "{sample}.all.bed"),
        reg_bed=os.path.join(base_outdir, "peka_beds", "{sample}.foreground.bed"),
        bg_bed=os.path.join(base_outdir, "peka_beds", "{sample}.background.bed")
    params:
        script="scripts/pas_to_peka_beds.py",
        foreground_key=lambda wildcards: get_param(wildcards, 'foreground_key'),
        output_prefix = os.path.join(base_outdir, "peka_beds", "{sample}")

    log:
        stdout=os.path.join(base_outdir, "logs", "pas_to_peka_beds.{sample}.stdout.txt"),
        stderr=os.path.join(base_outdir, "logs", "pas_to_peka_beds.{sample}.stderr.txt")

    conda:
        "envs/peka_fork.yaml"

    shell:
       """python {params.script} \
       -f {params.foreground_key} \
       {input.pas} \
       {params.output_prefix} \
       1> {log.stdout} \
       2> {log.stderr}
       """
       