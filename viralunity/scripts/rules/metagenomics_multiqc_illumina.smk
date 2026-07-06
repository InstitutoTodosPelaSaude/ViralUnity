rule generate_multiqc_report:
    input:
        expand(
            config["output"] + "qc/reports/trim.{sample}_fastp.json",
            sample=config["samples"]
        )
    output:
        config["output"] + "qc/reports/multiqc_report.html"
    params:
        temp = config["output"]
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        multiqc -f --cl-config "extra_fn_clean_exts: ['_R1']" -o {params.temp}/qc/reports/ {params.temp}/qc/reports/
        """
