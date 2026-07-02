rule sample_stats_table:
    input:
        config=workflow.source_path("report/datavzrd/sample_stats.yaml"),
        table="results/vcf/sample_stats.csv"
    output:
        report(
            directory("results/tables/sample_stats"),
            htmlindex="index.html",
            caption="report/sample_stats.rst",
            category="VCF",
            labels={"table": "Sample stats"}
        )
    log:
        "logs/report/sample_stats.log"
    wrapper:
        "v4.7.2/utils/datavzrd"