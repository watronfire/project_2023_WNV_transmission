rule preprocess_phylogeny:
    input:
        vcf = rules.convert_fasta_to_vcf.output.vcf,
        tree = rules.tree_building.output.tree
    output:
        protobuf = "generate_mat/coverage_masked.pb",
        log_dir = temp( directory( "generate_mat/usher_temp/" ) )
    threads: 8
    shell:
        """
        usher \
            --tree {input.tree} \
            --vcf {input.vcf} \
            --save-mutation-annotated-tree {output.protobuf} \
            --threads {threads} \
            --outdir {output.log_dir}
        """

#rule optimize_phylogeny:

rule generate_taxonium_input:
    input:
        protobuf = rules.preprocess_phylogeny.output.protobuf,
        genbank = "res/wnv.gb",
        reference = "res/wnv_reference.fasta",
        metadata = "res/wnv_metadata.csv"
    output:
        taxonium_input = "generate_mat/taxonium_input.jsonl.gz"
    params:
        columns = "strain,country,division,location,date,latitude,longitude,host"
    shell:
        """
        usher_to_taxonium \
            --input {input.protobuf} \
            --genbank {input.gtf_file} \
            --columns {params.columns} \
            --metadata {input.metadata} \
            --output {output.taxonium_input} 
        """
