rule generate_bootstrapped_tree:
    input:
        alignment = rules.mask_homoplasy_positions.output.alignment
    params:
        output_dir = "intermediates/selection/bootstrapped_tree",
        bootstraps = 100
    output:
        tree="intermediates/selection/bootstrapped_tree.treefile",
        report=temp( "intermediates/selection/bootstrapped_tree.iqtree" ),
        distances=temp( "intermediates/selection/bootstrapped_tree.mldist" ),
        log=temp( "intermediates/selection/bootstrapped_tree.log" ),
        consensus_tree = temp( "intermediates/selection/bootstrapped_tree.contree"),
        supports = temp( "intermediates/selection/bootstrapped_tree.splits.nex"),
        splits = temp( "intermediates/selection/bootstrapped_tree.splits")
    threads: 16
    shell:
        """
        iqtree -m GTR+G -T AUTO -af fasta \
            -s {input.alignment} \
            --prefix {params.output_dir} \
            -B {params.bootstraps}
        """

rule prepare_alignment:
    input:
        alignment = rules.mask_homoplasy_positions.output.alignment
    output:
        alignment = "intermediates/selection/masked.codon_removed.fasta"
    run:
        from Bio import AlignIO
        from Bio.Seq import Seq

        records = AlignIO.read( input.alignment, "fasta" )

        output_list = list()
        for record in records:
            seq_list = list( record.seq )
            for i in range( 3, len( seq_list ) + 1, 3 ):
                substring = seq_list[i-3:i]
                # Replace stop codons with ambiguous nucleotides if found.
                if substring in [["T", "A", "G"], ["T", "A", "A"], ["T", "G", "A"]]:
                    seq_list[i-1] = "N"
                    seq_list[i-2] = "N"
                    seq_list[i-3] = "N"
            record.seq = Seq( "".join( seq_list ) )
        AlignIO.write( records, output.alignment, "fasta" )

rule hyphy_wrapper:
    input:
        alignment = rules.prepare_alignment.output.alignment,
        tree = rules.generate_bootstrapped_tree.output.tree
    params:
        pvalue = 0.05,
        type = lambda wildcards: wildcards.type
    output:
        results = "intermediates/selection/homoplasy_masked.{type}.json"
    log: "intermediates/selection/homoplasy_masked.{type}.log"
    shell:
        """
        /gpfs/home/natem/scripts/hyphy/hyphy LIBPATH=/gpfs/home/natem/scripts/hyphy/res/ {params.type} \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --pvalue {params.pvalue} \
            --output {output.results} 2>&1 | tee {log}
        """