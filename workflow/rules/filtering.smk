def mask_by_positions( input_loc, output_loc, positions ):
    from Bio import AlignIO
    from Bio.Seq import Seq
    records = AlignIO.read( input_loc,"fasta" )
    for record in records:
        seq_list = list( record.seq )
        for position in positions:
            seq_list[position - 1] = "N"
        record.seq = Seq( "".join( seq_list ) )

    AlignIO.write( records, output_loc,"fasta" )


rule replace_ambiguous:
    message: "Remove all ambiguous characters, and retain only ['A', 'C', 'T', 'G', 'N', '-']."
    input:
        alignment = rules.rename_fasta.output.renamed_alignment
    output:
        replaced_alignment = temp( "intermediates/ambiguous_masked_alignment.fasta" )
    run:
        from Bio import AlignIO
        from Bio.Seq import Seq

        def replace_inverse( seq, values, replacement ):
            new_seq = ""
            for char in seq:
                if char not in values:
                    new_seq += replacement
                else:
                    new_seq += char
            return Seq( new_seq )
        
        records = AlignIO.read( input.alignment, "fasta" )
        for record in records:
            record.seq = replace_inverse( record.seq, ['A', 'C', 'T', 'G', 'N'], 'N' )

            # While we're at it we might as well fix the names.
            record.id = record.id.replace( " " , "" )
            record.name = ""
            record.description = ""

        AlignIO.write( records, output.replaced_alignment, "fasta" )


rule length_filter:
    message: "Filters sequences that are less than full length and contain less than {params.max_ns} Ns"
    input:
        alignment = rules.replace_ambiguous.output.replaced_alignment
    params:
        max_ns = 0.3,
        values = ['A', 'C', 'T', 'G'] 
    output:
        length_filtered_fasta = temp( "intermediates/length_filtered.fasta" ),
        length_filtered_log = "intermediates/length_stats.csv"
    run:
        from Bio import AlignIO, Align
        import pandas as pd

        return_df = {"strain" : [], "length" : [], "N" : [], "percent_N" : [], "filtered" : [] }
        for value in params.values:
            return_df[value] = []

        records = AlignIO.read( input.alignment, "fasta" )
        passing_records = []
        for record in records:
            if record.id == 'Undetermined':
                continue
            return_df["strain"].append( record.id )
            
            sequence_length = len( record.seq )
            return_df["length"].append( len( record.seq ) )
            
            values_count = 0
            for value in params.values:
                count = record.seq.count( value )
                values_count += count
                return_df[value].append( count )
            return_df["N"].append( sequence_length - values_count )
            percent_ns = 1.0 - ( float( values_count ) / sequence_length )
            return_df["percent_N"].append( percent_ns )

            if percent_ns > params.max_ns:
                return_df["filtered"].append( True )
            else:
                return_df["filtered"].append( False )
                passing_records.append( record )

        return_df = pd.DataFrame( return_df )
        print( f"Filtered {return_df.loc[return_df['filtered']].shape[0]} genomes ({return_df.loc[return_df['filtered']].shape[0]/return_df.shape[0]:.1%} of all genomes)" )

        return_df.to_csv( output.length_filtered_log, index=False )
        AlignIO.write( Align.MultipleSeqAlignment( passing_records ), output.length_filtered_fasta, "fasta" )


rule mask_low_coverage_bases:
    message: "Mask bases that aren't observed in {params.min_coverage} proportion of genomes."
    input:
        alignment=rules.length_filter.output.length_filtered_fasta
    params:
        min_coverage=0.9
    output:
        masked_alignment="intermediates/coverage_masked_alignment.fasta",
        masked_positions="intermediates/coverage_masked_info.csv"
    run:
        from Bio import AlignIO
        from Bio.Seq import Seq
        import pandas as pd

        records = AlignIO.read( input.alignment,"fasta" )

        pos = list()
        first = True
        entries = len( records )
        for record in records:
            for idx, position in enumerate( str( record.seq ) ):
                if first:
                    pos.append( [1] if position != "N" else [0] )
                else:
                    pos[idx].append( 1 if position != "N" else 0 )
            first = False

        masked_positions = list()
        output_df = { "position": [], "coverage": [], "percent_coverage": [], "masked": [] }
        for idx, coverage in enumerate( pos ):
            output_df["position"].append( idx )
            output_df["coverage"].append( sum( coverage ) )
            output_df["percent_coverage"].append( sum( coverage ) / entries )
            output_df["masked"].append( (sum( coverage ) / entries) < params.min_coverage )
            if (sum( coverage ) / entries) < params.min_coverage:
                masked_positions.append( idx )

        for record in records:
            seq_list = list( record.seq )
            for position in masked_positions:
                seq_list[position] = "N"
            record.seq = Seq( "".join( seq_list ) )

        output_df = pd.DataFrame( output_df )
        output_df.to_csv( output.masked_positions,index=False )
        AlignIO.write( records,output.masked_alignment,"fasta" )


rule calculate_minimum_distance:
    input:
        alignment = rules.mask_low_coverage_bases.output.masked_alignment
    output:
        raw_distances = temp( "intermediates/distances.csv" ),
        min_distances = "intermediates/qc_metrics/distances.csv"
    threads: 8
    run:
        import pandas as pd
        import numpy as np

        shell( "pairsnp -ct 8 {input.alignment} > {output.raw_distances}" )

        dist = pd.read_csv( output.raw_distances, header=None, index_col=0 )
        ndist = dist.to_numpy()
        np.fill_diagonal( ndist,1000 )
        min_dist = pd.DataFrame( { "sequence": dist.index, "min_distance": ndist.min( axis=1 ) } )
        min_dist.to_csv( output.min_distances, index=False )


rule filter_by_qc_metrics:
    message: """ 
        Removes sequences from alignment if they meet any of the following:
            - iSNV count > {params.max_count}
            - Average iSNV frequency > {params.max_aaf}
            - Minimum pairwise distance > {params.min_distance}
        """
    input:
        variants = rules.combine_variants.output.variants,
        distances = rules.calculate_minimum_distance.output.min_distances,
        alignment = rules.mask_low_coverage_bases.output.masked_alignment
    params:
        max_count = 100,
        max_aaf = 0.25,
        min_distance = 32
    output:
        qc_report = "intermediates/qc_metrics/qc_report.csv",
        alignment = "intermediates/qc_filtered.fasta"
    shell:
        """
        python3 workflow/scripts/qc_filtering.py \
            --alignment {input.alignment} \
            --variants {input.variants} \
            --distances {input.distances} \
            --max-count {params.max_count} \
            --max-aaf {params.max_aaf} \
            --min-distance {params.min_distance} \
            --report {output.qc_report} \
            --output {output.alignment}
        """


rule tree_building:
    message: "Generate a tree using iqtree"
    input:
        alignment = rules.filter_by_qc_metrics.output.alignment
    output:
        tree = "intermediates/qc_filtered.fasta.treefile",
        report = temp( "intermediates/qc_filtered.fasta.iqtree" ),
        distances = temp( "intermediates/qc_filtered.fasta.mldist" ),
        log = temp( "intermediates/qc_filtered.fasta.log" )
    threads: 9
    shell:
        """
        iqtree -fast -m GTR+G -T 9 -af fasta -s {input.alignment} 
        """


rule ancestral_reconstruction:
    message: "Add mutations to tree. We'll use to infer homoplasies."
    input:
        alignment = rules.filter_by_qc_metrics.output.alignment,
        tree = rules.tree_building.output.tree
    output:
        outout_directory = directory( "intermediates/ancestral/" ),
        labeled_tree = "intermediates/ancestral/annotated_tree.nexus",
        infered_alignement = "intermediates/ancestral/ancestral_sequences.fasta"
    shell:
        """
        treetime ancestral \
            --aln {input.alignment} \
            --tree {input.tree} \
            --outdir {output.outout_directory} \
            --reconstruct-tip-states
        """


rule clean_inferred_alignment:
    message: "Removes NODE sequences from alignment. Not necessary for Delta, but enables easier viewing later"
    input:
        alignment = rules.ancestral_reconstruction.output.infered_alignement
    output:
        alignment = "intermediates/ancestral/cleaned_ancestral_sequences.fasta"
    run:
        from Bio import AlignIO, Align
        records = AlignIO.read( input.alignment, "fasta" )
        passing_records = []
        for record in records:
            if record.id.startswith( "NODE_" ):
                continue
            passing_records.append( record )

        AlignIO.write( Align.MultipleSeqAlignment( passing_records ), output.alignment, "fasta" )


rule calculate_homoplasies:
    input:
        tree = rules.ancestral_reconstruction.output.labeled_tree
    output:
        mutations = "intermediates/ancestral/mutations.csv"
    run:
        import pandas as pd
        from dendropy import Tree
        from collections import Counter

        tree = Tree.get( path=input.tree, schema="nexus", extract_comment_metadata=False )
        mutations = []
        for node in tree.preorder_node_iter():
            if len( node.comments ) > 0:
                comment = node.comments[0]
                mutations.extend( comment.split( '"' )[1].split( "," ) )
        df = {"mutation" : [], "position" : [], "occurance" : [], "kind" : [] }
        for k, v in Counter( mutations ).items():
            df["mutation"].append( k )
            df["occurance"].append( v )
            df["position"].append( "".join( [i for i in k if i.isnumeric()] ) )
            df["kind"].append( "all" )

        mutations = []
        for leaf in tree.leaf_node_iter():
            if len( leaf.comments ) > 0:
                comment = leaf.comments[0]
                mutations.extend( comment.split( '"' )[1].split( "," ) )
        for k, v in Counter( mutations ).items():
            df["mutation"].append( k )
            df["occurance"].append( v )
            df["position"].append( "".join( [i for i in k if i.isnumeric()] ) )
            df["kind"].append( "terminal" )

        df = pd.DataFrame( df )
        df = df.sort_values( ["kind","position"] )
        df.to_csv( output.mutations, index=False )


rule determine_initial_mask:
    input:
        mutations = rules.calculate_homoplasies.output.mutations
    params:
        min_mutations = 10
    output:
        positions = "intermediates/homoplasy_positions.txt"
    run:
        import pandas as pd
        mut = pd.read_csv( input.mutations )
        mut = mut.loc[~mut["mutation"].isna()&(mut["kind"] == "all")]
        mut["position"] = mut["position"].astype( int )
        mut = mut.sort_values( "position" )
        positions = mut.loc[mut["occurance"] > params.min_mutations, "position"].unique()

        with open( output.positions, "w" ) as output_file:
            output_file.write( "\n".join( map( str, positions ) ) )


rule mask_homoplasy_positions:
    input:
        alignment = rules.filter_by_qc_metrics.output.alignment,
        mask = rules.determine_initial_mask.output.positions,
    output:
        alignment = "intermediates/homoplasy_masking/homoplasy_mask_alignment.fasta"
    run:
        from Bio import AlignIO
        from Bio.Seq import Seq

        with open( input.mask, "r" ) as mask_file:
            mask = [int( line.strip() ) for line in mask_file]

        mask_by_positions( input.alignment, output.alignment, mask )


rule tree_building_after_mask:
    message: "Generate a tree using iqtree"
    input:
        alignment = rules.mask_homoplasy_positions.output.alignment
    output:
        tree = "intermediates/homoplasy_masking/homoplasy_mask_alignment.fasta.treefile",
        report = temp( "intermediates/homoplasy_masking/homoplasy_mask_alignment.fasta.iqtree" ),
        distances = temp( "intermediates/homoplasy_masking/homoplasy_mask_alignment.fasta.mldist" ),
        log = temp( "intermediates/homoplasy_masking/homoplasy_mask_alignment.fasta.log" )
    threads: 9
    shell:
        """
        iqtree -fast -m GTR+G -T 9 -af fasta -s {input.alignment} 
        """

rule reroot_resolve_tree:
    input:
        tree = rules.tree_building_after_mask.output.tree
    output:
        tree = "intermediates/delta_statistic/masked.tree"
    shell:
        """
        echo 'AF196835|USA|NewYork-Bronx|1999-09-16|nan|nan' | \
        gotree reroot outgroup -i {input.tree} -l - | \
        gotree resolve -o {output.tree}
        """

rule calculate_delta:
    input:
        tree = rules.reroot_resolve_tree.output.tree,
        alignment = rules.clean_inferred_alignment.output.alignment,
        positions = rules.determine_initial_mask.output.positions
    output:
        delta = "intermediates/delta.csv"
    shell:
        """
        Rscript workflow/scripts/calculate_delta.R \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --positions {input.positions} \
            --output {output.delta}
        """


rule mask_low_signal_homoplasic_positions:
    input:
        delta = rules.calculate_delta.output.delta,
        alignment = rules.filter_by_qc_metrics.output.alignment,
    params:
        minimum_signal = 3750
    output:
        alignment = "output/masked_alignment.fasta"
    run:
        delta = pd.read_csv( input.delta )
        masked_positions = delta.loc[delta["delta"]<params.minimum_signal,"position"].to_list()

        mask_by_positions( input_loc=input.alignment, output_loc=output.alignment, positions=masked_positions )


rule tree_building_final:
    input:
        alignment = rules.mask_low_signal_homoplasic_positions.output.alignment
    params:
        raw_tree = "output/masked_alignment.fasta.treefile"
    output:
        tree = "output/final.tree",
        report = temp("output/masked_alignment.fasta.iqtree"),
        distances = temp("output/masked_alignment.fasta.mldist"),
        log = temp("output/masked_alignment.fasta.log"),
        unique = temp("output/masked_alignment.fasta.uniqueseq.phy"),
        bionj = temp("output/masked_alignment.fasta.bionj"),
        checkpoint = temp("output/masked_alignment.fasta.ckp.gz")
    threads: 9
    shell:
        """
        iqtree -fast -redo -m GTR+G -T {threads} -af fasta -s {input.alignment} && \
        mv {params.raw_tree} {output.tree}
        """

rule calculate_homoplasy_distribution:
    message: "Identify mutation events that occur more times than expected given a normal poisson distribution."
    input:
        alignment = rules.filter_by_qc_metrics.output.alignment,
        tree = rules.tree_building.output.tree
    output:
        logfile = "intermediates/homoplasies.txt"
    shell:
        """
        treetime homoplasy \
            --aln {input.alignment} \
            --tree {input.tree} \
            --detailed \
            --verbose 6 > {output.logfile}
        """

rule run_ClonalFrameML:
    message: "Identify putative recombinant events as these are often indicative of contamination"
    input:
        tree = rules.tree_building.output.tree,
        alignment = rules.mask_low_coverage_bases.output.masked_alignment
    output:
        report = "intermediates/cfml_output/importation_status.txt",
        direct = directory( "intermediates/cfml_output/" )
    shell:
        """
        ClonalFrameML {input.tree} {input.alignment} {output.direct}
        """

rule remove_wholly_recombinant_sequences:
    input:
        report = rules.run_ClonalFrameML.output.report,
        alignment = rules.mask_low_coverage_bases.output.masked_alignment
    output:
        sequences_list = "intermediates/cfml_wholly_recombinant_sequences.txt",
        pruned_alignment = "intermediates/recombinant_pruned_alignment.fasta"
    run:
        import pandas as pd
        from Bio import AlignIO, Align

        cf = pd.read_csv( input.report ,sep="\t" )

        # Only remove sequences that are wholly recombinant. But we might increase to other.
        odd_rc = cf.loc[(cf["Beg"]==1)&(cf["End"]==10302),"Node"].to_list()
        with open( output.sequences_list, "w" ) as output_file:
            output_file.write( "\n".join( odd_rc ) )


        records = AlignIO.read( input.alignment,"fasta" )
        new_records = list()
        count = 0
        removed = 0
        for record in records:
            if record.id not in odd_rc:
                new_records.append( record )
                count += 1
            else: removed += 1
        print( f"Added {count} sequences to new alignment. Removed {removed} sequences." )
        AlignIO.write( Align.MultipleSeqAlignment( new_records ), output.pruned_alignment, "fasta" )

rule recombinant_masked_tree:
    message: "Generate a tree using the alignment in which wholly recombinant sequences have been removed"
    input:
        alignment = rules.remove_wholly_recombinant_sequences.output.pruned_alignment
    output:
        tree = "intermediates/recombinant_pruned_alignment.fasta.treefile"
    threads: 8
    shell:
        """
        iqtree -m GTR+G -T 8 -af fasta -s {input.alignment} 
        """