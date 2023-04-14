rule get_dates:
    input:
        seqs = rules.mask_low_signal_homoplasic_positions.output.alignment
    output:
        dates = "intermediates/dates.csv"
    run:
        from datetime import datetime

        def date_to_dec( date_str ):
            date_obj = datetime( *map( int, date_str.split( "-" ) ) )
            return ( float( date_obj.strftime( "%j" ) ) - 1 ) / 365.25 + float( date_obj.strftime( "%Y" ) )

        with open( input.seqs, "r" ) as selected, open( output.dates, "w" ) as dates:
            dates.write( "name,date\n" )
            for line in selected:
                if line.startswith( ">" ):
                    name = line.strip()[1:]
                    date = name.split( "|" )[3]
                    dates.write( f"{name},{date_to_dec(date)}\n" )


rule timetree:
    input:
        selected_tree = rules.tree_building_final.output.tree,
        selected_alignment = rules.mask_low_signal_homoplasic_positions.output.alignment,
        dates = rules.get_dates.output.dates
    params:
        outdir = "intermediates/timetree_final/"
    output:
        shell_output = "intermediates/timetree_final/timetree.log",
        timetree = "intermediates/timetree_final/timetree.nexus",
        divergence_tree = "intermediates/timetree_final/divergence_tree.nexus"
    shell:
        """
        treetime \
            --tree {input.selected_tree} \
            --aln {input.selected_alignment} \
            --dates {input.dates} \
            --reroot 'AF196835|USA|NewYork-Bronx|1999-09-16|nan|nan' \
            --clock-rate 0.0005 \
            --outdir {params.outdir} > {output.shell_output}
        """


rule find_outliers_in_timetree:
    input:
        divergence_tree = rules.timetree.output.divergence_tree,
        dates = rules.get_dates.output.dates
    output:
        to_prune = "intermediates/timetree_final/prunes.txt"
    run:
        import pandas as pd
        from dendropy import Tree
        from scipy.stats import linregress

        tree = Tree.get( path=input.divergence_tree, schema="nexus" )
        distances = dict()
        for node in tree.leaf_node_iter():
            distances[node.taxon.label] = node.distance_from_root()

        dates = pd.read_csv( input.dates )
        dates["distance"] = dates['name'].map( distances )
        dates = dates.loc[~dates["distance"].isna()]

        res = linregress( dates["date"], dates["distance"] )
        dates["expected"] = res.intercept + res.slope * dates["date"]
        dates["residual"] = dates["distance"] - dates["expected"]
        incorrect_taxa = dates.loc[(dates["residual"]<-dates["residual"].mad() * 3)|(dates["residual"]>dates["residual"].mad() * 3), "name"].to_list()

        with open( output.to_prune, "w" ) as prunes:
            prunes.write( "taxon\n" )
            prunes.write( "\n".join( incorrect_taxa ) )


rule clock_filter_tree:
    input:
        selected_tree = rules.timetree.output.timetree,
        clock_violators = rules.find_outliers_in_timetree.output.to_prune
    output:
        pruned_tree_nwk = "intermediates/timetree_final/timetree_pruned.tree"
    shell:
        """
        jclusterfunk prune \
            -i {input.selected_tree} \
            -f newick \
            --taxon-file {input.clock_violators} \
            -o {output.pruned_tree_nwk} \
            -v
        """


rule resolve_tree:
    input:
        tree = rules.clock_filter_tree.output.pruned_tree_nwk
    output:
        resolved_tree = "intermediates/timetree_final/resolved.tree"
    shell:
        """
        gotree resolve --input {input.tree} --output {output.resolved_tree}
        """


rule generate_nexus:
    input:
        resolved_tree = rules.resolve_tree.output.resolved_tree,
        alignment = rules.mask_low_signal_homoplasic_positions.output.alignment
    output:
        nexus_file = "output/final.nexus"
    shell:
        """
        python3 workflow/scripts/generate_nexus.py -a {input.alignment} -t {input.resolved_tree} -o {output.nexus_file}
        """

rule generate_xml:
    input:
        nexus = rules.generate_nexus.output.nexus_file,
    params:
        template = "res/wnv.template",
        properties = 'chain_length=100000000,log_screen=100000,log_every=10000,stem=2023-04-12_wnv_skygrid_yearly_relaxed'
    output:
        xml = "output/2023-04-12_wnv_skygrid_yearly_relaxed.xml"
    shell:
        """
        python workflow/scripts/generate_beast_xml.py \
            --template {params.template} \
            --nexus {input.nexus} \
            -D {params.properties:q} \
            --order 3 \
            --output {output.xml}
        """

rule calculate_discrete_states:
    input:
        nexus = rules.generate_nexus.output.nexus_file
    params:
        additional = "res/additional_states.csv",
        min_sequences = 10
    output:
        traits = "intermediates/final_traits.csv"
    run:
        from dendropy import Tree
        from collections import Counter

        regions = {
            "NorthEast": ["Connecticut", "Massachusetts", "RhodeIsland", "NewJersey", "NewYork", "Pennsylvania",
                          "NewHampshire"],
            "South": ["DistrictofColumbia", "Florida", "Louisiana", "Georgia", "Maryland", "NorthCarolina", "Virginia",
                      "Alabama", "Kentucky", "Mississippi", "Tennessee", "Arkansas", "Oklahoma", "Texas",
                      "VirginIslands"],
            "Midwest": ["Illinois", "Michigan", "Ohio", "Wisconsin", "Iowa", "Kansas", "Minnesota", "Missouri",
                        "Nebraska", "Indiana", "NorthDakota", "SouthDakota"],
            "West": ["Arizona", "Idaho", "Colorado", "Montana", "Nevada", "NewMexico", "California", "Oregon",
                     "Washington", "Utah"]
        }


        def get_region( s, r ):
            for k, v in r.items():
                if s in v:
                    return k
            return None

        with open( params.additional, "r" ) as add:
            additional_states = {ls.split( "," )[0] : ls.split( "," )[1].strip() for ls in add}

        t = Tree.get( path=input.nexus, schema="nexus" )
        states = []
        for i in t.taxon_namespace:
            country = i.label.split( "|" )[1]
            state = i.label.split( "|" )[2].split( "-" )[0]
            if (country == "USA") & (state != "nan"):
                states.append( state )

        counts = Counter( states )
        accepted_states = [state for state in counts if counts[state]>=params.min_sequences]

        with open( output.traits, "w" ) as traits:
            traits.write( "taxon\tstate\n" )
            for i in t.taxon_namespace:
                state = i.label.split( "|" )[2].split( "-" )[0]
                if state in accepted_states:
                    write_state = state
                elif i.label.split( "|" )[0] in additional_states:
                    id = i.label.split( "|" )[0]
                    if additional_states[id] in accepted_states:
                        write_state = additional_states[id]
                    else:
                        write_state = get_region( additional_states[id], regions )
                else:
                    region = get_region( state, regions )
                    if region:
                        write_state = region
                    else:
                        write_state = "Other"

                traits.write( f"{i.label}\t{write_state}\n" )


rule generate_posthoc_xml:
    input:
        nexus = rules.generate_nexus.output.nexus_file,
        traits = rules.calculate_discrete_states.output.traits
    params:
        template = "res/wnv_discrete.template",
        properties = 'chain_length=2000000,log_screen=500,log_every=500,stem=2023-04-03_wnv_skygrid_yearly_relaxed_discrete'
    output:
        xml = "output/2023-04-03_wnv_skygrid_yearly_relaxed_discrete.xml"
    shell:
        """
        python workflow/scripts/generate_beast_xml.py \
            --template {params.template} \
            --nexus {input.nexus} \
            --traits {input.traits} \
            -D {params.properties:q} \
            --order 3 \
            --output {output.xml}
        """


### Note
### The correct run to start a beast docker is `docker run -it --gpus all --entrypoint /bin/bash -v <local>:<docker_dir>  andersenlabapps/beast-beagle`
### Once inside, use the command `/root/beast-beagle-docker/run <xml> <run_number> <gpu>`
