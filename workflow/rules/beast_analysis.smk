rule calculate_markov_jumps:
    input:
        jump_history = "beast_runs/discrete/2023-03-24_wnv_skygrid_yearly_relaxed_discrete.Location.history.trees",
        xml = rules.generate_posthoc_xml.output.xml
    params:
        burnin = int( 55000 / 500 ),
        mrsd = 2020.5,
    output:
        jumps = "beast_runs/discrete/2023-04-03_wnv_skygrid_yearly_relaxed_discrete.Location.csv"
    shell:
        """
        MRSD=$(fgrep date {input.xml} | cut -f2 -d\\" | sort | tail -n 2 | head -n 1) && \
        java -Xmx16G -cp ~/scripts/BEASTv1.10.5pre_thorney_0.1.2/lib/beast.jar dr.app.tools.TreeMarkovJumpHistoryAnalyzer \
            -burnin {params.burnin} \
            -mrsd $MRSD \
            {input.jump_history} \
            {output.jumps}
        """


# TODO: rule for making MCC tree
rule plot_mcc_tree:
    input:
        tree = "beast_runs/discrete/2023-03-24_wnv_skygrid_yearly_relaxed_discrete.mcc.tree"
    params:
        us_map = "res/shapefiles/cb_2018_us_state_20m.shp"
    log:
        notebook = "analyses/figureX_mcc-tree.ipynb"
    output:
        tree_figure = "analyses/plots/figureX_mcc-tree.pdf",
        tree_legend = "analyses/plots/figureX_mcc-tree-legend.pdf"
    notebook: "../notebooks/FigureX_mcc-tree_template.py.ipynb"


rule plot_transition_rates:
    input:
        log = "beast_runs/discrete/2023-03-24_wnv_skygrid_yearly_relaxed_discrete.log",
        traits = rules.calculate_discrete_states.output.traits
    params:
        us_map = "res/shapefiles/cb_2018_us_state_20m.shp",
        burnin = 55000,
        minimum_BF = 100
    log:
        notebook="analyses/figureX_transition-rate-geography.ipynb"
    output:
        map_figure = "analyses/plots/figureX_transition-rate-map.pdf",
        cc_figure= "analyses/plots/figureX_transition-rate-closeness-centrality.pdf",
        rates_figure= "analyses/plots/figureX_tranition-rate-totals.pdf"
    notebook: "../notebooks/figureX_transition-rate-geography_template.py.ipynb"


rule identify_transmission_lineages:
    input:
        beast_trees = "beast_runs/discrete/2023-04-03_wnv_skygrid_yearly_relaxed_discrete.down.trees",
        xml = rules.generate_posthoc_xml.output.xml
    params:
        burnin = int( 55000 / 500 ),
        trait = "Location"
    output:
        transmission_lineages = "data/whole_transmission_lineages.csv"
    shell:
        """
        states=(`fgrep "state code" {input.xml} | cut -f2 -d\\"`) && \
        python workflow/scripts/identify_lineages.py \
            --trees {input.beast_trees} \
            --burnin {params.burnin} \
            --states ${{states[*]}} \
            --trait {params.trait} \
            --output {output}
        """


rule plot_transition_lineage_metrics:
    input:
        transmission_lineages = rules.identify_transmission_lineages.output.transmission_lineages
    log:
        notebook="analyses/figureX_transmission-lineages.ipynb"
    output:
        tl_count_figure = "analyses/plots/figureX_transmission-lineage-count.pdf",
        tl_size_figure = "analyses/plots/figureX_transmission-lineage-size.pdf",
        tl_length_figure = "analyses/plots/figureX_transmission-lineage-length.pdf",
        tl_post_count_figure= "analyses/plots/figureX_transmission-lineage-count-post2002.pdf",
        tl_post_size_figure= "analyses/plots/figureX_transmission-lineage-size-post2002.pdf",
        tl_post_length_figure= "analyses/plots/figureX_transmission-lineage-length-post2002.pdf"
    notebook: "../notebooks/figureX_transmission-lineage-metrics_template.py.ipynb"