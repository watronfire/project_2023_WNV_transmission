rule merge_metadata:
    input:
        wn4k_md = "input/WNV_California_Samples.csv",
        other_md = "input/2017.11.07_WNV_download.csv"
    output:
        metadata = "input/merged_metadata.csv"
    run:
        import pandas as pd
        import numpy as np
        import unicodedata

        states = pd.read_csv( "https://raw.githubusercontent.com/jasonong/List-of-US-States/master/states.csv" )
        states = states.set_index( "Abbreviation" )["State"].to_dict()

        usecols = ["Scripps_ID", "Collection date", "Species", "Country", "State", "County", "Latitude", "Longitude"]
        md = pd.read_csv( input.wn4k_md, parse_dates=["Collection date"], usecols=usecols )
        md = md.loc[~md["Collection date"].isna()]
        md = md.rename( columns={ "Scripps_ID": "accession", "Species": "host", "Collection date": "collection_date" } )
        md.columns = map( lambda x: x.lower(), md.columns )
        md["state"] = md["state"].replace( states )
        md["scripps_num"] = md["accession"].apply( lambda x: int( x[1:] ) )

        usecols = ["accession", "host", "country", "state", "county", "collection_date", "latitude", "longitude"]
        od = pd.read_csv( input.other_md, parse_dates=["collection_date"], usecols=usecols )

        # Fix dates that were only given as a year.
        od.loc[(od["collection_date"].dt.month == 1) & (
                    od["collection_date"].dt.day == 1), "collection_date"] += pd.DateOffset( months=5 )
        od = od.loc[~od["collection_date"].isna()]

        total = pd.concat( [md, od] )
        total["country"] = total["country"].fillna( "USA" )
        total["latitude"] = total["latitude"].replace( { "x": np.nan, "(not provided)": np.nan } )
        total["longitude"] = total["longitude"].replace( { "x": np.nan, "(not provided)": np.nan } )
        total["latitude"] = total["latitude"].astype( float )
        total["longitude"] = total["longitude"].astype( float )
        total["latitude"] = total["latitude"].round( 4 )
        total["longitude"] = total["longitude"].round( 4 )
        total["county"] = total["county"].str.replace( "_", " " ).str.replace( ".", "" ).str.replace( "'", "" )

        total["location_str"] = "-" + total["county"].astype( str )
        total.loc[total["location_str"] == "-nan", "location_str"] = ""
        total["location_str"] = total["state"].astype( str ) + total["location_str"].astype( str )
        total["long_name"] = total["accession"] + "|" + total["country"] + "|" + total["location_str"] + "|" + total[
            "collection_date"].dt.strftime( "%Y-%m-%d" ) + "|" + total["latitude"].astype( str ) + "|" + total[
                                 "longitude"].astype( str )
        total["long_name"] = total["long_name"].str.replace( " ", "" ).str.replace( "'", "" )
        total["long_name"] = total[
            "long_name"].apply( lambda val: unicodedata.normalize( 'NFKD', val ).encode( 'ascii', 'ignore' ).decode() )
        total.drop( columns=["location_str"] )

        total.to_csv( output.metadata, index=False )


rule rename_fasta:
    input:
        alignment = "input/2023.04.12_alignment.fasta",
        metadata = rules.merge_metadata.output.metadata
    output:
        renamed_alignment = "input/formatted_alignment.fasta"
    run:
        import pandas as pd
        from Bio import AlignIO, Align

        md = pd.read_csv( input.metadata )
        accession_dict = md.set_index( "accession" )["long_name"].to_dict()
        num_dict = md.loc[~md["scripps_num"].isna()].set_index( "scripps_num" )["long_name"].to_dict()


        records = AlignIO.read( input.alignment, "fasta" )
        renamed_records = []
        for record in records:
            record_accession = record.id.split( "_" )[0]
            if record_accession in accession_dict:
                record.id = accession_dict[record_accession]
                record.name = ""
                record.description = ""
            elif record_accession.startswith( "W" ):
                num = int( record_accession[1:] )
                if num in num_dict:
                    record.id = num_dict[num]
                    record.name = ""
                    record.description = ""
                else:
                    print( f"{record.id} begins with W and is not found in metadata. Removing..." )
                    continue
            else:
                print( f"{record.id} is not found in metadata. Removing..." )
                continue
            renamed_records.append( record )

        AlignIO.write( Align.MultipleSeqAlignment( renamed_records ), output.renamed_alignment, "fasta" )


rule calculate_variants:
    input:
        reference = "res/WNV_REF_COAV997.fasta",
        bed_file = "res/WNV_400.bed"
    output:
        bam_file = temp( "intermediates/bams/{sample}.bam" ),
        variants = temp( "intermediates/variants/{sample}.tsv" )
    params:
        minimum_quality = 20,
        minimum_frequency = 0.03,
        minimum_coverage = 5,
        primer_offset = 3,
        prefix = "intermediates/variants/{sample}",
        cloud_file = lambda wildcards: CLOUD_FILES[wildcards.sample]
    shell:
        """
        gsutil cp {params.cloud_file} {output.bam_file} &&
        ivar trim \
            -i {output.bam_file} \
            -b {input.bed_file} \
            -x {params.primer_offset} | \
        samtools sort - | \
        samtools mpileup \
            -aa -A -d 0 -Q 0 \
            --reference {input.reference} \
            - | \
        ivar variants \
            -p {params.prefix} \
            -q {params.minimum_quality} \
            -t {params.minimum_frequency} \
            -m {params.minimum_coverage} \
            -r {input.reference}
        """


rule combine_variants:
    input:
        variants = expand( "intermediates/variants/{sample}.tsv",sample=CLOUD_FILES )
    output:
        variants = "intermediates/qc_metrics/variants.csv"
    run:
        import pandas as pd
        import os

        usecols = [
            "POS", "REF", "ALT",
            "REF_DP", "REF_RV", "REF_QUAL",
            "ALT_DP", "ALT_RV", "ALT_QUAL", 'ALT_FREQ',
            'TOTAL_DP', 'PVAL', 'PASS'
        ]

        return_df = list()
        for file in input.variants:
            name = os.path.basename( file ).replace( ".tsv", "" )
            df = pd.read_csv( file, sep="\t", usecols=usecols )
            df["SAMPLE"] = name
            return_df.append( df )

        return_df = pd.concat( return_df )
        return_df.to_csv( output.variants, index=False )