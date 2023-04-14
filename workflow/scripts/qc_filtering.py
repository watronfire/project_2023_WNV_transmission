import argparse
from subprocess import run

import pandas as pd
import numpy as np
import scipy.stats as st
from Bio import AlignIO, Align


def fit_expon( entry, func=st.expon ):
    fit = func.fit(entry["ALT_FREQ"])
    results = st.kstest(entry["ALT_FREQ"], func.cdf, fit)
    return results.pvalue > 0.01

def weighted_mean( entry ):
    return np.average( entry["ALT_FREQ"], weights=entry["TOTAL_DP"] )

def weighted_std( entry ):
    average = weighted_mean( entry )
    variance = np.average( ( entry["ALT_FREQ"] - average )**2, weights=entry["TOTAL_DP"] )
    return np.sqrt( variance )

def load_summarize_variants( variants_loc ):
    vr = pd.read_csv( variants_loc )
    vr = vr.loc[vr["PASS"]]
    vr = vr.loc[~vr["ALT"].str.contains( "\+|-" )]

    # Correct major and minor variants
    maj = vr.loc[vr["ALT_FREQ"] > 0.5]
    mino = vr.loc[vr["ALT_FREQ"] <= 0.5]

    rename_dict = {
        "REF": "ALT",
        "ALT": "REF",
        "REF_DP": "ALT_DP",
        "REF_RV": "ALT_RV",
        "REF_QUAL": "ALT_QUAL",
        "ALT_DP": "REF_DP",
        "ALT_RV": "REF_RV",
        "ALT_QUAL": "REF_QUAL",
    }

    maj = maj.rename( columns=rename_dict )
    maj["ALT_FREQ"] = 1 - maj["ALT_FREQ"]
    maj = maj.loc[(maj["ALT_FREQ"] > 0.03) & (maj["ALT_DP"] > 5)]

    vr = pd.concat( [mino, maj] )

    return  vr.groupby( "SAMPLE" ).apply( lambda r: pd.Series( {"mean": weighted_mean(r), "std": weighted_std(r), "count" : r.shape[0], "fit_expon" : fit_expon(r)} ) )


def load_distances( distances_loc ):
    return pd.read_csv( distances_loc )


def initialize_report( alignment ):
    alignment_names = run( f"fgrep '>' {alignment}", shell=True, text=True, capture_output=True  )
    report = { "name" : [], "strain" : [] }
    for name in alignment_names.stdout.split( "\n" ):
        report["name"].append( name[1:] )
        report["strain"].append( name[1:].split( "|" )[0] )
    report = pd.DataFrame( report )
    report = report.loc[report["strain"]!=""]
    return report


def qc_filter( alignment, variants_loc, distances_loc, max_count, max_aaf, min_distance, report_loc, output_loc ):
    summary = load_summarize_variants( variants_loc )
    distances = load_distances( distances_loc )

    report = initialize_report( alignment )
    report = report.merge( summary, left_on="strain", right_on="SAMPLE", how="left" )
    report = report.merge( distances, left_on="name", right_on="sequence", how="left" )

    report["drop_count"] = report["count"] > max_count
    report["drop_aaf"] = report["mean"] > max_aaf
    report["drop_distance"] = report["min_distance"] > min_distance
    report.to_csv( report_loc, index=False )

    # TODO: Given the low sampling, this should be replaced with a clock filter. We have one in place but its later in pipeline.
    #drops = report.loc[report[["drop_count", "drop_aaf", "drop_distance"]].any(axis=1),"name"]
    drops = report.loc[report[["drop_count", "drop_aaf"]].any(axis=1),"name"]

    records = AlignIO.read( alignment, "fasta" )
    passing_records = []
    for record in records:
        if record.id in drops:
            continue
        passing_records.append( record )

    AlignIO.write( Align.MultipleSeqAlignment( passing_records ), output_loc, "fasta" )


if __name__ == "__main__":
    parser = argparse.ArgumentParser( description="" )

    # Initialize optional arguments
    parser.add_argument( "--alignment" )
    parser.add_argument( "--variants" )
    parser.add_argument( "--distances" )
    parser.add_argument( "--max-count", type=int )
    parser.add_argument( "--max-aaf", type=float )
    parser.add_argument( "--min-distance", type=int )
    parser.add_argument( "--report" )
    parser.add_argument( "--output")

    args = parser.parse_args()

    qc_filter(
        args.alignment,
        args.variants,
        args.distances,
        args.max_count,
        args.max_aaf,
        args.min_distance,
        args.report,
        args.output
    )
