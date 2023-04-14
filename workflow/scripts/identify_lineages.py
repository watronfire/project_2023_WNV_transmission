import argparse
from dendropy import Tree
import time
from datetime import datetime as dt
import pandas as pd


def _toYearFraction( date, date_format="%Y-%m-%d" ):
    """ Converts datetime object to a decimal year
    Parameters
    ----------
    date: str
        date to be converted.

    Returns
    -------
    float
        date in decimal year format.
    """

    def sinceEpoch( d ):  # returns seconds since epoch
        return time.mktime( d.timetuple() )

    date = dt.strptime( date, date_format )

    year = date.year
    start_of_this_year = dt( year=year, month=1, day=1 )
    start_of_next_year = dt( year=year + 1, month=1, day=1 )

    year_elapsed = sinceEpoch( date ) - sinceEpoch( start_of_this_year )
    year_duration = sinceEpoch( start_of_next_year ) - sinceEpoch( start_of_this_year )
    fraction = year_elapsed / year_duration

    return date.year + fraction


def get_children( node, field ):
    return_list = []
    location = node.annotations.get_value( field )
    for j in node.child_node_iter( lambda x: x.annotations.get_value( field ) == location ):
        if j.is_leaf():
            return_list.append( j )
        else:
            return_list.extend( get_children( j, field ) )
    return return_list


def get_ancestral_state( node, field ):
    return next( node.ancestor_iter() ).annotations.get_value( field )


def get_introductions( tree, field, interest, root_height=0 ):
    return_list = list()

    def _identity_func( node ):
        if node.parent_node is None:
            return False
        return (node.annotations.get_value( field ) == interest) & (
                get_ancestral_state( node, field ) != interest)

    for i in tree.preorder_node_iter( filter_fn=_identity_func ):
        node_dict = i.annotations.values_as_dict()
        node_dict[f"par.{field}"] = get_ancestral_state( i, field )

        children = get_children( i, field )
        children_dates = [root_height + float( child.root_distance ) for child in children]
        node_dict["children"] = len( children )
        node_dict["children.names"] = [child.taxon.label for child in children]
        if len( children ) > 0:
            node_dict["earliest_date"] = min( children_dates )
            node_dict["latest_date"] = max( children_dates )
        node_dict["children.dates"] = children_dates
        node_dict["height"] = root_height + i.root_distance
        return_list.append( node_dict )
    return pd.DataFrame( return_list )


def parse_lineages( trees, burnin, trait, states, output ):
    yielder = Tree.yield_from_files( files=[trees], schema="nexus", preserve_underscores=True )
    recent_date = 0
    output_df = list()
    start_time = time.time()
    for tree_idx, t in enumerate( yielder ):
        if tree_idx == 0:
            recent_date = _toYearFraction( max( [i.label.split( "|" )[3] for i in t.taxon_namespace] ) )
        if tree_idx < burnin:
            continue

        print( f"Processing {tree_idx - burnin} tree...", end="" )

        root_height = recent_date - t.max_distance_from_root()
        t.calc_node_root_distances()

        for loc in states:
            temp_df = get_introductions( t, trait, loc, root_height=root_height )
            temp_df["tree"] = tree_idx
            output_df.append( temp_df )

        execution_time = time.time() - start_time
        print( f" Done in {execution_time:.1f} seconds" )
        start_time = time.time()

    output_df = pd.concat( output_df )
    output_df.to_csv( output, index=False )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="generates a list of transmission lineages from a posterior distribution of trees." )

    # Initialize optional arguments
    parser.add_argument( "--trees", help="path to trees file. Must be in Nexus format." )
    parser.add_argument( "--burnin", type=int, help="number of trees to discard as part of burnin" )
    parser.add_argument( "--states", nargs="+", help="location states to identify transmission lineages for" )
    parser.add_argument( "--trait", help="field in annotations representatin discrete location state" )
    parser.add_argument( "--output", help="path to save output CSV file." )

    args = parser.parse_args()

    print( f"Parsing lineages for states: {args.states}")
    parse_lineages( args.trees, args.burnin, args.trait, args.states, args.output )
