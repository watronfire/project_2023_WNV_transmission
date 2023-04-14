import argparse
from jinja2 import Environment, FileSystemLoader, select_autoescape
from dendropy import DataSet
from datetime import datetime as dt
import time
import numpy as np

def _toYearFraction( date ):
    """ Converts datetime object to a decimal year
    Parameters
    ----------
    date: datetime.datetime
        date to be converted.

    Returns
    -------
    float
        date in decimal year format.
    """
    def sinceEpoch( d ): # returns seconds since epoch
        return time.mktime( d.timetuple() )
    s = sinceEpoch

    year = date.year
    start_of_this_year = dt( year=year, month=1, day=1 )
    start_of_next_year = dt( year=year+1, month=1, day=1 )

    year_elapsed = s( date ) - s( start_of_this_year )
    year_duration = s( start_of_next_year ) - s( start_of_this_year )
    fraction = year_elapsed / year_duration

    return date.year + fraction

def parse_properties( properties ):
    """ Parses the properties argument of program. No checking or cleaning of input is performed so tread lightly.
    Parameters
    ----------
    properties: str
        Properties input of program. Should be a comma seperate list of 'key=value' pairs.

    Returns
    -------
    dict
        Dictionary containing parsed 'key=value' entries.

    """
    return_dict = dict()

    if properties is None:
        return return_dict
    else:
        for pair in properties.strip( "'" ).split( "," ):
            pair_split = pair.split( "=" )
            return_dict[pair_split[0].strip()] = pair_split[1].strip()

    return return_dict


def parse_date( label, sep="/", order=-1, date_format="%Y-%m-%d" ):
    """ Parses a date from a taxa label.
    Parameters
    ----------
    label: str
        taxa label containing a date
    sep: char
        character to split label by
    order: int
        index of item from split label that corresponds to a date
    date_format: str
        the format of the date string split from taxa label

    Returns
    -------
    float
        date parsed from taxa label in decimal year format
    """
    date_str = label.split( sep )[order]
    date = dt.strptime( date_str, date_format )
    return _toYearFraction( date )

def parse_nexus( nexus_loc, sep="|", order=3 ):
    """ Parses nexus file into a dictionary for template engine
    Parameters
    ----------
    nexus_loc: str
        location of nexus file
    sep: char
        deliminator in taxa name
    order: int
        index of date in taxa name after splitting by sep

    Returns
    -------
    dict
        dictionary containing tree in newick format, and a list of taxa containing their sequences, dates, and labels.
    """
    return_dict = { "taxa" : {} }

    print( "Parsing nexus file... ", end="" )
    nexus = DataSet.get( path=nexus_loc, schema="nexus" )

    return_dict["tree"] = nexus.tree_lists[0].as_string( "newick" )

    count = 0

    max_date = 2020.1
    min_date = 2020.1
    for k, v in nexus.char_matrices[0].items():
        count += 1

        date = parse_date( k.label, sep=sep, order=order )
        if date > max_date:
            max_date = date
        if date < min_date:
            min_date = date

        taxa_dict = { "id" : k.label,
                      "date" : date,
                      "sequence" : str( v ) }
        return_dict["taxa"][k.label] = taxa_dict
    print( "Done. {} taxa loaded".format( count ) )

    return_dict["site_count"] = nexus.char_matrices[0].max_sequence_size
    return_dict["taxon_count"] = len( return_dict["taxa"] )
    return_dict["cutoff"] = max_date - 2019.75
    return_dict["grid_points"] = round( return_dict["cutoff"] * 52 )

    return_dict["min_age"] = max_date - min_date

    return return_dict


def add_traits( values_dict, traits_loc, sep="\t", shuffle=False ):

    print( "Parsing traits file... ", end="" )
    count = 0

    taxa = []
    traits = []

    with open( traits_loc, "r" ) as traits_file:
        header = next( traits_file )
        value = header.split( sep )[1].strip()
        for line in traits_file:
            count += 1
            line_split = line.strip().split( sep )
            taxa.append( line_split[0] )
            traits.append( line_split[1] )

    if shuffle:
        for taxon, trait in zip( taxa, np.random.choice( traits, len( traits ), replace=False ) ):
            values_dict["taxa"][taxon][value] = trait
    else:
        for taxon, trait in zip( taxa, traits ):
            values_dict["taxa"][taxon][value] = trait
    trait_list = sorted( list( set( traits ) ) )    # Make the order deterministic
    traits_dict = dict()
    for idx, trait in enumerate( trait_list ):
        reward_array = [0.0] * len( trait_list )
        reward_array[idx] = 1.0
        traits_dict[trait] = {"id" : trait, "reward" : reward_array}

    values_dict["states"] = traits_dict
    values_dict["state_rates"] = ( ( len( trait_list ) * len( trait_list ) ) - len( trait_list ) )
    print( "Done. {} values added from {}".format( count, traits_loc ) )

def main( args ):
    # Get inputs if any.
    values = parse_properties( args.D )
    values.update( parse_nexus( args.nexus ) )

    if args.traits is not None:
        for i in args.traits:
            add_traits( values, i, shuffle=args.shuffle )

    env = Environment(
        loader=FileSystemLoader( searchpath="." ),
        autoescape=select_autoescape( ['xml'] )
    )
    template = env.get_template( args.template )
    with open( args.output, "w" ) as output:
        output.write( template.render( values ) )

if __name__ == "__main__" :
    parser = argparse.ArgumentParser( description="Basically a rewrite of BEASTGen..." )

    # Initialize  arguments
    parser.add_argument( "-t", "--template", required=True, help="template file containing kajiki formatting language" )
    parser.add_argument( "-n", "--nexus", required=True, help="nexus file containing alignment and tree" )
    parser.add_argument( "--traits", nargs="+", help="traits file for taxa in nexus file" )
    parser.add_argument( "--shuffle", action="store_true", help="Indicates whether to shuffle the discrete states or not" )
    parser.add_argument( "-D", help="properties for exchange in templates" )
    parser.add_argument( "-o", "--output", required=True, help="location to save formated xml" )
    parser.add_argument( "--sep", required=False, default="|", help="deliminator in sequence names" )
    parser.add_argument( "--order", required=False, default=3, help="index of date in taxa name after splitting by sep" )

    arguments = parser.parse_args()

    main( arguments )