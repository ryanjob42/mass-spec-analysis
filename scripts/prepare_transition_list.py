import glob
import logging
import os.path
import pandas as pd

from argparse import ArgumentParser

# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the very end detects which case we're in
# and either runs the "direct_main" or the "snakemake_main" function.
# These functions grab the input and output paths, then use the "run_skyline"
# function to actually run Skyline.

# If you're experiencing issues with this script, uncomment the line below
# to turn on debug logging.
#logging.basicConfig(level=logging.DEBUG)

def direct_main() -> None:
    '''Builds the transition list from data provided via command-line arguments.
    The first argument is the path to the "quant_full.csv" file produced by mzMine.
    The second argument is the path to output the transition list CSV file to.
    '''
    parser = ArgumentParser(
        prog='Transition List Builder',
        description='Builds a transition list CSV file from the "quant_full.csv" output of mzMine.')
    parser.add_argument('input_path', help='The path to the "quant_full.csv" file produced by mzMine.')
    parser.add_argument('output_path', help='The path to output the transition list CSV file to.')

    args = parser.parse_args()
    build_transition_list(args.input_path, args.output_path)

def snakemake_main() -> None:
    '''Builds the transition list from data provided by Snakemake.
    Assumes there is a single input, which is the path to mzMine's output directory.
    Assumes there is a single output, which is the path for the transition list CSV file.'''
    # As input, this Snakemake rule takes the entire folder of data
    # that mzMine produced. However, we only want the CSV file that
    # contains the text "quant_full", as that's the correct one
    # to build the transition list from.
    quant_full_path = find_quant_full_file(snakemake.input[0])
    build_transition_list(quant_full_path, snakemake.output[0])

def find_quant_full_file(mzmine_out_dir: str) -> str:
    '''Tries to find the correct "quant_full.csv" file within mzMine's output directory.
    If no such file is found, or if multiple are found, a message will be printed for the user
    and this script will immediately exit.

    @param mzmine_out_dir The path to mzMine's output directory.
    @returns The path to the "quant_full.csv" file.
    '''

    # Use "glob" to find all CSV files that contain the text "quant_full" in their name.
    # This returns a list of all matching files.
    found_files = glob.glob('*quant_full*.csv', root_dir=mzmine_out_dir)

    # We expect exactly one file to match. If no files match, or if multiple match,
    # let the user know and exit.
    if len(found_files) < 1:
        print('Did not find the "quant_full.csv" file within the mzMine output directory!')
        exit(1)
    if len(found_files) > 1:
        print('Multiple "quant_full.csv" files were found in the mzMine output directory, and it is unclear which to use!')
        exit(1)

    # Return the path to that file by merging the mzMine output directory
    # with the path that "glob" gave us.
    return os.path.join(mzmine_out_dir, found_files[0])

def build_transition_list(input_path: str, output_path: str) -> None:
    '''Builds a transition list CSV file from the output of mzMine,
    specifically the "quant_full.csv" file.

    The transition list will contain four columns as follows:
    1. "Precursor m/z": A copy of the "mz" column from the input file.
    2. "Molecule Name": A unique ID created from the "mz" and "rt" columns.
    3. "Precursor Charge": All values in this column are assumed to be "1".
    4. "Molecule List Name": This column is present, but empty.

    @param input_path The path to the "quant_full.csv" file produced by mzMine.
    @param output_path The path to output the transition list to.
    '''
    logging.debug('mzMine "quant_full.csv" file path: %s', input_path)
    logging.debug('Transition list output path: %s', output_path)

    # Read in the source data as a Pandas DataFrame.
    src_data = pd.read_csv(input_path)

    # Start a new DataFrame from scratch and copy the "mz" column over.
    transition_list = pd.DataFrame()
    transition_list['Precursor m/z'] = src_data['mz']

    # Round the "mz" and "rt" columns, then combine them into the new
    # "Molecule Name" column.
    mz_rounded = src_data['mz'].apply(round_mz_values)
    rt_rounded = src_data['rt'].apply(round_rt_values)
    transition_list['Molecule Name'] = 'MZ' + mz_rounded + '_' + rt_rounded

    # Set all the precursor charges to 1 and create the empty "Molecule List Name" column.
    transition_list['Precursor Charge'] = 1
    transition_list['Molecule List Name'] = ''

    # Output the transition list to the desired CSV file.
    # Do not include the index, as it causes Skyline to fail.
    transition_list.to_csv(output_path, index=False)

def round_mz_values(mz_value: float) -> str:
    '''Converts a precursor m/z value to a string rounded to 4 decimal places.'''
    return f'{mz_value:0.4f}'

def round_rt_values(rt_value: float) -> str:
    '''Converts a retention time to a string rounded to 2 decimal places.'''
    return f'{rt_value:0.2f}'

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # Otherwise, assume this is being run directly.
    if 'snakemake' in locals():
        logging.debug('This script was run by Snakemake.')
        snakemake_main()
    else:
        logging.debug('This script was not run by Snakemake.')
        direct_main()
