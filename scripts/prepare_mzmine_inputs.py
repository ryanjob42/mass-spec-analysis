import glob
import logging
import os

from argparse import ArgumentParser
from typing import List, Union

# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the very end detects which case we're in
# and either runs the "direct_main" or the "snakemake_main" function.
# These functions grab the input and output paths, then use the "run_skyline"
# function to actually run Skyline.

# If you're experiencing issues with this script, uncomment the line below
# to turn on debug logging.
# logging.basicConfig(level=logging.DEBUG)

##############################################
# Reading Data from Snakemake or Arguments
##############################################

def direct_main():
    '''Creates the mzML and spectral library file lists
    from data provided via command-line arguments.

    The first argument is the path to the directory containing all the mzML files.
    The second argument is the path to the directory containing all the spectral library files.
    The third argument is whether mzMine should analyze experimental samples (non-QC and non-blank samples).
    The fourth argument is whether mzMine should analyze QC samples.
    The fifth argument is whether mzMine should analyze blank samples.
    The sixth argument is the list of spectral libraries to import. Provide "None" to skip importing spectral libraries, "All" to import all discovered libraries, or a comma-separated list of strings for the spectral library files to include.
    The seventh argument is the path to output the mzML file list to.
    The eighth argument is the path to output the spectral library file list to.
    '''
    # Set up a parser to read in command line arguments.
    parser = ArgumentParser(
        prog='Prepare mzMine Inputs',
        description='Creates files listing all the mzML and specctral library files that mzMine should use for analsis.')
    parser.add_argument('mzml_dir', help='The path to the directory containing all the mzML files.')
    parser.add_argument('spectral_libs_dir', help='The path to the directory containing all the spectral library files.')
    parser.add_argument('analyze_experimental', help='"True" if the experimental (non-QC and non-blank) samples should be analyzed, or "False" if they should not be.')
    parser.add_argument('analyze_qc', help='"True" if the QC samples should be analyzed, or "False" if they should not be.')
    parser.add_argument('analyze_blanks', help='"True" if the blank samples should be analyzed, or "False" if they should not be.')
    parser.add_argument('spectral_libs_list', help='The list of spectral libraries to import. Provide "None" to skip importing spectral libraries, "All" to import all discovered libraries, or a comma-separated list of strings for the spectral library files to include.')
    parser.add_argument('mzml_list_path', help='The path to output the mzML file list to.')
    parser.add_argument('spectral_lib_list_path', help='The path to output the spectral library file list to.')
    args = parser.parse_args()

    # For creating the mzML file list, we need to provide a list of mzML files,
    # not just a directory containing the files.
    all_mzml_files = [
        os.path.join(args.mzml_dir, mzml_file)
        for mzml_file in os.listdir(args.mzml_dir)
        if mzml_file.lower().endswith('mzml')
    ]
    create_mzml_list(
        all_mzml_files,
        args.analyze_experimental,
        args.analyze_qc,
        args.analyze_blanks,
        args.mzml_list_path)

    # For creating the spectral library file list, we need to check if the user
    # specifed "None", "All", or a comma-separated list of files.
    # If it's the last, we need to convert that into a list.
    # Otherwise, we keep it as the "None" or "All" string.
    spectral_libs_list = args.spectral_libs_list
    if (spectral_libs_list.lower() != 'none') and (spectral_libs_list.lower() != 'all'):
        spectral_libs_list = spectral_libs_list.split(',')
    create_spectral_libs_list(
        args.spectral_libs_dir,
        spectral_libs_list,
        args.spectral_lib_list_path)

def snakemake_main():
    '''Creates the mzML and spectral library file lists from data provided by Snakemake.'''
    create_mzml_list(
        snakemake.input.mzml_files,
        snakemake.params.analyze_experimental,
        snakemake.params.analyze_qc,
        snakemake.params.analyze_blanks,
        snakemake.output.mzml_list)

    create_spectral_libs_list(
        snakemake.input.spectral_libs_dir,
        snakemake.params.spectral_libs_list,
        snakemake.output.spectral_libs_list)

##############################################
# Creating the mzML file list
##############################################

def create_mzml_list(all_mzml_files: List[str], analyze_experimental: bool, analyze_qc: bool, analyze_blanks: bool, output_path: str) -> None:
    '''Creates a text file which lists all the mzML files that mzMine should analyze (one file per line).
    @param all_mzml_files A list of paths to all the mzML files to consider.
    @param analyze_experimental "True" if the experimental (non-QC and non-blank) samples should be analyzed, or "False" if they shouldn't be.
    @param analyze_qc "True" if the QC samples should be analyzed, or "False" if they shouldn't be.
    @param analyze_blanks "True" if the blank samples should be analyzed, or "False" if they shouldn't be.
    @param output_path The path to output the mzML file list to.
    '''
    logging.debug('All mzML files: %s', ', '.join(all_mzml_files))
    logging.debug('Analyze experimental samples: %s', str(analyze_experimental))
    logging.debug('Analyze QC samples: %s', str(analyze_qc))
    logging.debug('Analyze blank samples: %s', str(analyze_blanks))
    logging.debug('mzML file list output path: %s', output_path)

    # Add the experimental, QC, and blank files to the list of files to analyze.
    to_analyze = []

    if analyze_experimental:
        experimental_files = [f for f in all_mzml_files if is_experimental_sample(f)]
        logging.debug('Experimental files found: %s', ', '.join(experimental_files))
        to_analyze.extend(experimental_files)

    if analyze_qc:
        qc_files = [f for f in all_mzml_files if is_qc_sample(f)]
        logging.debug('QC files found: %s', ', '.join(qc_files))
        to_analyze.extend(qc_files)

    if analyze_blanks:
        blank_files = [f for f in all_mzml_files if is_blank_sample(f)]
        logging.debug('Blank files found: %s', ', '.join(blank_files))
        to_analyze.extend(blank_files)

    # If we didn't find any files to analyze, let the user know
    # and exit with an error code.
    if len(to_analyze) == 0:
        print('ERROR: no mzML files were found to import into mzMine. Your settings were:')
        print(f'Include experimental samples: {analyze_experimental}')
        print(f'Include QC samples: {analyze_qc}')
        print(f'Include blank samples: {analyze_blanks}')
        print(f'The list of files to select from are:')
        for mzml_file in all_mzml_files:
            print(f'  - {mzml_file}')
        exit(1)

    # Finally, write all the discovered mzML files to the desired output.
    write_file(output_path, to_analyze)
    logging.debug('The mzML file list was successfully written.')

def is_qc_sample(mzml_file: str):
    '''Returns True if the mzML file is a QC sample, and False otherwise.'''
    return 'qc' in mzml_file.lower()

def is_blank_sample(mzml_file: str):
    '''Returns True if the mzML file is a blank sample, and False otherwise.'''
    return 'blank' in mzml_file.lower()

def is_experimental_sample(mzml_file: str):
    '''Returns True if the mzML file is an experimental sample (neither QC nor blank), and False otherwise.'''
    return (not is_qc_sample(mzml_file)) and (not is_blank_sample(mzml_file))


##############################################
# Creating the spectral libraries list
##############################################

def create_spectral_libs_list(spectral_libs_dir: str, spectral_libs_list: Union[str, List[str]], output_path: str) -> None:
    '''Creates a text file which lists all the mzML files that mzMine should analyze (one file per line).
    @param spectral_libs_dir The path to the directory containing all the spectral library files.
    @param spectral_libs_list The list of spectral libraries to import. Provide "None" to skip importing spectral libraries, "All" to import all discovered libraries, or a list of strings for the spectral library files to include.
    @param output_path The path to output the spectral library file list to.
    '''
    # If we want to skip the import, don't create a file at all.
    if str(spectral_libs_list).lower() == 'none':
        return

    # If we want to import all the spectral libraries, find all the files
    # in the spectral libs directory which are not Markdown (.md) files.
    elif str(spectral_libs_list).lower() == 'all':
        to_import = [
            f.path
            for f in os.scandir(spectral_libs_dir)
            if f.is_file and not (f.name.lower().endswith('md'))
        ]
        write_file(output_path, to_import)
        return

    # If this is reached, then the user specified a list of specific
    # spectral library files they want to import.
    else:
        to_import = [
            os.path.join(spectral_libs_dir, lib_file)
            for lib_file in spectral_libs_list
        ]
        write_file(output_path, to_import)
        return


##############################################
# Miscellaneous
##############################################

def write_file(path: str, lines: List[str]) -> None:
    '''Writes text to a file, with each piece of text being on a new line.'''
    with open(path, 'w') as out_file:
        # Note: writelines doesn't add newlines.
        out_file.writelines(line + '\n' for line in lines)

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # Otherwise, assume this is being run directly.
    if 'snakemake' in locals():
        logging.debug('This script was run by Snakemake.')
        snakemake_main()
    else:
        logging.debug('This script was not run by Snakemake.')
        direct_main()
