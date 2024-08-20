import logging
import os
import subprocess

from argparse import ArgumentParser
from tempfile import TemporaryDirectory
from typing import List

# This script runs MZMine.
# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the very end detects which case we're in
# and either runs the "direct_main" or the "snakemake_main" function.
# These functions grab the input and output paths, then use the "run_mzmine"
# function to actually run MZmine.

# If you're experiencing issues with this script, uncomment the line below
# to turn on debug logging.
#logging.basicConfig(level=logging.DEBUG)

def direct_main() -> None:
    '''Runs MZmine from data provided via command-line arguments.
    The first argument is the path to the batch file.
    The second argument is the path to the output directory.
    All remaining arguments are paths to the mzML files to process.
    At least one mzML file is required.
    '''
    parser = ArgumentParser(
        prog="MZmine",
        description="Runs MZmine."
    )
    parser.add_argument('batch_file', help='The path to the batch file.')
    parser.add_argument('output_dir', help='The path to the output directory.')
    parser.add_argument('mzml_files', nargs='+', help='A list of mzML files to process, separated by spaces.')
    args = parser.parse_args()
    run_mzmine(args.batch_file, args.mzml_files, args.output_dir)

def snakemake_main() -> None:
    '''Runs MZmine from data provided by Snakemake.
    Assumes there are two named inputs: "batch_file" and "mzml_files".
    The former is the path to the batch file, and the latter is a list of paths
    to the mzML files to process.
    Also assumes there is a single output which is marked as being a directory.
    '''
    # Grab the input and output paths from the "snakemake" object.
    # While this isn't declared in the current file, Snakemake will create it for us.
    run_mzmine(snakemake.input.batch_file, snakemake.input.mzml_files, snakemake.output[0])

def run_mzmine(batch_file: str, mzml_files: List[str], output_dir: str) -> None:
    '''Runs MZmine to process the mzML files.
    @param batch_file The batch file that specifies what MZmine needs to do.
    @param mzml_files The list of mzML files to use as inputs to MZmine.
    @param output_dir The path to the directory where outputs need to be stored.
    '''
    logging.debug('mzMine batch file: %s', batch_file)
    logging.debug('Input mzML files: %s', ', '.join(mzml_files))
    logging.debug('mzMine output directory: %s', output_dir)

    # Make sure the output directory exists first.
    os.makedirs(output_dir, exist_ok=True)

    with TemporaryDirectory(dir='.') as temp_dir:
        input_file = create_input_file_list(temp_dir, mzml_files)
        command = build_mzml_command(batch_file, input_file, output_dir, temp_dir)
        logging.debug('mzMine command: %s', ', '.join(command))

        # If mzMine returns an error code, the "check=True" argument here
        # will cause this script to throw an exception, letting Snakemake
        # know that it failed.
        _ = subprocess.run(command, check=True)

def create_input_file_list(temp_dir: str, mzml_files: List[str]) -> str:
    '''Creates a file that lists all the mzML files to input into MZmine.
    @param temp_dir The temporary directory for MZmine to use.
    @param mzml_files The list of mzML files to use as inputs to MZmine.
    @return Returns the path to the input file that lists all the mzML files.
    '''
    input_file = os.path.join(temp_dir, 'inputs.txt')
    with open(input_file, 'w') as file:
        # Note: writelines doesn't add newlines to each line.
        file.writelines(f + '\n' for f in mzml_files)
    return input_file

def build_mzml_command(batch_file: str, input_file: str, output_dir: str, temp_dir: str) -> List[str]:
    '''Builds the command needed to run MZmine.
    @param batch_file The batch file that specifies what MZmine needs to do.
    @param input_file The path to the file that lists all the input mzML files.
    @param output_dir The path to the directory where outputs need to be stored.
    @param temp_dir The temporary directory for MZmine to use.
    @return Returns the list of strings that make up the command to run.
    '''
    # For some reason, MZmine requires the output directory to be formatted as
    # "/output/path/{}", with extra "{}" characters added in.
    output_arg = os.path.join(output_dir, '{}')

    # Run the given batch file in MZmine.
    # The "--input" argument takes in a file which lists all the mzML files to read,
    # and the files listed in the batch file will be ignored.
    # Likewise, the "--output" argument overwrites where the outputs are stored.
    # The "--temp" argument gives MZmine a temporary directory to work with,
    # as the actual temporary directory might fill up on our cluster.
    # Finally, the "--memory all" argument tells it to use as much memory as it needs
    # to speed up the process, as we have plenty.
    return [
        'mzmine',
        '--batch', batch_file,
        '--input', input_file,
        '--output', output_arg,
        '--temp', temp_dir,
        '--memory', 'all'
    ]

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # Otherwise, assume this is being run directly.
    if 'snakemake' in locals():
        logging.debug('This script was run by Snakemake.')
        snakemake_main()
    else:
        logging.debug('This script was not run by Snakemake.')
        direct_main()
