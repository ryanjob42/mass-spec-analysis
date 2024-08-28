import logging
import os
import subprocess

from argparse import ArgumentParser
from tempfile import TemporaryDirectory
from typing import List, Optional, Union

# This script runs mzMine.
# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the very end detects which case we're in
# and either runs the "direct_main" or the "snakemake_main" function.
# These functions grab the input and output paths, then use the "run_mzmine"
# function to actually run mzMine.

# A valid mzMine batch file must start with a data import step.
# That step would normally define which mzML files to process and the spectral libraries to import.
# However, mzMine will replace these with the contents of the mzML and spectral library file lists provided.
# This is useful, as you probably wrote the batch file on your own computer
# but are trying to run it on an HPC cluster where all the paths are different.

# For the spectral libraries, if you pass in an empty file,
# no libraries will be imported, even if the batch file asks for them.
# On the other hand, if you don't provide a spectral library file but your
# batch file tries to import them, mzMine will still try to import them.
# This will cause it to fail as mzMine won't be able to find the files you wanted.

# If you're experiencing issues with this script, uncomment the line below
# to turn on debug logging.
# logging.basicConfig(level=logging.DEBUG)

def direct_main() -> None:
    '''Runs MZmine from data provided via command-line arguments.
    The first argument is the path to the batch file.
    The second argument is the number of threads to use, or "auto" to let mzMine pick.
    The third argument is the path to the output directory.
    The fourth argument is the path to a ".txt" file listing the mzML files to process (one per line).
    The fifth arugment is an optional path to a ".txt" file listing the spectral libraries to import (one per line).
    '''
    parser = ArgumentParser(
        prog="MZmine",
        description="Runs MZmine."
    )
    parser.add_argument('batch_file', help='The path to the batch file.')
    parser.add_argument('thread_count', help='The number of threads to use (or "auto" to let mzMine pick).')
    parser.add_argument('output_dir', help='The path to the output directory.')
    parser.add_argument('mzml_file_list', help='The path to a ".txt" file listing the mzML files to process (one per line).')
    parser.add_argument('spectral_lib_list', nargs='?', help='Optional. The path to a ".txt" file listing the spectral libraries to import (one per line).')
    args = parser.parse_args()

    run_mzmine(
        args.batch_file,
        args.thread_count,
        args.output_dir,
        args.mzml_file_list,
        args.spectral_lib_list)

def snakemake_main() -> None:
    '''Runs MZmine from data provided by Snakemake.
    Assumes there are two named inputs: "batch_file" and "mzml_files".
    The former is the path to the batch file, and the latter is a list of paths
    to the mzML files to process.
    Also assumes there is a single output which is marked as being a directory.
    '''
    # Determine how many threads to use. If defined, use the number of cpus per task.
    # Otherwise, use the number of threads alloted.
    if hasattr(snakemake.resources, 'cpus_per_task'):
        thread_count = snakemake.resources.cpus_per_task
    else:
        thread_count = snakemake.threads

    # Grab the input and output paths from the "snakemake" object.
    # While this isn't declared in the current file, Snakemake will create it for us.
    run_mzmine(
        snakemake.input.batch_file,
        thread_count,
        snakemake.output[0],
        snakemake.input.mzml_file_list,
        snakemake.input.spectral_lib_list)

def run_mzmine(batch_file: str, thread_count: Union[int, str], output_dir: str, mzml_file_list: str, spectral_lib_list: Optional[str]) -> None:
    '''Runs MZmine to process the mzML files.
    @param batch_file The batch file that specifies what MZmine needs to do.
    @param thread_count The number of threads to use, or "auto" to let mzMine pick.
    @param output_dir The path to the directory where outputs need to be stored.
    @param mzml_file_list The list of mzML files to use as inputs to MZmine.
    @param spectral_lib_list The directory containing spectra libraries to import. If you do not want to import any, provide an empty folder or one that does not exist.
    '''
    logging.debug('mzMine batch file: %s', batch_file)
    logging.debug('mzMine thread count: %s', str(thread_count))
    logging.debug('mzMine output directory: %s', output_dir)
    logging.debug('Input mzML files: %s', ', '.join(mzml_file_list))
    logging.debug('Spectral libraries folder: %s', ', '.join(spectral_lib_list))

    # Make sure the output directory exists first.
    os.makedirs(output_dir, exist_ok=True)

    with TemporaryDirectory(dir='.') as temp_dir:
        command = build_mzml_command(batch_file, thread_count, output_dir, mzml_file_list, spectral_lib_list, temp_dir)
        logging.debug('mzMine command: %s', ' '.join(command))

        # If mzMine returns an error code, the "check=True" argument here
        # will cause this script to throw an exception, letting Snakemake
        # know that it failed.
        _ = subprocess.run(command, check=True)

    # Sometimes, mzMine fails but doesn't return an error code.
    # To work around this, we'll check that the output directory,
    # and if it's empty, let the user know and return an error code.
    if not os.listdir(output_dir):
        logging.error('mzMine did not produce any outputs, which is probably an error.')
        logging.error('This may have occurred for the following reasons:')
        logging.error('The batch file may have a spectral library import step, but no libraries were found by this workflow.')
        logging.error('The batch file might specify blank files incorrectly. If you are using a pattern, make sure you have the correct "*" characters. For example, "*Blank*" finds all files with "Blank" in the name. Alternatively, try specifying the blank files individually. Note: you do not need to worry about the folders being different between your system and this one, mzMine is smart enough to figure that out.')
        logging.error('If you do not expect any output files and want to skip running mzMine, you can comment out or delete the "mzmine_output" line from "rule all" in the "Snakefile". Note: you may also have to delete the comma that comes before it if there is one.')
        exit(1)

def build_mzml_command(batch_file: str, thread_count: Union[int, str], output_dir: str, mzml_file_list: str, spectral_lib_list: Optional[str], temp_dir: str) -> List[str]:
    '''Builds the command needed to run MZmine.
    @param batch_file The batch file that specifies what MZmine needs to do.
    @param thread_count The number of threads to use, or "auto" to let mzMine pick.
    @param output_dir The path to the directory where outputs need to be stored.
    @param mzml_file_list The path to the file that lists all the input mzML files.
    @param spectral_lib_list The path to the file that lists all the spectral library files, or None if the import is to be skipped.
    @param temp_dir The temporary directory for MZmine to use.
    @return Returns the list of strings that make up the command to run.
    '''
    # For some reason, MZmine requires the output directory to be formatted as
    # "/output/path/{}", with extra "{}" characters added in.
    output_arg = os.path.join(output_dir, '{}')

    # Run the given batch file in MZmine.
    # The "--memory all" argument tells it to use as much memory as it needs
    # to speed up the process, as we have plenty.
    # The "--input" argument takes in a file which lists all the mzML files to read,
    # and the files listed in the batch file will be ignored.
    # Likewise, the "--output" argument overwrites where the outputs are stored.
    # The "--temp" argument gives MZmine a temporary directory to work with.
    # We want to specify this, as the actual temporary directory might fill up on our cluster.
    command = [
        'mzmine',
        '--batch', batch_file,
        '--threads', str(thread_count),
        '--memory', 'all',
        '--input', mzml_file_list,
        '--output', output_arg,
        '--temp', temp_dir,
    ]

    # If a file was provided that lists the spectral libraries to import,
    # add that as an argument. Otherwise, leave it off.
    if spectral_lib_list:
        logging.debug('A spectral library import step is being added to the mzMine command.')
        command.extend(['--libraries', spectral_lib_list])
    else:
        logging.debug('Spectral library import is not being added to the mzMine command.')

    return command

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # Otherwise, assume this is being run directly.
    if 'snakemake' in locals():
        logging.debug('This script was run by Snakemake.')
        snakemake_main()
    else:
        logging.debug('This script was not run by Snakemake.')
        direct_main()
