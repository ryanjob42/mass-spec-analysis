import logging
import os
import subprocess

from argparse import ArgumentParser
from tempfile import TemporaryDirectory
from typing import List, Optional

# This script runs mzMine.
# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the very end detects which case we're in
# and either runs the "direct_main" or the "snakemake_main" function.
# These functions grab the input and output paths, then use the "run_mzmine"
# function to actually run mzMine.

# Two temporary files are created by this script.
# One contains a list of all the mzML files input.
# The other contains a list of all spectral libraries found in the provided directory.
# These files are given to mzMine via the "--inputs" and "--libraries" arguments.
# This makes it easier (and more consistent) to pass this information to mzMine.
# Note: the spectral libraries are optional. If none are found in the provided directory,
# or if said directory doesn't exist, the list file won't be created, and mzMine won't
# be called with the "--libraries" option. This way, if no libraries were prevented but
# the batch file tries to import them, mzMine will fail. This gives the user a chance to
# either add the library files or to update the batch file to skip the import.

# If you're experiencing issues with this script, uncomment the line below
# to turn on debug logging.
logging.basicConfig(level=logging.DEBUG)

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
    parser.add_argument('spectral_libraries_dir', help='The directory containing spectral libraries to import. If you do not want to import any, provide an empty folder or one that does not exist.')
    parser.add_argument('output_dir', help='The path to the output directory.')
    parser.add_argument('mzml_files', nargs='+', help='A list of mzML files to process, separated by spaces.')
    args = parser.parse_args()
    run_mzmine(args.batch_file, args.mzml_files, args.spectral_libraries_dir, args.output_dir)

def snakemake_main() -> None:
    '''Runs MZmine from data provided by Snakemake.
    Assumes there are two named inputs: "batch_file" and "mzml_files".
    The former is the path to the batch file, and the latter is a list of paths
    to the mzML files to process.
    Also assumes there is a single output which is marked as being a directory.
    '''
    # Grab the input and output paths from the "snakemake" object.
    # While this isn't declared in the current file, Snakemake will create it for us.
    run_mzmine(snakemake.input.batch_file, snakemake.input.mzml_files, snakemake.input.spectral_libs_dir, snakemake.output[0])

def run_mzmine(batch_file: str, mzml_files: List[str], spectral_libs_dir: str, output_dir: str) -> None:
    '''Runs MZmine to process the mzML files.
    @param batch_file The batch file that specifies what MZmine needs to do.
    @param mzml_files The list of mzML files to use as inputs to MZmine.
    @param spectral_libs_dir The directory containing spectra libraries to import. If you do not want to import any, provide an empty folder or one that does not exist.
    @param output_dir The path to the directory where outputs need to be stored.
    '''
    logging.debug('mzMine batch file: %s', batch_file)
    logging.debug('Input mzML files: %s', ', '.join(mzml_files))
    logging.debug('Spectral libraries folder: %s', ', '.join(spectral_libs_dir))
    logging.debug('mzMine output directory: %s', output_dir)

    # Make sure the output directory exists first.
    os.makedirs(output_dir, exist_ok=True)

    with TemporaryDirectory(dir='.') as temp_dir:
        input_file = create_input_file_list(temp_dir, mzml_files)
        spectral_libs_file = create_spectral_libraries_file_list(temp_dir, spectral_libs_dir)
        command = build_mzml_command(batch_file, input_file, spectral_libs_file, output_dir, temp_dir)
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

def create_input_file_list(temp_dir: str, mzml_files: List[str]) -> str:
    '''Creates a file that lists all the mzML files to input into mzMine.
    @param temp_dir The temporary directory for mzMine to use.
    @param mzml_files The list of mzML files to use as inputs to mzMine.
    @return Returns the path to the input file that lists all the mzML files.
    '''
    input_file = os.path.join(temp_dir, 'inputs.txt')
    with open(input_file, 'w') as file:
        # Note: writelines doesn't add newlines to each line.
        file.writelines(f + '\n' for f in mzml_files)
    return input_file

def create_spectral_libraries_file_list(temp_dir: str, spectral_libs_dir: str) -> str:
    '''Creates a file that lists all the spectral library files to import into mzMine.
    @param temp_dir The temporary directory for mzMine to use.
    @param spectral_libs_dir The directory containing the spectral library files to import into mzMine. To skip importing spectral libraries, provide 'None', a folder that does not exist, or a folder with no spectral libraries in it.
    @return If the spectral libraries directory is None, does not exist, or doesn't contain any libraries, returns None. Otherwise, returns the path to a file that lists all of the spectral libraries.
    '''
    if (spectral_libs_dir is None) or (not os.path.isdir(spectral_libs_dir)):
        return None

    # Assume all files in the directory, aside from the README, are libraries.
    library_files = [lib.path for lib in os.scandir(spectral_libs_dir) if lib.name != 'README.md']
    if not library_files:
        return None
    logging.debug('Discovered spectral libraries: %s', ', '.join(library_files))

    libraries_list_file = os.path.join(temp_dir, 'spectral-libraries.txt')
    with open(libraries_list_file, 'w') as file:
        # Note: writelines doesn't add newlines to each line.
        file.writelines(f + '\n' for f in library_files)
    return libraries_list_file

def build_mzml_command(batch_file: str, input_file: str, spectral_libs_file: Optional[str], output_dir: str, temp_dir: str) -> List[str]:
    '''Builds the command needed to run MZmine.
    @param batch_file The batch file that specifies what MZmine needs to do.
    @param input_file The path to the file that lists all the input mzML files.
    @param spectral_libs_file The path to the file that lists all the spectral library files, or None if the import is to be skipped.
    @param output_dir The path to the directory where outputs need to be stored.
    @param temp_dir The temporary directory for MZmine to use.
    @return Returns the list of strings that make up the command to run.
    '''
    # For some reason, MZmine requires the output directory to be formatted as
    # "/output/path/{}", with extra "{}" characters added in.
    output_arg = os.path.join(output_dir, '{}')

    # Determine how many threads to use. If defined, use the number of cpus per task
    # as indicated in the resources. Otherwise, use the number of threads indicated
    # by Snakemake.
    if hasattr(snakemake.resources, 'cpus_per_task'):
        thread_count = snakemake.resources.cpus_per_task
    else:
        thread_count = snakemake.threads

    # Run the given batch file in MZmine.
    # The "--input" argument takes in a file which lists all the mzML files to read,
    # and the files listed in the batch file will be ignored.
    # Likewise, the "--output" argument overwrites where the outputs are stored.
    # The "--temp" argument gives MZmine a temporary directory to work with,
    # as the actual temporary directory might fill up on our cluster.
    # Finally, the "--memory all" argument tells it to use as much memory as it needs
    # to speed up the process, as we have plenty.
    command = [
        'mzmine',
        '--batch', batch_file,
        '--input', input_file,
        '--output', output_arg,
        '--temp', temp_dir,
        '--memory', 'all',
        '--threads', str(thread_count)
    ]

    if spectral_libs_file:
        logging.debug('A spectral library import step is being added to the mzMine command.')
        command.extend(['--libraries', spectral_libs_file])
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
