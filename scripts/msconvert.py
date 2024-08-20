import os.path
import subprocess

from argparse import ArgumentParser
from tempfile import TemporaryDirectory
from typing import List

# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the very end detects which case we're in
# and either runs the "direct_main" or the "snakemake_main" function.
# These functions grab the input and output paths, then use the "run_msconvert"
# function to actually run MSConvert.

# MSConvert's native Linux program is not able to convert the proprietary
# formats (as of Aug-2024). Only the Windows version is able to do that.
# However, they provide a Docker container which can run the Windows
# version of MSConvert using a compatibility layer called "Wine".
# We convert the Docker container to Singularity so it'll run on Slipstick.
# However, this conversion results in a slight permissions issue.
# To fix this, they provide a "mywine" script which creates some files
# into a "/mywineprefix" folder. This folder doesn't normally exist inside
# the container, so we create a temporary directory on Slipstick and mount
# it into the container at that location.
# Note: normally, MSConvert is able to do many conversions at once.
# Unfortunately, the Docker container may run into issues if you do so.
# Thus, this script should only be called for a single conversion.
# However, you can run this script many times in parallel, which effectively
# works around the entire issue.

def direct_main() -> None:
    '''Runs MSConvert from data provided via command-line arguments.
    The first argument is a single input path (to the raw data).
    The second argument is the output file path (to the mzML file expected).
    '''
    # Set up a parser to read in command line arguments, then run MSConvert.
    parser = ArgumentParser(
        prog="MSConvert",
        description="Runs the official MSConvert Docker container through Singularity.")
    parser.add_argument('input_path', help='The path to the raw data to convert.')
    parser.add_argument('output_path', help='The mzML file to output.')
    args = parser.parse_args()
    run_msconvert(args.input_path, args.output_path)

def snakemake_main() -> None:
    '''Runs MSConvert from data provided by Snakemake.
    Assumes there is only a single input path (to the raw data)
    and a single output file path (to the mzML file expected).
    '''
    # Grab the input and output paths from the "snakemake" object.
    # While this isn't declared in the current file, Snakemake will create it for us.
    run_msconvert(snakemake.input[0], snakemake.output[0])

def run_msconvert(input_path: str, output_path: str) -> None:
    '''Runs MSConvert to convert the raw data into an mzML file.
    @param input_path The path to the raw data to convert.
    @param output_path The path to the mzML file to create.
    '''
    # Create a temporary directory for the "mywine" script to use.
    # Python will automatically delete it once it exits the "with" block.
    # We create it in the current directory instead of the default temp directory
    # so it doesn't run out of room and cause this script to crash.
    with TemporaryDirectory(dir='.') as temp_dir:
        command = build_msconvert_command(temp_dir, input_path, output_path)
        _ = subprocess.run(command)

def build_msconvert_command(temp_dir: str, input_path: str, output_path: str) -> List[str]:
    '''Builds the command needed to run MSConvert.
    @param temp_dir The path to a temporary directory MSConvert can use.
    @param input_path The path to the raw data to convert.
    @param output_path The path to the mzML file to create.
    @return Returns the list of strings that make up the command to run.
    '''

    print('RYANRYANRYAN!!!!!!!!!!!!!!!!!!!!!!')
    print('tempdir: ' + temp_dir)
    print('input_p: ' + input_path)
    print('output:  ' + output_path)
    print(f'-B "{temp_dir}:/mywineprefix"')
    print('RYANRYANRYAN!!!!!!!!!!!!!!!!!!!!!!')

    # MSConvert requires that we split the output path into a folder and a file name.
    # If the output folder name is empty, use the current directory (i.e., ".").
    (output_folder, file_name) = os.path.split(output_path)
    if not output_folder:
        output_folder = '.'

    # For detailed descriptions, see:
    # Singularity: https://docs.sylabs.io/guides/latest/user-guide/cli/singularity_exec.html
    # MSConvert Container: https://hub.docker.com/r/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses
    # MSConvert Program: https://proteowizard.sourceforge.io/tools/msconvert.html
    return [
        # We're using Singularity to run the container.
        'singularity', 'exec',

        # Mount the temporary directory into the container as "/mywineprefix".
        '-B', f'{temp_dir}:/mywineprefix',

        # The Docker container for MSConvert.
        'docker://proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses',

        # Use their "mywine" script to run MSConvert.
        'mywine', 'msconvert',

        # Point MSConvert to the input and output.
        input_path,
        '--outdir', output_folder,
        '--outfile', file_name,

        # The output should be an mzML file.
        '--mzML',

        # Use 64-bit precision.
        '--64',

        # Binary data is compressed using "zlib".
        '--zlib',

        # We want to filter by applying peak picking, then use the vendor-supplied
        # centroid algorithms, and only select spectra with an "msLevel" of at least 1.
        # Note: for MSMS data, we could set the "msLevel" to "1-2". However, this only
        # speeds up the process by a couple of seconds and doesn't change the results,
        # so we will skip it to make the script easier to use.
        '--filter', 'peakPicking vendor msLevel=1-',
    ]

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # Otherwise, assume this is being run directly.
    if 'snakemake' in locals():
        snakemake_main()
    else:
        direct_main()
