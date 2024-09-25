import logging
import subprocess

from argparse import ArgumentParser
from tempfile import TemporaryDirectory
from typing import List, Union

# This script can either be run directly, or be run by Snakemake.
# The "if" statement at the very end detects which case we're in
# and either runs the "direct_main" or the "snakemake_main" function.
# These functions grab the input and output paths, then use the "run_skyline"
# function to actually run Skyline.

# If you're experiencing issues with this script, uncomment the line below
# to turn on debug logging.
logging.basicConfig(level=logging.DEBUG)

def direct_main() -> None:
    '''Runs Skyline from data provided via command-line arguments.
    The first argument is the path to the transition list CSV file.
    The second argument is the path to the data to analyze. This should contain the ".d" folders or ".mzML" files.
    The third argument is the type of data to analyze. It should either indicate "raw" for the ".d" folders or "mzML" for ".mzML" files.
    The fourth argument is the type of samples to analyze. It should either indicate "all" to process all files, or "QC" to only process QC and blanks.
    The fifth argument is the number of threads that Skyline should use.
    The sixth argument is the path to output the Skyline document to.
    The seventh argument is the path to output the CSV report to.
    '''
    # Set up a parser to read in command line arguments, then run Skyline.
    parser = ArgumentParser(
        prog='Skyline',
        description='Runs the official Skyline Docker container through Singularity.')
    parser.add_argument('transition_list', help='The path to the transition list CSV file.')
    parser.add_argument('input_dir', help='The path to the data to analyze. This should contain the ".d" folders or ".mzML" files.')
    parser.add_argument('is_centroided', help='If True, centroided data will be assumed (mzML files). If False, TOF profile data will be assumed (raw ".d" folders).')
    parser.add_argument('sample_type', help='The type of samples to analyze. It should either indicate "all" to process all files, or "QC" to only process QC and blanks.')
    parser.add_argument('thread_count', help='The number of threads that Skyline should use.')
    parser.add_argument('doc_path', help='The path to output the Skyline document to.')
    parser.add_argument('report_path', help='The path to output the CSV report to.')

    args = parser.parse_args()
    qc_only = (args.sample_type.lower() == 'qc')
    run_skyline(
        args.transition_list,
        args.input_dir,
        args.is_centroided,
        qc_only,
        args.thread_count,
        args.doc_path,
        args.report_path)

def snakemake_main() -> None:
    '''Runs Skyline from data provided by Snakemake.
    Assumes there is only a single output path (to the CSV report).
    '''
    run_skyline(
        snakemake.input.transition_list,
        snakemake.params.input_dir,
        snakemake.params.is_centroided,
        snakemake.params.qc_only,
        snakemake.params.skyline_threads,
        snakemake.output.skyline_doc,
        snakemake.output.report_file)

def run_skyline(transition_list: str, input_dir: str, is_centroided: bool, qc_only: bool, thread_count: Union[int,str], doc_path: str, report_path: str) -> None:
    '''Runs MSConvert to convert the raw data into an mzML file.
    @param transition_list The path to the transition list CSV file.
    @param input_dir The path to the data to analyze. This should contain the ".d" folders or ".mzML" files.
    @param is_centroided If True, centroided data will be assumed (mzML files). If False, TOF profile data will be assumed (raw ".d" folders).
    @param qc_only If True, only QC and blank samples will be processed. Otherwise, all samples will be processed.
    @param thread_count The number of threads that Skyline should use.
    @param doc_path The path to output the Skyline document to.
    @param report_path The path to output the CSV report to.
    '''
    # Create a temporary directory for the "mywine" script to use.
    # Python will automatically delete it once it exits the "with" block.
    # We create it in the current directory instead of the default temp directory
    # so it doesn't run out of room and cause this script to crash.
    logging.debug('Skyline transition list: %s', transition_list)
    logging.debug('Directory for Skyline inputs: %s', input_dir)
    logging.debug('Analyzing centroided data: %s', str(is_centroided))
    logging.debug('Processing QC/blank samples only: %s', str(qc_only))
    logging.debug('Skyline thread count: %s', str(thread_count))
    logging.debug('Skyline document output path: %s', doc_path)
    logging.debug('Skyline report output path: %s', report_path)

    with TemporaryDirectory(dir='.') as temp_dir:
        command = build_skyline_command(temp_dir, transition_list, input_dir, is_centroided, qc_only, thread_count, doc_path, report_path)
        logging.debug('Skyline command: %s', ' '.join(command))

        # If Skyline returns an error code, the "check=True" argument here
        # will cause this script to throw an exception, letting Snakemake
        # know that it failed.
        _ = subprocess.run(command, check=True)

def build_skyline_command(temp_dir: str, transition_list: str, input_dir: str, is_centroided: bool, qc_only: bool, thread_count: Union[int,str], doc_path: str, report_path: str) -> List[str]:
    '''Builds the command needed to run Skyline.
    @param temp_dir A temporary directory to use for our Singularity container.
    @param transition_list The path to the transition list CSV file.
    @param input_dir The path to the data to analyze. This should contain the ".d" folders or ".mzML" files.
    @param is_centroided If True, centroided data will be assumed (mzML files). If False, TOF profile data will be assumed (raw ".d" folders).
    @param qc_only If True, only QC and blank samples will be processed. Otherwise, all samples will be processed.
    @param thread_count The number of threads that Skyline should use.
    @param doc_path The path to output the Skyline document to.
    @param report_path The path to output the CSV report to.
    @return Returns the list of strings that make up the command to run.
    '''
    # For detailed descriptions, see:
    # Singularity: https://docs.sylabs.io/guides/latest/user-guide/cli/singularity_exec.html
    # Skyline Container: https://hub.docker.com/r/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses

    # Start building the command and adding the options we need.
    command = [
        # We're using Singularity to run the container.
        'singularity', 'exec',

        # Mount the temporary directory into the container as "/mywineprefix".
        '-B', f'{temp_dir}:/mywineprefix',

        # The Docker container for MSConvert.
        'docker://proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses',

        # Use their "mywine" script to run MSConvert.
        'mywine', 'SkylineCMD',

        # Tell Skyline to begin a new document, overwriting a previous one if needed.
        f'--new={doc_path}', '--overwrite',

        # Set the transition settings.
        '--full-scan-precursor-isotopes=Count',
        '--full-scan-precursor-threshold=1',
        '--full-scan-rt-filter=none',
        '--integrate-all=True',
        '--tran-product-ion-types=y,p',
        '--instrument-max-mz=1700',

        # Import the transition list and desired input files.
        f'--import-transition-list={transition_list}',
        f'--import-all={input_dir}',
        
        # Use the desired number of threads for importing data.
        # Note: as of September 2024, using multiple processes from the Singularity
        # container causes the program to stall indefinitely, so we won't use them.
        f'--import-threads={str(thread_count)}',

        # Specify the output report.
        '--report-name=Molecule Transition Results',
        f'--report-file={report_path}',
        '--report-format=csv',
    ]

    # There are some additional transition settings which depend on
    # the type of data we're processing (either mzML or raw data).
    if is_centroided:
        logging.debug('Adding Skyline options for processing centroided data (mzML files).')
        command.extend(centroided_transition_options())
    else:
        logging.debug('Adding Skyline options for processing TOF profile data (raw .d folders).')
        command.extend(tof_profile_transition_options())

    # If we only want the QC and blank files, we need to add more options.
    # Otherwise, if we want all files, we don't need any more options.
    if qc_only:
        logging.debug('Only importing blank and QC files into Skyline.')
        command.extend(qc_only_options())
    else:
        logging.debug('Importing all files into Skyline.')

    return command

def qc_only_options() -> List[str]:
    '''Returns the list of command line options needed if only processing QC samples.'''
    # If we only want the QC and blank samples, use the "--import-filename-pattern" option
    # to only import the files which match the following regular expression.
    # For a visual explanation of the regular expression, paste it into "regexr.com".
    # In short, if we've selected the "qc_only" option,
    # we only want to use files with either "blank" or "QC" in their name.
    # However, we don't want to make assumptions for capitalization, so each letter
    # is written like "[Bb]", which matches either "B" or "b".
    # Each word is in parentheses, with "|" between them.
    # This is what effectively says we need to have either "blank" or "QC".
    # Finally, we start and end with ".*" to indicate we don't care about the rest of the file name.
    import_regex = '.*(([Bb][Ll][Aa][Nn][Kk])|([Qq][Cc])).*'

    # We only need the one option to specify this.
    # However, we'll set it up this way in case that changes at a later date.
    return [f'--import-filename-pattern={import_regex}']

def centroided_transition_options() -> List[str]:
    '''Returns the list of command line options needed for processing centroided data (mzML files).'''
    return [
        '--full-scan-precursor-analyzer=centroided',
        '--full-scan-precursor-res=10',
    ]

def tof_profile_transition_options() -> List[str]:
    '''Returns the list of command line options needed for processing TOF profile data (raw ".d" folders).'''
    return [
        '--full-scan-precursor-analyzer=tof',
        '--full-scan-precursor-res=35000',
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
