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
# logging.basicConfig(level=logging.DEBUG)

##############################################
# Reading Data from Snakemake or Arguments
##############################################

def direct_main() -> None:
    '''Runs Skyline from data provided via command-line arguments.
    The first argument is the path to the transition list CSV file.
    The second argument is the path to the data to analyze. This should contain the ".d" folders or ".mzML" files.
    The third argument is the type of data to analyze. It should either indicate "raw" for the ".d" folders or "mzML" for ".mzML" files.
    The fourth argument is whether mzMine should analyze experimental samples (non-QC and non-blank samples).
    The fifth argument is whether mzMine should analyze QC samples.
    The sixth argument is whether mzMine should analyze blank samples.
    The seventh argument is the number of threads that Skyline should use.
    The eighth argument is the path to output the Skyline document to.
    The ninth argument is the path to output the CSV report to.
    '''
    # Set up a parser to read in command line arguments, then run Skyline.
    parser = ArgumentParser(
        prog='Skyline',
        description='Runs the official Skyline Docker container through Singularity.')
    parser.add_argument('transition_list', help='The path to the transition list CSV file.')
    parser.add_argument('input_dir', help='The path to the data to analyze. This should contain the ".d" folders or ".mzML" files.')
    parser.add_argument('is_centroided', help='If True, centroided data will be assumed (mzML files). If False, TOF profile data will be assumed (raw ".d" folders).')
    parser.add_argument('analyze_experimental', help='"True" if the experimental (non-QC and non-blank) samples should be analyzed, or "False" if they should not be.')
    parser.add_argument('analyze_qc', help='"True" if the QC samples should be analyzed, or "False" if they should not be.')
    parser.add_argument('analyze_blanks', help='"True" if the blank samples should be analyzed, or "False" if they should not be.')
    parser.add_argument('thread_count', help='The number of threads that Skyline should use.')
    parser.add_argument('doc_path', help='The path to output the Skyline document to.')
    parser.add_argument('report_path', help='The path to output the CSV report to.')

    args = parser.parse_args()
    run_skyline(
        args.transition_list,
        args.input_dir,
        as_boolean(args.is_centroided),
        as_boolean(args.analyze_experimental),
        as_boolean(args.analyze_qc),
        as_boolean(args.analyze_blanks),
        args.thread_count,
        args.doc_path,
        args.report_path)

def as_boolean(value: Union[str, bool]) -> bool:
    '''Converts a string containing a boolean to an actual boolean.'''
    return str(value).lower() == 'true'

def snakemake_main() -> None:
    '''Runs Skyline from data provided by Snakemake.
    Assumes there is only a single output path (to the CSV report).
    '''
    run_skyline(
        snakemake.input.transition_list,
        snakemake.params.input_dir,
        snakemake.params.is_centroided,
        snakemake.params.analyze_experimental,
        snakemake.params.analyze_qc,
        snakemake.params.analyze_blanks,
        snakemake.params.skyline_threads,
        snakemake.output.skyline_doc,
        snakemake.output.report_file)


##############################################
# Running Skyline
##############################################

def run_skyline(
        transition_list: str,
        input_dir: str,
        is_centroided: bool,
        analyze_experimental: bool,
        analyze_qc: bool,
        analyze_blanks: bool,
        thread_count: Union[int,str],
        doc_path: str,
        report_path: str) -> None:
    '''Runs MSConvert to convert the raw data into an mzML file.
    @param transition_list The path to the transition list CSV file.
    @param input_dir The path to the data to analyze. This should contain the ".d" folders or ".mzML" files.
    @param is_centroided If True, centroided data will be assumed (mzML files). If False, TOF profile data will be assumed (raw ".d" folders).
    @param analyze_experimental "True" if the experimental (non-QC and non-blank) samples should be analyzed, or "False" if they shouldn't be.
    @param analyze_qc "True" if the QC samples should be analyzed, or "False" if they shouldn't be.
    @param analyze_blanks "True" if the blank samples should be analyzed, or "False" if they shouldn't be.
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
    logging.debug('Analyzing experimental samples: %s', str(analyze_experimental))
    logging.debug('Analyzing QC samples: %s', str(analyze_qc))
    logging.debug('Analyzing blank samples: %s', str(analyze_blanks))
    logging.debug('Skyline thread count: %s', str(thread_count))
    logging.debug('Skyline document output path: %s', doc_path)
    logging.debug('Skyline report output path: %s', report_path)

    with TemporaryDirectory(dir='.') as temp_dir:
        command = build_skyline_command(temp_dir, transition_list, input_dir, is_centroided, analyze_experimental, analyze_qc, analyze_blanks, thread_count, doc_path, report_path)
        logging.debug('Skyline command: %s', ' '.join(command))

        # If Skyline returns an error code, the "check=True" argument here
        # will cause this script to throw an exception, letting Snakemake
        # know that it failed.
        _ = subprocess.run(command, check=True)


##############################################
# Building the Skyline Command
##############################################

def build_skyline_command(
        temp_dir: str,
        transition_list: str,
        input_dir: str,
        is_centroided: bool,
        analyze_experimental: bool,
        analyze_qc: bool,
        analyze_blanks: bool,
        thread_count: Union[int,str],
        doc_path: str,
        report_path: str) -> List[str]:
    '''Builds the command needed to run Skyline.
    @param temp_dir A temporary directory to use for our Singularity container.
    @param transition_list The path to the transition list CSV file.
    @param input_dir The path to the data to analyze. This should contain the ".d" folders or ".mzML" files.
    @param is_centroided If True, centroided data will be assumed (mzML files). If False, TOF profile data will be assumed (raw ".d" folders).
    @param analyze_experimental "True" if the experimental (non-QC and non-blank) samples should be analyzed, or "False" if they shouldn't be.
    @param analyze_qc "True" if the QC samples should be analyzed, or "False" if they shouldn't be.
    @param analyze_blanks "True" if the blank samples should be analyzed, or "False" if they shouldn't be.
    @param thread_count The number of threads that Skyline should use.
    @param doc_path The path to output the Skyline document to.
    @param report_path The path to output the CSV report to.
    @return Returns the list of strings that make up the command to run.
    '''
    # For detailed descriptions, see:
    # Singularity: https://docs.sylabs.io/guides/latest/user-guide/cli/singularity_exec.html
    # Skyline Container: https://hub.docker.com/r/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses

    # To control which types of files to import into Skyline, we need a "regular expression"
    # which will be compared against the input files to determine which to include.
    file_name_pattern = build_file_pattern(analyze_experimental, analyze_qc, analyze_blanks)
    logging.debug('File name pattern: %s', file_name_pattern)

    # Start building the command and adding the options we need.
    # These are all the options which are used no matter what type of analysis is performed.
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
        f'--import-filename-pattern={file_name_pattern}',
        
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

    return command

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


##############################################
# Building the File Name Pattern
##############################################

def build_file_pattern(analyze_experimental: bool, analyze_qc: bool, analyze_blanks: bool) -> str:
    # These are the regular expression pieces for matching either "QC" or "blank"
    # where each letter can be either uppercase or lowercase.
    qc = '[Qq][Cc]'
    blank = '[Bb][Ll][Aa][Nn][Kk]'

    # These are all eight possible combinations of experimental, QC, and blank samples.
    # Each one gets its own regular expression so that Skyline can find appropriate files.
    match analyze_experimental, analyze_qc, analyze_blanks:
        # Can't run Skyline with no inputs.
        case False, False, False:
            print('ERROR: You have selected to run Skyline, but indicated that no files should be processed. This is not allowed.')
            exit(1)

        # Only blank samples.
        case False, False, True:
            return regexp_contains(blank)

        # Only QC samples.
        case False, True, False:
            return regexp_contains(qc)

        # Both blank and QC samples.
        case False, True, True:
            return regexp_contains(regexp_or(blank, qc))

        # Only experimental sampes. That is, samples which are neither QC nor blank samples.
        case True, False, False:
            return regexp_not_contains(regexp_or(blank, qc))

        # Experimental and blank samples. That is, no QC samples.
        case True, False, True:
            return regexp_not_contains(qc)

        # Experimental and QC samples. That is, no blank samples.
        case True, True, False:
            return regexp_not_contains(blank)

        # If we reach this point, that means the user wants all samples.
        case True, True, True:
            return '.*'

def regexp_or(exp1: str, exp2: str) -> str:
    '''Returns a regular expression which will match for either of the given regular expressions.'''
    return f'({exp1})|({exp2})'

def regexp_contains(exp: str) -> str:
    '''Returns a regular expression which matches any text containing the given regular expression.'''
    return f'.*({exp}).*'

def regexp_not_contains(exp: str) -> str:
    '''Returns a regular expression which rejects any text containing the given regular expression.'''
    return f'^(.(?!{exp}))*$'


##############################################
# Miscellaneous
##############################################

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # Otherwise, assume this is being run directly.
    if 'snakemake' in locals():
        logging.debug('This script was run by Snakemake.')
        snakemake_main()
    else:
        logging.debug('This script was not run by Snakemake.')
        direct_main()
