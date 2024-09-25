import glob
import os.path

##############################################
# Configuration
##############################################

# Read the configuration file and create some shorthand variables
# to make it easier to use, and so there is only one place to make
# changes if the names of the configurations change.
configfile: 'config.yaml'

METADATA_FILE = config['metadata-file']
RAW_DATA_DIR = config['raw-data-directory']
MZML_OUT_DIR = config['msconvert-output-directory']
BATCH_FILE = config['mzmine-batch-file']
MZMINE_IMPORT_TYPE = config['mzmine-files-to-analyze']
SPECTRAL_LIBS_DIR = config['mzmine-spectral-libraries-directory']
SPECTRAL_LIBS_LIST = config['mzmine-spectral-libraries-import']
MZMINE_OUT_DIR = config['mzmine-output-directory']
SKYLINE_OUT_DIR = config['skyline-output-directory']
SKYLINE_THREAD_COUNT = config['skyline-thread-count']
SKYLINE_PROCESS_RAW = config['skyline-process-raw-data']
SKYLINE_PROCESS_MZML = config['skyline-process-mzml-files']
SKYLINE_IMPORT_TYPE = config['skyline-files-to-analyze']

# We only want the mzMine rules to require a spectral libraries list file
# if the user specified that they want spectral libraries.
USE_SPECTRAL_LIBS = str(SPECTRAL_LIBS_LIST).lower() != 'none'

# Make sure the raw data directory exists and has at least one ".d" file.
# This prevents issues in the event that none are found.
if not glob.glob(os.path.join(RAW_DATA_DIR, '*.d')):
    print(f'No raw data was found in the raw data folder: {RAW_DATA_DIR}')
    print('Please double check you uploaded them and that the configuration file is correct.')
    exit(1)


##############################################
# Top-Level Rule
##############################################

# This rule must come first, as it ensures the metadata file exists
# and makes sure MSConvert and mzMine are run as needed.
rule all:
    input:
        # Require that the metadata file exists.
        metadata = METADATA_FILE,

        # Require all the ".mzML" files produced by MSConvert.
        # The file names are the same as one of the raw data ".d" directories,
        # but end in ".mzML" and are found in a different parent directory.
        mzml_files = expand(
            os.path.join(MZML_OUT_DIR, '{name}.mzML'),
            name=glob_wildcards(os.path.join(RAW_DATA_DIR, '{name}.d')).name
        ),

        # Require the results directory from running mzMine.
        mzmine_output = MZMINE_OUT_DIR

        # If desired, reuqire the Skyline analysis for centroided (mzML)
        # and TOF profile data (raw ".d" folders).
        skyline_centroided_report = branch(SKYLINE_PROCESS_MZML, f'{SKYLINE_OUT_DIR}/centroided_report.csv'),
        skyline_tof_profile_report = branch(SKYLINE_PROCESS_RAW, f'{SKYLINE_OUT_DIR}/tof_profile_report.csv')


##############################################
# MSConvert
##############################################

# This rule runs MSConvert on a single ".d" directory to produce
# a single ".mzML" file.
rule run_msconvert:
    input: f'{RAW_DATA_DIR}/{{experiment}}.d'
    output: f'{MZML_OUT_DIR}/{{experiment}}.mzML'
    script: 'scripts/msconvert.py'


##############################################
# mzMine
##############################################

# This rule prepares a list of mzML files and (optionally) a list of spectral library files
# that mzMine will import. These lists are written to temporary .txt files.
rule prepare_mzmine_inputs:
    input:
        mzml_files = rules.all.input.mzml_files,
        spectral_libs_dir = branch(USE_SPECTRAL_LIBS, SPECTRAL_LIBS_DIR)
    params:
        mzml_import_type = MZMINE_IMPORT_TYPE,
        spectral_libs_list = SPECTRAL_LIBS_LIST
    output:
        mzml_list = temp('mzmine_mzml_list.txt'),
        spectral_libs_list = branch(USE_SPECTRAL_LIBS, temp('mzmine_spectral_libs_list.txt'))
    script: 'scripts/prepare_mzmine_inputs.py'

# This rule runs mzMine on all of the ".mzML" files.
# The list of all ".mzML" files and the directory to put the mzMine results into
# come from the "all" rule. This makes it easier for Snakemake to plan the workflow.
rule run_mzmine:
    input:
        batch_file = BATCH_FILE,
        mzml_file_list = rules.prepare_mzmine_inputs.output.mzml_list,
        spectral_lib_list = rules.prepare_mzmine_inputs.output.spectral_libs_list
    output: directory(MZMINE_OUT_DIR)
    script: 'scripts/mzmine.py'


##############################################
# Skyline
##############################################

rule prepare_transition_list:
    input: MZMINE_OUT_DIR
    output: f'{SKYLINE_OUT_DIR}/transition_list.csv'
    script: 'scripts/prepare_transition_list.py'

rule skyline_centroided_analysis:
    input:
        transition_list = rules.prepare_transition_list.output,
        mzml_files = rules.all.output.mzml_files
    params:
        input_dir = MZML_OUT_DIR,
        is_centroided = True,
        qc_only = (SKYLINE_IMPORT_TYPE.lower() == 'qc'),
        skyline_threads = SKYLINE_THREAD_COUNT
    output:
        skyline_doc = f'{SKYLINE_OUT_DIR}/centroided_analysis.sky',
        report_file = f'{SKYLINE_OUT_DIR}/centroided_report.csv'
    script: 'scripts/skyline.py'

rule skyline_tof_profile_analysis:
    input:
        transition_list = rules.prepare_transition_list.output,
        mzml_files = rules.all.output.mzml_files
    params:
        input_dir = MZML_OUT_DIR,
        is_centroided = False,
        qc_only = (SKYLINE_IMPORT_TYPE.lower() == 'qc'),
        skyline_threads = SKYLINE_THREAD_COUNT
    output:
        skyline_doc = f'{SKYLINE_OUT_DIR}/tof_profile_analysis.sky',
        report_file = f'{SKYLINE_OUT_DIR}/tof_profile_report.csv'
    script: 'scripts/skyline.py
