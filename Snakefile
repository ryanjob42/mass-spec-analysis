import glob
import os.path

# Read the configuration file and create some shorthand variables
# to make it easier to use, and so there is only one place to make
# changes if the names of the configurations change.
configfile: 'config/config.yaml'

metadata_file = config['metadata-file']
raw_data_dir = config['raw-data-directory']
mzml_out_dir = config['msconvert-output-directory']
batch_file = config['mzmine-batch-file']
spectral_libs_dir = config['mzmine-spectral-libraries-directory']
mzmine_out_dir = config['mzmine-output-directory']

# Make sure the raw data directory exists and has at least one ".d" file.
if not glob.glob(os.path.join(raw_data_dir, '*.d')):
    print(f'No raw data was found in the raw data folder: {raw_data_dir}')
    exit(1)

# This rule must come first, as it ensures the metadata file exists
# and makes sure MSConvert and mzMine are run as needed.
rule all:
    input:
        # Require that the metadata file exists.
        metadata = metadata_file,

        # Require all the ".mzML" files produced by MSConvert.
        # The file names are the same as one of the raw data ".d" directories,
        # but end in ".mzML" and are found in a different parent directory.
        mzml_files = expand(
            os.path.join(mzml_out_dir, '{name}.mzML'),
            name=glob_wildcards(os.path.join(raw_data_dir, '{name}.d')).name
        ),

        # Require the results directory from running mzMine.
        mzmine_output = mzmine_out_dir

# This rule runs MSConvert on a single ".d" directory to produce
# a single ".mzML" file.
rule run_msconvert:
    input: f'{raw_data_dir}/{{experiment}}.d'
    output: f'{mzml_out_dir}/{{experiment}}.mzML'
    script: 'scripts/msconvert.py'

# This rule runs mzMine on all of the ".mzML" files.
# The list of all ".mzML" files and the directory to put the mzMine results into
# come from the "all" rule. This makes it easier for Snakemake to plan the workflow.
rule run_mzmine:
    input:
        batch_file = batch_file,
        mzml_files = rules.all.input.mzml_files,
        spectral_libs_dir = spectral_libs_dir
    output: directory(mzmine_out_dir)
    script: 'scripts/mzmine.py'
