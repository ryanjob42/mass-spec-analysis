# Read the configuration file and create some shorthand variables
# to make it easier to use.
configfile: 'config/config.yaml'

raw_data_dir = config['raw-data-directory']
mzml_dir = config['msconvert-directory']
mzmine_dir = config['mzmine-directory']

# This rule must come first, as it ensures the metadata file exists
# and makes sure MSConvert and mzMine are run as needed.
rule all:
    input:
        # Require that the metadata file exists.
        metadata = config['metadata-file'],

        # Require all the ".mzML" files produced by MSConvert.
        # The file names are the same as one of the raw data ".d" directories,
        # but end in ".mzML" and are found in a different parent directory.
        mzml_files = expand(
            os.path.join(mzml_dir, '{name}.mzML'),
            name=glob_wildcards(os.path.join(raw_data_dir, '{name}.d')).name
        ),

        # Require the results directory from running mzMine.
        mzmine_output = config['mzmine-directory']

# This rule runs MSConvert on a single ".d" directory to produce
# a single ".mzML" file.
rule run_msconvert:
    input: f'{raw_data_dir}/{{experiment}}.d'
    output: f'{mzml_dir}/{{experiment}}.mzML'
    script: 'scripts/msconvert.py'

# This rule runs mzMine on all of the ".mzML" files.
# The list of all ".mzML" files and the directory to put the mzMine results into
# come from the "all" rule. This makes it easier for Snakemake to plan the workflow.
rule run_mzmine:
    input:
        batch_file = config['mzmine-batch-file'],
        mzml_files = rules.all.input.mzml_files
    output: rules.all.input.mzmine_output
    script: 'scripts/mzmine.py'
