# Mass Spectrometry Analysis
A Snakemake workflow for analyzing mass spectrometry data for John Belisle's lab.

## Setup
There is a one-time setup process for this workflow.
See the [Automatic First Time Setup](./docs/AutomaticFirstTimeSetup.md) instructions for more details.

## Usage
There are four important pieces to using this workflow.

First, you will want to update the `config.yaml` configuration file.
This is a text file that you can use to tell the workflow what you want to do.
The file is commented with some instructions, including a list of settings that are commonly changed.
These include:

- `metadata-file`: The path to the metadata file for your experiment.
- `mzmine-batch-file`: The path to the batch file you are using for running mzMine.
- `mzmine-spectral-libraries-import`: Indicate which spectral libraries to use when running mzMine.

Second, upload your metadata file.
By default, it is expected to be named `metadata.csv` and put in the main folder for the workflow (the same place where the "Snakemake" folder is).
If you would like to change the name or location, you can change the `metadata-file` setting in the configuration file.

Third, upload your raw mass spectrometry data to the raw data directory.
For Agilent systems, this will be your ".d" folders.
By default, this is the `raw-data` folder here, but you can configure it with the `raw-data-directory` setting in the configuration file.
When the workflow runs, it will automatically convert everything there into ".mzML" files.

Fourth, upload your mzMine batch file.
You can create this on any other system with mzMine installed.
By default, it's expected to be named `batch.xml` and put in the main folder for the workflow (the same place where the "Snakemake" folder is).
If you would like to change the name or location, you can change the `mzmine-batch-file` setting in the configuration file.

Fifth, upload your spectral libary files to the spectral libraries directory.
These will either be ".json" or ".mgf" files.
By default, this is the `spectral-libraries` folder here, but you can configure it with the `mzmine-spectral-libraries-directory` setting in the configuration file.
Unlike the raw data folder, you can manually specify (in the configuraiton file) which libraries, if any, to import.

Finally, once everything is done, you can simply run the `start-workflow.sh` script.
This contains the couple of commands needed to run Snakemake.
While it's running, you must leave the terminal window open.
The workflow typically takes about 10-15 minutes to complete.
