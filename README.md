# Mass Spectrometry Analysis
A Snakemake workflow for analyzing mass spectrometry data for John Belisle's lab.

## Requirements
These are the requirements for running this workflow.

### Singularity
This workflow requires Singularity.
This is already installed on our HPC clusters at CSU.

TODO: add instructions for pulling the necessary container.
```shell
singularity TODO
```

If you are running on a different cluster that has Apptainer instead,
you can change the commands inside the Python scripts (found in the `scripts` directory)
to use that instead.

### Python and our Conda Environment
This workflow requires you to install Python with Conda.
We recommend using Miniconda, as it is easy to use from our HPC systems.
From the link below, follow the commands under the "Linux" tab to install it.

https://docs.anaconda.com/miniconda/#quick-command-line-install

There is a Conda virtual environment already set up for this workflow.
You can install it using the command below.
If this does not work, or if you need updated software,
see [Updating the Conda Environment](#updating-the-conda-environment) below.

```shell
conda env create -f environment.yml
```

### mzMine
TODO: add instructions for downloading the portable mzMine.

First-time install:
```shell
mkdir -p ~/mzmine
cd ~/mzmine
wget TODO
unzip *.zip
chmod +x bin/mzmine
```

Add the following to the end of your .bashrc file.
TODO: can I add it to the above instructions?
```
export PATH="$HOME/bin/mzmine:$PATH"
```

Updating:
```shell
cd ~/mzmine
rm -rf ./*
wget TODO
unzip *.zip
chmod +x bin/mzmine
```

## Usage
These are the steps to use for running this pipeline.

1. From the command line, download a copy of this workflow using the command below, replacing `my-project` with the name of your project (this will be the name of the folder created).
   1. `git clone https://github.com/ryanjob42/mass-spec-analysis.git my-project`
2. Upload your metadata file and the mzMine batch file to this folder.
   1. You can use FileZilla, or commands like `scp` or `rsync` to do this.
3. Upload your raw data files (the ".d" folders) to the "raw-data" folder here.
   1. If these are on your local computer, you can use FileZilla, or command like `scp` or `rsync`.
   2. If they're on the Berre server, you can use the command below, replacing `source-folder` with the name of the folder that contains the ".d" folders.
      1. `rsync -r berre:/raid/belisle/source-folder raw-data`
4. Update the `config/config.yaml` file as follows:
   1. `metadata-file` should be the name of your metadata file. If it's in a folder, include the folders as well.
   2. `mzmine-batch-file` should be the name of your mzMine batch file. If it's in a folder, include the folders as well.
   3. `raw-data-folder` should be the name of the folder that contains all the raw data (the ".d" folders).
5. Run `start-workflow.sh` to start the workflow.
   1. Note: if you changed the name of the Conda environment, you will need to update it in the script.

## Updating the Conda Environment
The Conda environment file provided here lists the exact version numbers of the Python packages to use.
If there are updates to Snakemake, or any of the other packages, this environment file will not use them.

If you would like to update an existing Conda environment, use the command below.
Note: if you changed the name of the environment, you will need to update it in this command.

```shell
conda update --name mass-spec-analysis --all
```

If you would like to create the environment from scratch, use the commands below.

```shell
conda create -n mass-spec-analysis -c bioconda -c conda-forge -y snakemake
conda activate mass-spec-analysis
pip install snakemake-executor-plugin-slurm
```

After updating (or re-creating) the conda environment, you can update the environment file as follows.
Note: if you changed the name of the environment, you will need to update it in this command.

```shell
conda env export --name mass-spec-analysis --file environment.yml
```
