# Manual First Time Setup

## Requirements
These are the requirements for running this workflow.

### Singularity
This workflow requires Singularity.
This is already installed on our HPC clusters at CSU.

First, you will need to pull the Singularity container for MSConvert and ProteoWizard.
To do this, run the commands below on the cluster.
Note: if you are wondering why we aren't using `singularity build` or `singularity pull`,
it's because we don't actually need the .sif file created, and just want to make sure it's cached.

```shell
module load singularity
singularity exec docker://proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses echo SUCCESS
```

Notes if you're running on a different cluster than Slistick or Riviera:
* Your system administrators may want you to run these commands within a Slurm job, and not directly on the login node.
* If your cluster also uses Singularity, the `module` command may not be needed.
  * If it gives you an error, try running the `singularity exec` command on its own.
  * If that still gives you an error, ask your system administrator how to use Singularity.
* If your cluster usess Apptainer instead of singularity, you will need to make some changes.
  * First, you'll need to change the command above that caches the container.
  * Second, the Python scripts (in the `scripts` directory) will need to be updated to use Apptainer instead of Singularity.
  * As of September 2024, these changes should be as simple as replacing `singularity` with `apptainer` in the commands.

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
Our workflow uses mzMine for processing.

To install it, or update to the latest version:
1. Go to the website: https://github.com/mzmine/mzmine/releases
2. Skip the "Latest Development Build", and get the first one with a version number.
   1. There should be a small bubble that says "Latest" next to it.
   2. If you want a different version than the latest, find that version instead.
3. If it's not already open, open the "Assets" drop-down for the version you want.
4. Find the one that says (or is similar to) `mzmine_Linux_portable-X.X.X.zip`.
   1. Here, `X.X.X` is the version number.
5. Right-click it (two-finger click on Mac) and select "Copy Link" (or the most similar option depending on your browser).
6. Run either the "First-time install" or "Update an existing install" commands, replacing the URL of the `wget` command with what you've copied.

#### First-time install:
First, run these commands if you're installing mzMine for the first time.
Note: replace the URL for the `wget` command with the URL you copied above.

```shell
mkdir -p ~/mzmine
mkdir -p ~/.mzmine
pushd ~/mzmine
wget https://github.com/mzmine/mzmine/releases/download/v4.3.0/mzmine_Linux_portable-4.3.0.zip -O "$HOME/mzmine/mzmine.zip"
unzip *.zip
chmod +x bin/mzmine
echo export PATH="$HOME/bin/mzmine:$PATH" >> $HOME/.bashrc
source ~/.bashrc
popd
```

Next, you'll need to log in to mzMine.
To avoid duplicating the instructions, see the instructions in the [Automatic First-Time Setup](./AutomaticFirstTimeSetup.md#log-in-to-mzmine) document.


#### Update an existing install:
Run these commands if you're updating an existing installation of mzMine.
Note: replace the URL for the `wget` command with the URL you copied above.

```shell
pushd ~/mzmine
rm -rf ./bin ./lib ./*.zip
wget https://github.com/mzmine/mzmine/releases/download/v4.3.0/mzmine_Linux_portable-4.3.0.zip
unzip ./*.zip
chmod +x bin/mzmine
popd
```
