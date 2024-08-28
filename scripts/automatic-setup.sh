#!/bin/bash

# Load the user's .bashrc file and the Singularity module.
source "$HOME/.bashrc"
module load singularity

# Use Singularity's "exec" command to cache the MSConvert image.
# We are doing it this way (instead of using the "pull" command)
# so it doesn't create a ".sif" file in the current directory,
# as the workflow doesn't strictly need one.
singularity exec docker://proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses echo SUCCESS

# Install miniconda if it hasn't been yet.
if [ ! -d "$HOME/miniconda3" ]; then
    mkdir -p "$HOME/miniconda3"
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$HOME/miniconda3/miniconda.sh"
    bash "$HOME/miniconda3/miniconda.sh" -b -u -p "$HOME/miniconda3"
    rm "$HOME/miniconda3/miniconda.sh"
    source "$HOME/.bashrc"
    "$HOME/miniconda3/bin/conda init bash"
fi

# Install the conda environment.
conda env create -f environment.yml

# Install mzMine if it hasn't been yet.
if [ ! -d "$HOME/mzmine" ]; then
    mkdir -p "$HOME/.mzmine"
    mkdir -p "$HOME/mzmine"
    wget https://github.com/mzmine/mzmine/releases/download/v4.2.0/mzmine_Linux_portable-4.2.0.zip -O "$HOME/mzmine/mzmine.zip"
    unzip "$HOME/mzmine/mzmine.zip" -d "$HOME/mzmine/"
    rm "$HOME/mzmine/mzmine.zip"
    chmod +x "$HOME/mzmine/bin/mzmine"
fi

echo
echo Automatic installation complete.
echo Please log out and log back in to finalize the install.
echo After that, you will still need to manually log in to mzMine.
echo See the instructions in the README.md document.
