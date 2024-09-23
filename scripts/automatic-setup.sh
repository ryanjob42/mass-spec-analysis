#!/bin/bash

# Load the user's .bashrc file and the Singularity module.
source "$HOME/.bashrc"
module load singularity

# Use Singularity's "exec" command to cache the MSConvert image.
# We are doing it this way (instead of using the "pull" command)
# so it doesn't create a ".sif" file in the current directory,
# as the workflow doesn't strictly need one.
echo "Ensuring the ProteoWizard (msconvert and Skyline) Singularity image is installed."
singularity exec docker://proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses echo SUCCESS

# Install miniconda if it hasn't been yet.
if [ -d "$HOME/miniconda3" ]; then
    echo "miniconda3 is already installed. Skipping."
else
    # Create a directory for miniconda to be put in.
    mkdir -p "$HOME/miniconda3"

    # Download the miniconda installation script, run it,
    # then delete it once we don't need it.
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O "$HOME/miniconda3/miniconda.sh"
    bash "$HOME/miniconda3/miniconda.sh" -b -u -p "$HOME/miniconda3"
    rm "$HOME/miniconda3/miniconda.sh"

    # Tell miniconda to initialize our Bash shell.
    "$HOME/miniconda3/bin/conda" init bash
fi

# Install mzMine if it hasn't been yet.
if [ -d "$HOME/mzmine" ]; then
    echo "mzMine is already installed. Skipping."
else
    # Create maMine's main folder and its configuration data folder.
    mkdir -p "$HOME/mzmine"
    mkdir -p "$HOME/.mzmine"

    # Download the portable version of mzMine.
    # Note: to update this URL, go to "https://github.com/mzmine/mzmine/releases/latest",
    # then copy the link for the Linux Portable version and paste it here.
    wget https://github.com/mzmine/mzmine/releases/download/v4.3.0/mzmine_Linux_portable-4.3.0.zip -O "$HOME/mzmine/mzmine.zip"

    # Unzip the portable version into the main mzMine folder, delete the zip,
    # and make the mzmine program executable.
    unzip "$HOME/mzmine/mzmine.zip" -d "$HOME/mzmine/"
    rm "$HOME/mzmine/mzmine.zip"
    chmod +x "$HOME/mzmine/bin/mzmine"

    # Add mzMine's bin folder to your path.
    echo 'export PATH="$HOME/mzmine/bin/:$PATH"' >> ~/.bashrc
fi

# Download a copy of this repo if none exists yet.
if [ -d "mass-spec-analysis" ]; then
    echo "The mass-spec-analysis repository is already installed. Skipping."
else
    git clone https://github.com/ryanjob42/mass-spec-analysis.git
fi

# Install the conda environment for this repo if it hasn't been.
"$HOME/miniconda3/bin/conda" env list | grep "mass-spec-analysis" &>/dev/null
if [ $? -eq 0 ]; then
    echo "The mass-spec-analysis conda environment is already installed. Skipping."
else
    "$HOME/miniconda3/bin/conda" env create -f "mass-spec-analysis/environment.yml"
fi

echo
echo Automatic installation complete.
echo Please log out and log back in to finalize the install.
echo After that, you will still need to manually log in to mzMine.
echo See the instructions in the README.md document.
