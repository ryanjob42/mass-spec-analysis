# Automatic First Time Setup
These are the instructions for setting up the workflow for the first time.
Most of the process has been automated, making it easier to use.
If it does not work for you, please open an Issue on the GitHub page for this workflow, then try the Manual First Time Setup instructions.

This workflow is designed for the Slipstick and Riviera clusters owned by the Data Science Research Institute (DSRI) at Colorado State University (CSU).
For other clusters, some minor changes will likely be needed.

## Run the Automatic Setup Script
These are the instructions for running the automatic setup script.

1. Log in to your HPC cluster (e.g., Slipstick).
2. If you would like the workflow downloaded somewhere specific, use the `cd` command to go there.
3. Run the command below by copying it and pasting it into your terminal and pressing Enter. This will download and run the automatic setup script.
   1. Note: for some terminals, the "Ctrl+V" shortcut will not work. Instead, use the right-click menu (two-finger click on Mac) to paste the command.

```shell
wget -O - https://raw.githubusercontent.com/ryanjob42/mass-spec-analysis/main/scripts/automatic-setup.sh | bash
```

If the script ran successfully, you should see something like this in the output:

```
Automatic installation complete.
Please log out and log back in to finalize the install.
After that, you will still need to manually log in to mzMine.
Please continue following the instructions in the setup document.
```

Once you see that, log off of Slipstick, then continue the instructions below.

## Log In to mzMine
As of August 2024, the command-line way to log in to mzMine does not work.
To see if it has been fixed, you can check the link below.
If it is working, all you should need to do is run `mzmine --login` and follow the instructions.

https://github.com/mzmine/mzmine/issues/2044

If it is not fixed yet, or if you are getting errors, you will need to follow a more complex process.
The link below is for the official instructions.
To make it easier, we have listed them out in more detail here.

https://mzmine.github.io/mzmine_documentation/services/users.html#offline-use

1. Download mzMine to a different device that has an actual screen.
   1. Either Mac or Windows is recommended.
   2. Visit the download page: [https://github.com/mzmine/mzmine/releases/latest](https://github.com/mzmine/mzmine/releases/latest)
   3. Under the "Assets" tab, select the version you want to download.
      1. It appears that pre-M1 Macs are not supported anymore.
      2. There is an "installer" and a "portable" version for each operating system. Use "installer" to install as normal, or "portable" to get a zip file that contains the mzMine program.
   4. Install (or unzip) mzMine like any other program.
2. Launch mzMine.
3. If you have not signed in before:
   1. From the top, select "User > Sign in/Sign up".
   2. Create an account, or log in to your existing account.
4. Run the command below to copy the files to Slipstick, replacing `username` with your Slipstick username.
   1. Mac or Linux: `scp -r "$HOME/.mzmine" username@slipstick.colostate.edu:./`
   2. Windows: `scp -r "%HOMEPATH%\.mzmine" username@slipstick.colostate.edu:./`
5. If that command fails, you can use a program like FileZilla to copy the files. These files should come from the locations below (replacing `username` with your username on the device you're using), and should be copied into the `.mzmine` folder (notice the dot at the beginning) on Slipstick. This will include a `users` folder and a `.mzconfig` file.
   1. If it doesn't appear (or if the program crashes), you can find it manually per below, replacing `username` with your own username.
   2. Windows: `C:\Users\username\.mzmine\`
   3. Mac: `/Users/username/.mzmine/`
   4. Linux: `/home/users/username/.mzmine/`
