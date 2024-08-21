# Spectral Libraries
This is where spectral libraries will go to use within mzMine.
These will be your ".mgf" or ".json" files.
The "config/config.yaml" file must point to where this is for it to work correctly.

All library files here will be imported, even if your batch file doesn't have an import step.
Furthermore, if you have an import step, the files you specified will be ignored, and all the files here will be used.

If you don't put any files here, the workflow will not try to import them.
However, if your batch file has a step to import spectral libraries, since there are none here, mzMine will fail.
