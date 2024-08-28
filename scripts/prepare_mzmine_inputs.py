import logging
import os

from typing import List, Optional

def snakemake_main() -> None:
    # Write a list of all the mzML files to import.
    mzml_list = get_mzml_list()
    write_file(snakemake.output.mzml_list, mzml_list)

    # Get a list of the spectral library files to import.
    # If the list is "None", do nothing.
    # Otherwise, write everything out to a file.
    spec_lib_list = get_spec_lib_list()
    if spec_lib_list is not None:
        write_file(snakemake.output.spectral_libs_list, spec_lib_list)

def get_mzml_list() -> List[str]:
    '''Returns a list of mzML files that mzMine should import.'''
    # "snakemake.inputs.mzml_files" contains a list of all mzML files available.
    # "snakemake.params.mzml_import_type" determines if we're importing all of them
    # or only the pool QC and blanks.
    mzml_list = snakemake.inputs.mzml_files
    import_type = snakemake.params.mzml_import_type

    if import_type.lower() == 'all':
        return mzml_list
    if import_type.lower() == 'qc':
        return [f for f in mzml_list if is_qc_file(f)]

    print(f'Unrecognized mzML import type: {import_type}')
    exit(1)

def is_qc_file(mzml_file: str) -> bool:
    '''Returns True if the given mzML file is a pool QC or blank sample.'''
    if 'qc' in mzml_file.lower():
        return True
    if 'blank' in mzml_file.lower():
        return True
    return False

def write_file(path: str, lines: List[str]) -> None:
    '''Writes text to a file, with each piece of text being on a new line.'''
    with open(path, 'w') as out_file:
        # Note: writelines doesn't add newlines.
        out_file.writelines(line + '\n' for line in lines)

def get_spec_lib_list() -> Optional[List[str]]:
    '''
    Returns a list of spectral libraries that mzMine should import.
    If "None" is returned, don't create a file.
    If an empty list is returned, create an empty file.
    If a non-empty list is returned, create a file with each element on its own line.
    '''
    spec_lib_dir = snakemake.inputs.spectral_libs_dir
    import_type = snakemake.params.spectral_libs_list

    # The import type can either be "None", "All" or a list of strings.
    # To make sure ".lower()" is available, we need to first ensure the
    # import type is a string. Then we can reliably check if it's "none" or "all".
    import_type_as_str = str(import_type).lower()
    if str(import_type).lower() == 'none':
        return None
    elif str(import_type).lower() == 'all':
        return [
            f.path
            for f in os.scandir(spec_lib_dir)
            if f.is_file() and (f.name != 'README.md')
        ]
    else:
        return [
            os.path.join(spec_lib_dir, f)
            for f in import_type
        ]

if __name__ == '__main__':
    # If the "snakemake" variable is defined, assume this is being run by Snakemake.
    # This file is currently not set up to run without Snakemake.
    if 'snakemake' in locals():
        logging.debug('This script was run by Snakemake.')
        snakemake_main()
    else:
        print('This file must be run through Snakemake.')
        exit(1)
