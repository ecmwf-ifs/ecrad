#!/usr/bin/env python3

import argparse
from datetime import datetime
from pathlib import Path
import shutil
import sys

def pprint_update(ecrad_path, ecrad_files, ifs_path, ifs_files):
    updates = [
        f'{file.relative_to(ifs_path)} -> {ecrad_files[file.name].relative_to(ecrad_path)}'
        for file in ifs_files.values()
    ]
    return '\n'.join(updates)

# CLI arguments
parser = argparse.ArgumentParser(prog=Path(__file__).name)
parser.add_argument('path_to_ifs_source', help='Path to ifs-source worktree')
parser.add_argument('--apply', help='Copy files from ifs-source', action='store_true')
args = parser.parse_args()

# Base paths
ecrad_path = Path(__file__).parent.parent.resolve()
ifs_path = Path(args.path_to_ifs_source)
logpath = Path(__file__).with_suffix('.log')

# Build list of files to update
files_to_update = [ecrad_path/'radiation/radiation_ifs_rrtm.F90']
files_to_update += list(ecrad_path.glob('ifsrrtm/*.F90'))
files_to_update += list(ecrad_path.glob('ifs/*.F90'))
ecrad_files = {file.name: file for file in files_to_update}

# Find files in ifs-source
ifs_files = {
    ifs_file.name: ifs_file
    for ifs_file in ifs_path.glob('**/*.F90')
    if ifs_file.name in ecrad_files
}
assert all(f in ifs_files for f in ecrad_files)

# Print the planned update
print('Files to update:')
print('')
print(pprint_update(ecrad_path, ecrad_files, ifs_path, ifs_files))

# Apply changes
if args.apply:
    print('Copying files...')
    for file in ifs_files.values():
        shutil.copy(file, ecrad_files[file.name])
    print('Done.')

    with logpath.open('w') as logfile:
        logfile.write(f'Updated files on {datetime.now()}:\n')
        logfile.write('\n')
        logfile.write(pprint_update(ecrad_path, ecrad_files, ifs_path, ifs_files))

    print(f'Logfile written to {logpath}.')

sys.exit(0)
