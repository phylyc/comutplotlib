#!/bin/bash

# This script is to be executed within the folder.
# First: download and install e.g. miniforge: https://github.com/conda-forge/miniforge
# Then: run this script.
#conda create -n comutplotlib -y python=3.12
#conda activate comutplotlib

# If you want to install the package in an existing environment, just run this line:
python -m pip install --editable .
pyensembl install --release 75 --species human
