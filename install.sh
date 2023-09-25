#!/bin/bash

# This script is to be executed within the folder.
# First: download and install e.g. miniconda: https://docs.conda.io/projects/miniconda/en/latest/
# Then: run this script.
conda create -n comutplotlib -y python=3.11
conda activate comutplotlib
python -m pip install --editable .
