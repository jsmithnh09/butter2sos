#!/bin/bash
#
# Using an external script to do the pip setup and install the package.
# Once python is setup, we'll copy the DLL into the package data.
#
# Note that this script must be at the root level!
#========================================================================
cp bazel-bin/lib/butterlib.dll src/pybutter
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
pip install .
