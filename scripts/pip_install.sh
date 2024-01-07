#!/bin/bash
#
# Using an external script to do the pip setup and install the package.
# Once python is setup, we'll copy the DLL into the package data.
#
# Note that this script must be at the root level!
#========================================================================

pip install --upgrade pip setuptools wheel
pip install .
