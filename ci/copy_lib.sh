#!/bin/bash
#
# Copy the DLL to the python package for bundling and testing.
# ===========================================================
cp bazel-bin/lib/butterlib.dll src/pybutter/
if [ ! -f src/pybutter/butterlib.dll ]; then
    echo "DLL is missing.";
    exit 1;
fi