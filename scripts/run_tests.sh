#!/bin/bash
#
# Using this as the entry point for testing since bazel's
# py_rules are a bit convoluted and not package friendly...
#==========================================================
pytest test/
test_outputs=$()
exit_code=$?
if [ "exit_code" == "0" ]; then
    echo "Test complete!";
    exit 0;
else
    echo "Testing failed."
    echo "{$test_outputs}"
    exit 64;
fi
