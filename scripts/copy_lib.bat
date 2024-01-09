:: this should be used for copying the DLL into the source directory
:: and for kicking off the tests.
::
:: Similar to the other scripts that are present in this directory,
:: this needs to be executed from the projects root folder.
@echo off
set CI_ROOT_DIR=%cd%
set DLL_SRC_LOC="%CI_ROOT_DIR%\lib\bazel-bin\lib\butterlib.dll"
if not exist %DLL_SRC_LOC% (
    echo "The DLL artifact is missing. Aborting."
    exit 1
)
set DLL_DEST_LOC="%CI_ROOT_DIR%\src\pybutter"
xcopy %DLL_SRC_LOC% %DLL_DEST_LOC% /c /x
