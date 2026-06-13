#!/bin/bash

# Executable bash script to build (via the Makefile) and run the EXHALE code

echo '=============================================='
echo '|                                            |'
echo '|                                            |'
echo '|             WELCOME TO EXHALE              |'
echo '|                                            |'
echo '|                                            |'
echo '=============================================='
echo ''
echo 'EXHALE is a heavily extended fork of ATES (Caldiroli et al.).'
echo 'For usage, see docs/EXHALE_user_manual.pdf or'
echo ' https://github.com/seoncafe/EXHALE'

#------- Directories -------#

# Main program directory
DIR_MAIN="$(pwd)"

# Source / utilities directories
DIR_SRC="$DIR_MAIN/src"
DIR_UTILS="$DIR_SRC/utils"

#------- Select Fortran compiler (passed through to the Makefile) -------#
#   ./run_EXHALE.sh            gfortran (default)
#   ./run_EXHALE.sh --ifort    Intel classic
#   ./run_EXHALE.sh --ifx      Intel LLVM
case "$1" in
   --ifort) MAKE_FC="FC=ifort" ;;
   --ifx)   MAKE_FC="FC=ifx"   ;;
   *)       MAKE_FC=""         ;;
esac

#------- Call python interface to create input file -------#

TABLE_FILE="$DIR_UTILS/params_table.txt"

if [ ! -f "$TABLE_FILE" ]; then
   python3 -W ignore "$DIR_UTILS/gen_file.py"
   mv "params_table.txt" "$DIR_UTILS/params_table.txt"
fi

python3 -W ignore "$DIR_UTILS/EXHALE_interface_main.py"

#------- Check input parameters -------#

INPUT_FILE="$DIR_MAIN/input.inp"

# Check if input parameters file exists
echo "Searching for the input parameters file..."
if [ -f "$INPUT_FILE" ]; then
   echo "Input file found. Proceeding..."
else  # Abort if no input.inp exists
   echo "Unable to find a valid input file."
   echo "Please create one through the EXHALE interface."
   exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$DIR_MAIN/output" ]; then
   mkdir "$DIR_MAIN/output"
fi

#------- Build EXHALE.x via the Makefile (incremental) -------#
# Only the sources you edited (and their dependents) are recompiled.
# For a clean rebuild from scratch:  make clean
echo "Building EXHALE.x ..."
if ! make -C "$DIR_MAIN" $MAKE_FC; then
   echo "Build failed. Aborting."
   exit 1
fi

#------- ... and execute -------#
"$DIR_MAIN/EXHALE.x"

# Print when execution is over
echo "
----- EXHALE shutdown -----"
echo '=============================================='
