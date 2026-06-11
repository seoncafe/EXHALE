#!/bin/bash

# Executable bash script to build (via the Makefile) and run the ATES code

echo '=============================================='
echo '|                                            |'
echo '|                                            |'
echo '|            WELCOME TO ATES-2.0             |'
echo '|                                            |'
echo '|                                            |'
echo '=============================================='
echo ''
echo 'For instruction on how to use, please consult'
echo ' https://github.com/AndreaCaldiroli/ATES-Code'

#------- Directories -------#

# Main program directory
DIR_MAIN="$(pwd)"

# Source / utilities directories
DIR_SRC="$DIR_MAIN/src"
DIR_UTILS="$DIR_SRC/utils"

#------- Select Fortran compiler (passed through to the Makefile) -------#
#   ./run_ATES.sh            gfortran (default)
#   ./run_ATES.sh --ifort    Intel classic
#   ./run_ATES.sh --ifx      Intel LLVM
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

python3 -W ignore "$DIR_UTILS/ATES_interface_main.py"

#------- Check input parameters -------#

INPUT_FILE="$DIR_MAIN/input.inp"

# Check if input parameters file exists
echo "Searching for the input parameters file..."
if [ -f "$INPUT_FILE" ]; then
   echo "Input file found. Proceeding..."
else  # Abort if no input.inp exists
   echo "Unable to find a valid input file."
   echo "Please create one through ATES interface."
   exit 1
fi

# Create output directory if it doesn't exist
if [ ! -d "$DIR_MAIN/output" ]; then
   mkdir "$DIR_MAIN/output"
fi

#------- Build ATES.x via the Makefile (incremental) -------#
# Only the sources you edited (and their dependents) are recompiled.
# For a clean rebuild from scratch:  make clean
echo "Building ATES.x ..."
if ! make -C "$DIR_MAIN" $MAKE_FC; then
   echo "Build failed. Aborting."
   exit 1
fi

#------- ... and execute -------#
"$DIR_MAIN/ATES.x"

# Print when execution is over
echo "
----- ATES shutdown -----"
echo '=============================================='
