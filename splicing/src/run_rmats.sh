#!/bin/bash

# Author: Daniel F. Guevara-Diaz

# Usage: ./run_rmats.sh --b1 path/to/b1.txt [--b2 path/to/b2.txt] --nT 15 --oD path/to/output [--paired true|false] [--variable true|false] [--onesample true|false]

# Default values
paired=false
variable=false
onesample=false

# Parse arguments
while [[ "$#" -gt 0 ]]; do
	case $1 in
		--b1-file|--b1) b1_file="$2"; shift ;;
		--b2-file|--b2) b2_file="$2"; shift ;;
		--nThreads|--nT) n_threads="$2"; shift ;;
		--output-dir|--oD) output_dir="$2"; shift ;;
		--paired|-p) paired="$2"; shift ;;
		--variable|-v) variable="$2"; shift ;;
		--onesample|--oS) onesample="$2"; shift;;
		*) echo "Unknown parameters: $1"; exit 1 ;;
	esac
	shift
done

# Print usage if no arguments given
if [[ -z "$b1_file" || -z "$n_threads" || -z "$output_dir" ]]; then
	echo "Usage: $0 --b1 path/to/b1.txt [--b2 path/to/b2.txt] --nT <threads> --oD <path> [--paired true|false] [--variable true|false] [--onesample true|false]"
	exit 1
fi

if [[ "$onesample" != true && -z "$b2_file" ]]; then
	echo "Error: --b2 is required unless --onesample true is specified"
	exit 1
fi

# Fixed paths and variables
gtf_file="/nfs1/Reference/Human_Genome/GRCh38.103/GRCh38.103_150/GRCh38.103.gtf"
read_length=150

# Construct optional flags
paired_flag=""
variable_flag=""
onesample_flag=""
b2_flag=""

if [[ "$paired"  == true ]]; then
	paired_flag="--paired-stats"
fi

if [[ "$variable" == true ]]; then
	variable_flag="--variable-read-length"
fi

if [[ "$onesample" == true ]]; then
	onesample_flag="--statoff"
else
	b2_flag="--b2 $b2_file"
fi

# Run rmats.py
rmats.py \
	--b1 "$b1_file" \
	$b2_flag \
	--gtf "$gtf_file" \
	-t paired \
	--readLength "$read_length" \
	--nthread "$n_threads" \
	--od "$output_dir/output" \
	--tmp "$output_dir/tmp" \
	$paired_flag \
	$variable_flag \
	$onesample_flag
