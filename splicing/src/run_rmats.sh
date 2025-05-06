#!/bin/bash

# Author: Daniel F. Guevara-Diaz

# Usage: ./run_rmats.sh --b1 path/to/b1.txt --b2 path/to/b2.txt --rL 150 --nT 8 --oD path/to/output [--paired true|false] [--variable true|false]

# Default values
paired=false
variable=false

# Parse arguments
while [[ "$#" -gt 0 ]]; do
	case $1 in
		--b1-file|--b1) b1_file="$2"; shift ;;
		--b2-file|--b2) b2_file="$2"; shift ;;
		--readLength|--rL) read_length="$2"; shift ;;
		--nThreads|--nT) n_threads="$2"; shift ;;
		--output-dir|--oD) output_dir="$2"; shift ;;
		--paired|-p) paired="$2"; shift ;;
		--variable|-v) variable="$2"; shift ;;
		*) echo "Unknown parameters: $1"; exit 1 ;;
	esac
	shift
done

# Print usage if no arguments given
if [[ -z "$b1_file" || -z "$b2_file" || -z "$read_length" || -z "$n_threads" || -z "$output_dir" ]]; then
	echo "Usage: $0 --b1 path/to/b1.txt --b2 path/to/b2.txt --rL <read_length> --nT <threads> --oD <path> [--paired true|false] [--variable true|false]"
	exit 1
fi

# Fixed paths
gtf_file="/nfs1/Reference/Human_Genome/GRCh38.103/GRCh38.103_150/GRCh38.103.gtf"

# Construct optional flags
paired_flag=""
variable_flag=""

if [[ "$paired" == true ]]; then
	paired_flag="-t paired"
fi

if [[ "$variable" == true ]]; then
	variable_flag="--variable-read-length"
fi

# Run rmats.py
rmats.py \
	--b1 "$b1_file" \
	--b2 "$b2_file" \
	--gtf "$gtf_file" \
	--readLength "$read_length" \
	--nthread "$n_threads" \
	--od "$output_dir/output" \
	--tmp "$output_dir/tmp" \
	$paired_flag \
	$variable_flag
