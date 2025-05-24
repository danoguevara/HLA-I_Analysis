# Splicing Analysis - Type I HLA

## Samples
Tumor samples are from PDX (98% similarity with FLC tumors) and Normal samples come from patients.
Testing both rMATS pipelines with .clean.fastq or .bam files as inputs.

- RU51 (PDX)
- RU57 (PDX)
- RU59 (PDX)
- RU62 (PDX) - Only tumor samples
- RU63 (PDX)
- RU88 (PDX) - Only tumor samples
- RU117 (PDX) - Only tumor samples
- RU123 (PDX) - Only tumor samples
- RU148 (T + N)

## Scripts ran for each sample and all samples
./src/run_rmats.sh --b1 inputs/117T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/117_tumor_variable -v true --oS true &
./src/run_rmats.sh --b1 inputs/123T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/123_tumor_variable -v true --oS true &
./src/run_rmats.sh --b1 inputs/148N_bam.txt --b2 inputs/148T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/148 &
./src/run_rmats.sh --b1 inputs/148N_bam.txt --b2 inputs/148T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/148_variable -v true &
./src/run_rmats.sh --b1 inputs/51N_bam.txt --b2 inputs/51T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/51_variable -v true &
./src/run_rmats.sh --b1 inputs/57N_bam.txt --b2 inputs/57T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/57 &
./src/run_rmats.sh --b1 inputs/57N_bam.txt --b2 inputs/57T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/57_paired -p true &
./src/run_rmats.sh --b1 inputs/59N_bam.txt --b2 inputs/59T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/59_variable -v true &
./src/run_rmats.sh --b1 inputs/62T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/62_tumor --oS true &

./src/run_rmats.sh --b1 inputs/63N_bam.txt --b2 inputs/63T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/63_variable -v true &
./src/run_rmats.sh --b1 inputs/63N_bam.txt --b2 inputs/63T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/63_variable_paired -p true -v true &
./src/run_rmats.sh --b1 inputs/88T_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/88_tumor_variable -v true --oS true &
./src/run_rmats.sh --b1 inputs/all_normal_bam.txt --b2 inputs/all_tumor_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/all_samples &
./src/run_rmats.sh --b1 inputs/all_normal_bam.txt --b2 inputs/all_tumor_bam.txt --nT 25 --oD /nfs1/Splicing_HLA/all_samples_variable -v true

## Internal (renamed) Libraries 
- RU-B (RU-B): 51, 63, 88
- RU-C (RU-B): 59
- RU-G (RU-D): 57, 62
- RU-I (RU-E): 123
- RU-L (RU-F): 117
- RU-N (RU-F): 148

## Meeting 1 (1000 commits to github)
- Run splicing events per samples. (Find recurrent events across samples)
- Run splicing events per sequencing batch. (Fins recurrent events across batches).
- Statistical analysis for all samples and paired.

## Plotting
- Top spliced genes with stacked bars by ASE category.
- Percentage per category.
- Bar plots or dispersion plots with regression.

- Penetrance within patients.
- Consistency between patients.
- Controls for normal tissue against genome.
