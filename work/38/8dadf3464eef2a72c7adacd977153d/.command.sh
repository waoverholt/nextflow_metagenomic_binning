#!/bin/bash -ue
spades.py --meta -o KK01_B_1a_err_correct -1 KK01_B_1a_err_correct_R1.fastq -2 KK01_B_1a_err_correct_R2.fastq -t 4     -m 220 --tmp-dir /work/overholt/temp --only-assembler
