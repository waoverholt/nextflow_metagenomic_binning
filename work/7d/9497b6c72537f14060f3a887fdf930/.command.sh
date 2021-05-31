#!/bin/bash -ue
spades.py --meta -o KK_M_7b_err_correct -1 KK_M_7b_err_correct_R1.fastq -2 KK_M_7b_err_correct_R2.fastq -t 4     -m 220 --tmp-dir /work/overholt/temp --only-assembler
