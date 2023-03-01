#!/bin/bash
#BSUB -J join_lite_test
#BSUB -n 144
#BSUB -R span[ptile=72]
#BSUB -R rusage[mem=4]
#BSUB -W 04:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr

source ~/.bashrc
conda activate simenv2
python do_something.py