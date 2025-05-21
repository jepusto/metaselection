#!/bin/bash


k=$(($1 + $2))

echo $k
Rscript simulation-step-function-revised/6_run_sim_study.R -batch $k -iterations 100
