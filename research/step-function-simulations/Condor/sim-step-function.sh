#!/bin/bash

echo $1
Rscript research/step-function-simulations/5_run_sim_study.R -batch $1
