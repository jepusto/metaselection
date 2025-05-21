#!/bin/bash

echo $1
Rscript research/beta-function-simulations/6_run_sim_study.R -batch $1
