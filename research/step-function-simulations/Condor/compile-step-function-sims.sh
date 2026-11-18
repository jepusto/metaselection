#!/bin/bash

echo $1
Rscript research/step-function-simulations/6-1_compile_bootstrap_sims.R -batch $1
