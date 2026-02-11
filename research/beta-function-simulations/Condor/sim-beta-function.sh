#!/bin/bash

echo $1
Rscript research/beta-function-simulations/5-run-sim-study.R -batch $1
