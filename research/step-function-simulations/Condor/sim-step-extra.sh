#!/bin/bash

echo $1
Rscript research/step-function-simulations/7_big_B_simulations.R -batch $1
