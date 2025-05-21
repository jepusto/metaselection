FROM rocker/r-ver:4.4.0
MAINTAINER "James Pustejovsky" jepusto@gmail.com

RUN install2.r --error \
    tidyr \
    dplyr \
    purrr \
    tictoc \
    Formula \
    MASS \
    mvtnorm \
    optimx \
    nleqslv \
    rootSolve \
    metafor \
    clubSandwich \ 
    remotes
    
RUN R -q -e 'remotes::install_github("meghapsimatrix/simhelpers")'
