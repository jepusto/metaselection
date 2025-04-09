FROM rocker/r-ver:4.4.1
MAINTAINER "James Pustejovsky" jepusto@gmail.com

RUN install2.r --error \
    tidyr \
    dplyr \
    purrr \
    scales \
    ggplot2 \
    tictoc \
    Formula \
    MASS \
    mvtnorm \
    optimx \
    nleqslv \
    future.apply \
    progressr \
    clubSandwich \ 
    remotes
    
RUN R -q -e 'remotes::install_github("meghapsimatrix/simhelpers")'
RUN R -q -e 'remotes::install_github("wviechtb/metafor")'
RUN R -q -e 'remotes::install_github("jepusto/metaselection")'
