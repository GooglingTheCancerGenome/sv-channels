#!/usr/bin/env bash

set -xe

MY_ENV=sv-channels
MY_ENV2=utils

eval "$(conda shell.bash hook)"
printenv
conda env create -n $MY_ENV -f environment.yaml
conda activate $MY_ENV
conda list

conda env create -n $MY_ENV2 -f scripts/utils/environment.yaml
conda activate $MY_ENV2
conda list
