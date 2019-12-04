#!/usr/bin/env bash

set -x
conda env update -n $CONDA_DEFAULT_ENV --file scripts/genome_wide/environment.yaml
