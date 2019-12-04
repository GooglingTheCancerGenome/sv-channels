#!/usr/bin/env bash

set -xe
printenv
conda env update -n $MY_ENV --file scripts/genome_wide/environment.yaml
