#!/bin/bash
source "/home/dbulgakov/intel/oneapi/vtune/latest/env/vars.sh" && vtune-gui --project-path "/home/dbulgakov/Study/hpc/vtune_no_omp" --app-path "/home/dbulgakov/Study/hpc/src/build/src_no_omp"