#!/bin/bash
export CUDA_VISIBLE_DEVICES=$(( SLURM_LOCALID % 4 ))
exec "$@"
