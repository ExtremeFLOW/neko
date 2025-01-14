#!/usr/bin/bash
export MPICH_GPU_SUPPORT_ENABLED=1
export CUDA_VISIBLE_DEVICES=1
./neko $1
