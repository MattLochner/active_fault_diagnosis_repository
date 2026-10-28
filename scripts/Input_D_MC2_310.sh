#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=20GB

module load app/juliaup/julia
cd ~/full_simulation/scripts
julia -t 1 --project=/home/25262718/full_simulation Input_D_MC_310.jl