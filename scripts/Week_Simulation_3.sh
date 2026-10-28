#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:mem=20GB

module load app/juliaup/julia
cd ~/full_simulation/scripts
julia -t 10 --project=/home/25262718/full_simulation Week_Simulation_3.jl