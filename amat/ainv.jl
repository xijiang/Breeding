#!/usr/bin/env julia
using DelimitedFiles, SparseArrays, LinearAlgebra
D   = readdlm("D.vec")[:,1]
dat = readdlm("T.mat")
T   = sparse(Int.(dat[:,1]), Int.(dat[:,2]), dat[:,3])
dat = nothing

Ai  = T'Diagnal(D)T
