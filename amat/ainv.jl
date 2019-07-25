#!/usr/bin/env julia
using DelimitedFiles, SparseArrays, LinearAlgebra, Serialization

D   = readdlm("D.vec")[:,1]
dat = readdlm("T.mat")
T   = sparse(Int.(dat[:,1]), Int.(dat[:,2]), dat[:,3])
dat = nothing

Ai  = T'Diagonal(D)T

# An efficient way for data I/O
# write the inverse to a file, say inverse.bin
file = open("inverse.bin", "w")
serialize(file, Ai)
close(file)

# read them back
Ai = deserialize(open("inverse.bin"))
