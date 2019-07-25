#!/usr/bin/env bash

make

echo {1..7} >id.lst

# to calculate inbreeding values for a given list
cat henderson1976.ped |
    ./amat F id.lst >diag.val

# to calculate the inverse of A matrix
cat henderson1976.ped |
    ./amat T id.lst

# the result is in inverse.bin
# one can modify below code for other julia program
./ainv.jl
