#!/usr/bin/env bash

echo {1..7} >id.lst

cat henderson1976.ped |
    ./amat F id.lst >diag.val
