#!/bin/bash

gfortran -O3 -o gridgen.exe gridgen.f

./gridgen.exe '10' '10000' '1.0d0' '30'
