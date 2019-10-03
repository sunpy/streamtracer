#!/bin/sh

f2py -c Streamtracer.f90 --f90flags='-fopenmp' -lgomp -m streamtracer
#f2py -c Streamtracer.f90 -m streamtracer
