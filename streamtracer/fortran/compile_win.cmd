call F2py -c streamtracer.f90 --f90flags='-fopenmp' -lgomp --fcompiler=gnu95 -m streamtracer

::call F2py -c Streamtracer.f90 --f90flags='/Qopenmp' --fcompiler=intelvem -m streamtracer
