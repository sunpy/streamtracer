#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

exts = [Extension(name='streamtracer.fortran.streamtracer',
                  sources=['streamtracer/fortran/Streamtracer.f90'],
                  # f2py_options=['-gomp'],
                  extra_f90_compile_args=['-fopenmp'],
                  extra_link_args=['-fopenmp']
                  ),
        ]

if __name__ == "__main__":
    setup(name='streamtracer',
          version='0.1',
          description='Python wrapped fortran to calculate streamlines',
          url='https://github.com/dstansby/streamtracer',
          author='David Stansby, Lars Mejnertsen',
          author_email='dstansby@gmail.com',
          install_requires=['numpy', 'scipy'],
          python_requires='>=3.6',
          packages=['streamtracer'],
          ext_modules=exts,
          include_package_data=True,
          )
