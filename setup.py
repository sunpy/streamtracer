#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

exts = [Extension(name='streamtracer.fortran.streamtracer',
                  sources=['streamtracer/fortran/Streamtracer.f90']
                  ),
        ]

if __name__ == "__main__":
    setup(name='streamtracer',
          version='0.1',
          description='Python wrapped fortran to caclulate streamlines',
          author='Lars Mejnertsen & David Stansby',
          author_email='dstansby@gmail.com',
          install_requires=['numpy', ],
          python_requires='>=3.6',
          packages=['streamtracer', ],
          ext_modules=exts,
          include_package_data=True,
          )
