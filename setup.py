import os
from numpy.distutils.core import setup, Extension

# Don't want to build the fortran on readthedocs
exts = []
if not os.environ.get('READTHEDOCS', None):
    extra_args = [
        '-fopenmp',
        '-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib'
    ]
    exts += [Extension(name='streamtracer.fortran.streamtracer',
                       sources=['streamtracer/fortran/Streamtracer.f90'],
                       extra_f90_compile_args=extra_args,
                       extra_link_args=extra_args
                       ),
             ]

if __name__ == "__main__":
    setup(ext_modules=exts)
