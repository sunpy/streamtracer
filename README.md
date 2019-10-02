### Building

The fortran sources will need to be built first; for this, you will need a fortran compiler.
To see the compilers that are available, run the following command:

```sh
f2py -c --help-fcompiler
```
If the compiler you want is not available, but is installed on the system, you may need to load that particular compiler.

Intel and Windows seems notorious for its hard to load/find compilers, so it may be easier to use gfortran.
However, if Intel is required, find and run `compilervars.bat` in your command line.
This is usually found at `<intel_install_dir>\compilers_and_libraries\windows\bin\compilervars.bat`

Once the compiler required is available, run one of the following commands:

| OS | Compiler | Command |
| :-- | :-------- | :------- |
| Windows | Intel | `python setup.py config_fc --fcompiler intelvem --f90flags "/Qopenmp" build_src build_ext` |
| Windows | gfortran | `python setup.py config_fc --fcompiler gnu95 --f90flags "-fopenmp" build_src build_ext -lgomp` |
| Linux | Intel | `python setup.py config_fc --fcompiler intelem --f90flags "-fopenmp" build_src build_ext` |
| Linux | gfortran | `python setup.py config_fc --fcompiler gnu95 --f90flags "-fopenmp" build_src build_ext -lgomp` |

*Note, for intel, the fcompiler may be a slightly different name. Change to whichever is available.*

### Installation

Once built, it can simply be installed by running:

```sh
python setup.py install
```

Once installed, the package can be used by importing:

```python
import streamtracer
```

See the examples for more use cases
