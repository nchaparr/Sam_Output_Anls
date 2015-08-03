from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension('fastfit', ['fastfit.pyx'],
                       include_dirs = [numpy.get_include()],
                       extra_compile_args=['-O3', '-fPIC'],
                       library_dirs=['.'],
                       language="c++")]

setup(name        = 'fastfit',
      cmdclass    = {'build_ext': build_ext},
      ext_modules = ext_modules
      )    

