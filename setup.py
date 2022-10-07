from distutils.core import setup, Extension
from Cython.Build import cythonize
from numpy import get_include 

ext = Extension("get_ham", sources=["get_ham.pyx"], include_dirs=['.', get_include()])
setup(name="get_ham", ext_modules=cythonize([ext]))
