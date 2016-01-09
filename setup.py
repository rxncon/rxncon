from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules=cythonize("rxncon/test/simulation/rule_based/cython_cpp/hello_world.pyx",
                            sources=["rxncon/test/simulation/rule_based/cython_cpp/HelloWorld.cpp"],
                            language="c++"))
