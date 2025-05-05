from setuptools import setup, Extension, find_packages
import os

from Cython.Build import cythonize

PREFIX=os.getenv("CONDA_PREFIX")
include_dirs = [os.path.join(PREFIX, 'include', 'igraph'),
                os.path.join(PREFIX, 'include', 'eigen3'),
                os.path.join(PREFIX, 'include', 'boost'),
                ]
setup(
        name="MCMC",
        ext_modules=cythonize([Extension(name="MCMC._graph_cast",
                                       sources=["_graph_cast.pyx"],
                                       include_dirs=include_dirs,
                                       language="c++",
                                       extra_objects=["-ligraph",],
                                       extra_compile_args=["-Wl,-rpath,/Users/alaink/miniconda3/envs/MCMC/lib "]),]),
        packages = find_packages(),
        install_requires=[
            'numpy',
            'Cython',
            'python-igraph',
            'igraph',
            'eigen',
            'cmake',
            'networkx',
            'scipy',
        ],
        setup_requires = ["cython"],
      )
