#!/usr/bin/env python
import os, sys
import numpy
from os.path import join as pjoin
import shutil


try:
    from setuptools import setup, Extension, Command
    from setuptools.command.build_ext import build_ext as _build_ext
    from setuptools.command.build import build
except ImportError:
    from distutils.core import setup, Extension, Command
    from distutils.command.build_ext import build_ext as _build_ext
    from distutils.command.build import build
    from distutils.command.sdist import sdist as _sdist

class CleanCommand(Command):
    """Custom distutils command to clean the .so and .pyc files."""

    user_options = [("all", "a", "")]

    def initialize_options(self):
        self.all = True
        self._clean_me = []
        self._clean_trees = []
        self._clean_exclude = ["mpfit.c"]

        for root, dirs, files in list(os.walk('kmos2d')):
            for f in files:
                if f in self._clean_exclude:
                    continue
                if os.path.splitext(f)[-1] in ('.pyc', '.so', '.o',
                                               '.pyo',
                                               '.pyd', '.c', '.orig'):
                    self._clean_me.append(pjoin(root, f))
            for d in dirs:
                if d == '__pycache__':
                    self._clean_trees.append(pjoin(root, d))

        for d in ('build',):
            if os.path.exists(d):
                self._clean_trees.append(d)

    def finalize_options(self):
        pass

    def run(self):
        for clean_me in self._clean_me:
            try:
                os.unlink(clean_me)
            except Exception:
                pass
        for clean_tree in self._clean_trees:
            try:
                import shutil
                shutil.rmtree(clean_tree)
            except Exception:
                pass




try:
    from Cython.Build import cythonize
except ImportError:
    raise ImportError("pySPARK requires cython to install")

##compile the cython modules as necessary
#class sdist(_sdist):
#    def run(self):
#        try:
#            from Cython.Build import cythonize
#        except ImportError:
#            raise ImportError("kmos2d_drp requires cython to install")
#        else:
#            cythonize(['kmos2d/cbroaden.pyx', 
#                       'kmos2d/fit_sky.pyx',
#                       'kmos2d/cmpfit.pyx',
#                       'kmos2d/sincrebin.pyx'])
#            _sdist.run(self)

class build_ext(_build_ext):
    def build_extension(self, ext):
        _build_ext.build_extension(self, ext)
        

if __name__ == "__main__":

    include_dirs = [
            "include",
            numpy.get_include(),
            ]

    #cb_ext = Extension("kmos2d.cbroaden",   sources=["kmos2d/cbroaden.pyx"], include_dirs=include_dirs)
    #fs_ext = Extension("kmos2d.fit_sky",    sources=["kmos2d/fit_sky.pyx"], include_dirs=include_dirs)
    #mp_ext = Extension("kmos2d.cmpfit",     sources=["kmos2d/cmpfit.pyx","kmos2d/mpfit.c"], include_dirs=include_dirs)
    #sr_ext = Extension("kmos2d.sincrebin",  sources=["kmos2d/sincrebin.pyx"], include_dirs=include_dirs)

    #ext_modules = cythonize([cb_ext, fs_ext, mp_ext, sr_ext])

    #scripts = ['scripts/'+file for file in os.listdir('scripts/')]  
    scripts = []

    cmdclass = {'clean': CleanCommand,
                'build_ext': build_ext}

    if sys.version_info[0] < 3:
        import __builtin__ as builtins
    else:
        import builtins
    builtins.__KMOS2D_SETUP__ = True
    import pyspark

    setup(
        name = "pyspark",
        version=pyspark.__version__,
        author="Trevor Mendel",
        #ext_modules = ext_modules,
        cmdclass = cmdclass,
        scripts = scripts, 
        packages=["pyspark",
                  "pyspark.detector",
                  "pyspark.utils",
                  "pyspark.calibration",
                  "pyspark.cubes",
                  "pyspark.sky"],
        package_data={"pyspark": ["data/superdark_coeffs_mean.fits.fz", "data/superbad_60s.fits.fz"]},
        include_package_data=True,
        zip_safe=False
    )

