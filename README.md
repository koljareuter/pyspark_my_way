# pySPARK
Reduction scripts for processing KMOS IFU data.  This will eventually include routines for detector-level processing

Current version: v0.9.0b0

# References
Eventually there will be a paper describing some of this...

# Installation
If you have git installed, then the reduction routines can be obtained using the following command:
```
git clone https:/github.com/jtmendel/pyspark.git pyspark
cd pyspark
python setup.py install
```
This will compile the necessary routines and install them to their appropriate locations.  There are various options controlling *where* those locations might be, with perhaps the most germane being `--prefix=/path/to/installation` and `--install-scripts=/directory/in/your/path`.  Additional install options can be perused using:

```
python setup.py install --help
```
Note that you also need a number of additional bits installed in order to run the reduction procedures:
* [SPARK](https://www.eso.org/sci/software/pipelines/kmos/kmos-pipe-recipes.html) - tested against v1.3.5 of the ESO pipeline.
* [Cython](http://cython.org/) - used for speeding up the interpolation routines

As well as the following Python modules:
* [bottleneck](https://pypi.python.org/pypi/Bottleneck) - fast methods for dealing with NaN values in numpy arrays
* [astropy](http://www.astropy.org/) - Python tools for astronomy
* [scikit-learn](http://scikit-learn.org/) - Machine learning tools.  Used for PCA in ZAPv2.

as well as numpy and scipy.



# Contents
* `scripts`: the main point of entry to using pySPARK.  This currently includes a relatively generic example script which can be driven using the example configuration file found in `examples', e.g.

```
reduce_kmos.py -p configuration.cfg
```
Currently this will execute all steps necesary to produce sky-subtracted, illumination correct data cubes.  Combining is not currently implemented, but will be available as an option soon.

* `pyspark`: core processing routines used by the reduction script.

* `doc`: contains (eventually) a manual, installation instructions, and version history, blah blah blah.
