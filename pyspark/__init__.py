try:
    from ._version import __version__
except(ImportError):
    pass

from . import detector
from . import utils
from . import calibration
from . import cubes
from . import sky


__all__ = ['detector', 'utils', 'calibration', 'cubes', 'sky']
