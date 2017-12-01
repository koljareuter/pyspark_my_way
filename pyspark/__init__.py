try:
    from ._version import __version__
except(ImportError):
    pass

from . import detector
from . import utils
from . import calibration
from . import cubes
from . import sky
