from importlib.metadata import PackageNotFoundError, version

from .core import Barcodes

try:
    __version__ = version('barcodeEZ')
except PackageNotFoundError:      # running from a source tree, not installed
    __version__ = 'unknown'

__all__ = ['Barcodes', '__version__']
