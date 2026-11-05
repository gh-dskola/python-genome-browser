"""pygbrowse: plots of genomic datasets in a genome-browser-like format."""

from importlib.metadata import PackageNotFoundError, version

from . import datasources, plots, utilities

try:
    __version__ = version("pygbrowse")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"

name = "pygbrowse"
