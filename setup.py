from distutils.command import clean
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration


def configuration(parent_package='', top_path='', package_path='halftones'):
    config = Configuration('halftones', parent_package, top_path, package_path)
    config.add_extension('compiled',
                          sources = ['halftones/halftones_compiled.f90'])
    config.add_include_dirs('halftones')
    return config

if __name__ == "__main__":
    setup(**configuration().todict())
