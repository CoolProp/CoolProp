"""
CPWeb - A collection of commonly used routines to produce CoolProp's online documentation
=====

"""
from __future__ import division, absolute_import, print_function

import codecs
import csv
import cStringIO


def get_version():
    return 5.0


if __name__ == "__main__":
    print('You are using version %s of the Python package for creating CoolProp\' online documentation.'%(get_version()))
    print()
