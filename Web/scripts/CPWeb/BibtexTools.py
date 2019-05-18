'''
Created on 26 Sep 2014

@author: jowr
'''
import os
from CoolProp.BibtexParser import BibTeXerClass


def getPath(filename, search=True):
    # Path to root
    coolprop_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    # Test for file
    fname = filename
    if os.path.exists(fname): return os.path.abspath(fname)
    # Test relative to this file
    fname = os.path.join(os.path.dirname(__file__), filename)
    if os.path.exists(fname): return os.path.abspath(fname)
    # Test relative to root notation
    fname = os.path.join(coolprop_dir, filename)
    if os.path.exists(fname): return os.path.abspath(fname)
    # Search in root tree
    fname = os.path.basename(filename)
    if search:
        result = []
        for root, dirs, files in os.walk(coolprop_dir):
            if fname in files:
                result.append(os.path.join(root, fname))
        if len(result) == 1:
            return os.path.abspath(result[0])
        elif len(result) > 1:
            print("Found multiple files with the name {0}. Try to specify the path as well.".format(fname))
            print(result)
            return os.path.abspath(result[0])

    raise ValueError("Found no file with the name {0}. Try to specify the path as well.".format(fname))


def getBibtexParser(filename='../../../CoolPropBibTeXLibrary.bib'):
    """Create a parser object that can be used to extract entries from a
    library in Bibtex format."""
    fpath = getPath(filename)
    bibtexer = BibTeXerClass(fpath)
    return bibtexer


def getCitationOrAlternative(bibtexer, bibtex_key):
    """Find the key in the library and convert to a citation, if it is not found,
    we return a footnote string for sphinx."""
    bibtex_key = bibtex_key.strip()

    if bibtex_key in bibtexer.library.entries:
        return u':cite:`{0}`'.format(bibtex_key)
    else:
        return u':raw-html:`<span title="{0}">Source</span>`'.format(bibtex_key)
