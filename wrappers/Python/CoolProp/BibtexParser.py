#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import division, absolute_import, print_function
from __future__ import generators
from pybtex.database.input import bibtex
from pybtex.plugin import find_plugin
import io
import latexcodec.codec
import copy
import string
from pybtex.database import FieldDict


class BibTeXerClass(object):
    """
    A base class that defines all the variables needed
    to print a very basic bibliography.
    """

    def __init__(self, fName = u'../../../CoolPropBibTeXLibrary.bib'):
        self.DEBUG = False

        latexcodec.codec.register()

        self.sanitizeDict = {}
        self.sanitizeDict['year']        = "no year"
        self.sanitizeDict['title']       = "no title"
        self.sanitizeDict['journal']     = "no journal"
        self.sanitizeDict['volume']      = "no volume"
        self.sanitizeDict['pages']       = "no pages"
        self.sanitizeDict['booktitle']   = "no book title"
        self.sanitizeDict['school']      = "no school"
        self.sanitizeDict['note']        = "no note"
        self.sanitizeDict['publisher']   = "no publisher"
        self.sanitizeDict['institution'] = "no institution"
        self.loadLibrary(fName)


    def loadLibrary(self, path, keys=[]):
        if len(keys)>0:
            bib_parser = bibtex.Parser(wanted_entries=keys)
        else:
            bib_parser = bibtex.Parser()

        self.library = bib_parser.parse_file(path)

        #filename = path.splitext(aux_filename)[0]
        #aux_data = auxfile.parse_file(aux_filename, output_encoding)
        #bib_parser = find_plugin('pybtex.database.input', "bibtex")
        #bib_data = bib_parser(
        #    encoding=bib_encoding,
        #    wanted_entries=aux_data.citations,
        #    min_crossrefs=min_crossrefs,
        #).parse_files(aux_data.data, bib_parser.default_suffix)



    ########################################
    # Find entry in library and make unicode
    ########################################
    def getUnicodeEntry(self, key):
        if self.library is None:
            raise ValueError("Library has not been loaded.")
        entry = copy.deepcopy(self.library.entries[key])

        for key, value in entry.fields.iteritems():
            entry.fields[key] = self.fromLatex(value)

        for key in entry.persons.keys():
            for i in range(len(entry.persons[key])):
                entry.persons[key][i]._first = self.fromLatex(entry.persons[key][i].first())
                entry.persons[key][i]._middle = self.fromLatex(entry.persons[key][i].middle())
                entry.persons[key][i]._prelast = self.fromLatex(entry.persons[key][i].prelast())
                entry.persons[key][i]._last = self.fromLatex(entry.persons[key][i].last())
                entry.persons[key][i]._lineage = self.fromLatex(entry.persons[key][i].lineage())

        return entry


    ########################################
    # Basic text processing functions
    ########################################

    def stripCurls(self, text):
        table = { ord(u'{'): None, ord(u'}'): None }
        stripped = text.translate(table)
        if self.DEBUG: print(u"From {0} to {1}".format(text, stripped))
        return stripped


    def fromLatex(self, latexStr):
        """Converter for Latex strings

        A function that always returns unicode. It also tries to recurse into
        lists and dicts, but be careful this is not tested.
        """
        if isinstance(latexStr, str): # ordinary string
            latexStr = unicode(latexStr)
        elif isinstance(latexStr, unicode): # unicode string
            pass
        else:
            try:
                for k in latexStr.keys():
                    latexStr[k] = self.fromLatex(latexStr[k])
            except AttributeError:
                pass
            # Check for list
            if type(latexStr)==type([]):
                for i in range(len(latexStr)):
                    latexStr[i] = self.fromLatex(latexStr[i])
            # return the plain object
            return latexStr

        if self.DEBUG: print(u"From {0} ".format(latexStr),end=u"")
        decoded = latexStr.decode('latex')
        if self.DEBUG: print(u"to {0} ".format(decoded))
        return decoded


    #######################################
    # Custom output routines
    #######################################

    def prepareSingleEntry(self, key):

        # TODO: What is this for?
        if key.startswith('__'):
            return None,None

        entry = self.getUnicodeEntry(key)
        authors = u'; '.join( unicode(author) for author in entry.persons['author'] )
        authors = self.stripCurls(authors)

        # Strip off the opening and closing brackets
        for key in entry.fields:
            if entry.fields[key].startswith('{') and entry.fields[key].endswith('}'):
                entry.fields[key] = entry.fields[key][1:len(entry.fields[key])-1]

        # Use the basic dictionary and update it with existing values
        #entry.fields = FieldDict(self.sanitizeDict).update(entry.fields)
        for key in self.sanitizeDict:
            if not key in entry.fields:
                entry.fields[key] = self.sanitizeDict[key]

        return authors,entry


    def entry2txt(self, key):

        authors,entry = self.prepareSingleEntry(key)
        if authors is None:
            return None
        else:
            year        = entry.fields['year']
            title       = entry.fields['title']
            journal     = entry.fields['journal']
            volume      = entry.fields['volume']
            pages       = entry.fields['pages']
            booktitle   = entry.fields['booktitle']
            school      = entry.fields['school']
            note        = entry.fields['note']
            publisher   = entry.fields['publisher']
            institution = entry.fields['institution']

            if entry.type == 'article':
                return authors + ', ' + year + ', ' + title + ', ' + journal + ', ' + volume + ':' + pages
            elif entry.type == 'conference':
                return authors + ', ' + year + ', ' + title + ', ' + booktitle
            elif entry.type == 'mastersthesis' or entry.type == 'phdthesis':
                return authors + ', ' + year + ', ' + title + ', ' + school
            elif entry.type == 'unpublished':
                return authors + ', ' + year + ', ' + title + ', note: ' + note
            elif entry.type == 'book':
                return authors + ', ' + year + ', ' + title + ', ' + publisher
            elif entry.type == 'techreport':
                return authors + ', ' + year + ', ' + title + ', ' + institution
            else:
                print(u"There was an error processing: {0}".format(entry))
                return None


    def entry2rst(self, key):

        authors,entry = self.prepareSingleEntry(key)
        if authors is None:
            return None
        else:
            year        = entry.fields['year']
            title       = entry.fields['title']
            journal     = entry.fields['journal']
            volume      = entry.fields['volume']
            pages       = entry.fields['pages']
            booktitle   = entry.fields['booktitle']
            school      = entry.fields['school']
            note        = entry.fields['note']
            publisher   = entry.fields['publisher']
            institution = entry.fields['institution']

            if entry.type == 'article':
                return authors + ', ' + year + ', ' + title + ', *' + journal + '*, ' + volume + ':' + pages
            elif entry.type == 'conference':
                return authors + ', ' + year + ', ' + title + ', *' + booktitle + '*'
            elif entry.type == 'mastersthesis' or entry.type == 'phdthesis':
                return authors + ', ' + year + ', ' + title + ', *' + school + '*'
            elif entry.type == 'unpublished':
                return authors + ', ' + year + ', ' + title + ', note: ' + note
            elif entry.type == 'book':
                return authors + ', ' + year + ', *' + title + '*, ' + publisher
            elif entry.type == 'techreport':
                return authors + ', ' + year + ', *' + title + '*, ' + institution
            else:
                print(u"There was an error processing: {0}".format(entry))
                return None


    def entry2html(self, key):

        authors,entry = self.prepareSingleEntry(key)
        if authors is None:
            return None
        else:
            year        = entry.fields['year']
            title       = entry.fields['title']
            journal     = entry.fields['journal']
            volume      = entry.fields['volume']
            pages       = entry.fields['pages']
            booktitle   = entry.fields['booktitle']
            school      = entry.fields['school']
            note        = entry.fields['note']
            publisher   = entry.fields['publisher']
            institution = entry.fields['institution']

            if entry.type == 'article':
                return authors + ', ' + year + ', ' + title + ', <i>' + journal + '</i>, ' + volume + ':' + pages
            elif entry.type == 'conference':
                return authors + ', ' + year + ', ' + title + ', <i>' + booktitle + '</i>'
            elif entry.type == 'mastersthesis' or entry.type == 'phdthesis':
                return authors + ', ' + year + ', ' + title + ', <i>' + school + '</i>'
            elif entry.type == 'unpublished':
                return authors + ', ' + year + ', ' + title + ', note: ' + note
            elif entry.type == 'book':
                return authors + ', ' + year + ', <i>' + title + '</i>, ' + publisher
            elif entry.type == 'techreport':
                return authors + ', ' + year + ', <i>' + title + '</i>, ' + institution
            else:
                print(u"There was an error processing: {0}".format(entry))
                return None


    def entries2All(self,

        keys=None,
        output_encoding=None,
        output_backend='plaintext',
        style='unsrt'):

        style_cls = find_plugin('pybtex.style.formatting', style)
        style = style_cls()
        bib_data = self.library #[ self.getEntry(key) for key in keys ]
        formatted_bibliography = style.format_bibliography(bib_data,citations=keys)

        output_backend = find_plugin('pybtex.backends', output_backend)
        #output_filename = filename + output_backend.default_suffix
        #output_backend(output_encoding).write_to_file(formatted_bibliography, output_filename)

        stream = io.StringIO()
        output_backend(output_encoding).write_to_stream(formatted_bibliography, stream)
        # Retrieve contents
        contents = stream.getvalue()
        # Close object and discard memory buffer --
        # .getvalue() will now raise an exception.
        stream.close()
        return contents



if __name__=='__main__':
#     B = BibTeXerClass()
#     B.loadLibrary('../../../Web/fluid_properties/Incompressibles.bib')
#     print(B.entry2rst('Cesar2013'))
#     print(B.entry2html('Cesar2013'))
#     print(B.entries2All(keys=['Cesar2013']))
    B = BibTeXerClass('../../../CoolPropBibTeXLibrary.bib')
    print(B.entry2rst('Mulero-JPCRD-2012'))
    print(B.entry2html('Mulero-JPCRD-2012'))
    print(B.entries2All(keys=['Mulero-JPCRD-2012']))
