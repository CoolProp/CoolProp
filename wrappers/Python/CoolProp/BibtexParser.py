#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import division, absolute_import, print_function
from pybtex.database.input import bibtex
from pybtex.plugin import find_plugin
import io



class BibTeXerClass(object):
    """
    A base class that defines all the variables needed
    to print a very basic bibliography.
    """

    def __init__(self, fName = '../../../CoolPropBibTeXLibrary.bib'):
        self.DEBUG = False

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



    def getEntry(self, key):
        if self.library is None:
            raise ValueError("Library has not been loaded.")
        if key in self.library.entries.keys():
            return self.library.entries[key]
        else:
            return "Key {0} not in library.".format(key)

    #######################################
    # Some of Ian's old functions
    #######################################
    def accent_substitutions(self, name):
        mapping = [('{\\~n}','\xf1'), # �
                   ('{\\`e}','\xe8'), # �
                   ("{\\'e}",'\xe9'), # �
                   ("{\\'a}",'\xe1'), # �
                   ("{\\`a}",'\xe0'), # �
                   ("{\\'i}",'\xed'), # �
                   ("{\\'i}",'\xed'), # �
                   ('{\\\"o}','\xf6'), # �
                   ('{\\\"u}','\xfc'), # �
                   ('{\\v s}','\x161'), # �
                   ]
        for old, new in mapping:
            name = name.replace(old, new)
        return name

    def count_substr(self, s, ss):
        c = 0
        for e in s:
            if e == ss:
                c += 1
        return c

    def DE(self, s):
        try:
            return s.decode('ascii').encode('utf-8')
        except UnicodeEncodeError:
            print('Decoding error for',s)


    def entry2rst(self, key):

        if key.startswith('__'):
            return ''

        entry = self.getEntry(key)

        if entry is None:
            return ''

        try:
            authors = '; '.join([self.accent_substitutions(unicode(author).decode('ascii').encode('utf-8')) for author in entry.persons['author']])
        except UnicodeEncodeError:
            print('Decoding error for',[author for author in entry.persons['author']])

        if authors.find('{') > -1 or authors.find('}') > -1:
            print(authors)
            raise ValueError("authors [{authors:s}] may not have '{{' or '}}' character".format(authors = authors))

        fields = entry.fields

        # Strip off the opening and closing brackets
        fields['title'] = fields['title'].strip()
        if fields['title'].startswith('{') and fields['title'].endswith('}'):
            fields['title'] = fields['title'][1:len(entry.fields['title'])-1]

        f = fields
        for key in f:
            f[key] = self.DE(f[key])
        authors = str(authors)

        if entry.type == 'article':
            if 'journal' not in fields: fields['journal'] = ''
            if 'volume' not in fields: fields['volume'] = ''
            if 'pages' not in fields: fields['pages'] = ''

            return authors + ', ' + f['year'] + ', ' + f['title'] + ', *' + f['journal'] + '*, ' + f['volume'] + ':' + f['pages']

        elif entry.type == 'conference':
            if 'journal' not in f: f['journal'] = ''
            return authors + ', ' + f['year'] + ', ' + f['title'] + ', *' + f['booktitle'] + '*'

        elif entry.type == 'mastersthesis':
            return authors + ', ' + f['year'] + ', ' + f['title'] + ', *' + f['school'] + '*'

        elif entry.type == 'unpublished':
            return authors + ', ' + f['year'] + ', ' + f['title'] + ', note: ' + f['note']

        elif entry.type == 'book':
            return authors + ', ' + f['year'] + ', *' + f['title'] + '*, ' + f['publisher']

        elif entry.type == 'techreport':
            return authors + ', ' + f['year'] + ', *' + f['title'] + '*, ' + f['institution']

        else:
            print(entry)

    def entry2HTML(self, key):

        if key.startswith('__'):
            return ''

        entry = self.getEntry(key)

        if entry is None:
            return ''

        try:
            authors = '; '.join([self.accent_substitutions(unicode(author).decode('ascii').encode('utf-8')) for author in entry.persons['author']])
        except UnicodeEncodeError:
            print('Decoding error for',[author for author in entry.persons['author']])

        if authors.find('{') > -1 or authors.find('}') > -1:
            print(authors)
            raise ValueError("authors [{authors:s}] may not have '{{' or '}}' character".format(authors = authors))

        fields = entry.fields

        # Strip off the opening and closing brackets
        fields['title'] = fields['title'].strip()
        if fields['title'].startswith('{') and fields['title'].endswith('}'):
            fields['title'] = fields['title'][1:len(entry.fields['title'])-1]

        f = fields
        for key in f:
            f[key] = self.DE(f[key])
        authors = str(authors)

        if entry.type == 'article':
            if 'journal' not in fields: fields['journal'] = ''
            if 'volume' not in fields: fields['volume'] = ''
            if 'pages' not in fields: fields['pages'] = ''

            return authors + ', ' + f['year'] + ', ' + f['title'] + ', <i>' + f['journal'] + '</i>, ' + f['volume'] + ':' + f['pages']

        elif entry.type == 'conference':
            if 'journal' not in f: f['journal'] = ''
            return authors + ', ' + f['year'] + ', ' + f['title'] + ', <i>' + f['booktitle'] + '</i>'

        elif entry.type == 'mastersthesis':
            return authors + ', ' + f['year'] + ', ' + f['title'] + ', <i>' + f['school'] + '</i>'

        elif entry.type == 'unpublished':
            return authors + ', ' + f['year'] + ', ' + f['title'] + ', note: ' + f['note']

        elif entry.type == 'book':
            return authors + ', ' + f['year'] + ', <i>' + f['title'] + '</i>, ' + f['publisher']

        elif entry.type == 'techreport':
            return authors + ', ' + f['year'] + ', <i>' + f['title'] + '</i>, ' + f['institution']

        else:
            print(entry)


    def entries2All(self,
        keys=[],
        output_encoding=None,
        output_backend='plaintext',
        style='unsrt'):

        style_cls = find_plugin('pybtex.style.formatting', style)
        style = style_cls()
        bib_data = self.library #[ self.getEntry(key) for key in keys ]
        formatted_bibliography = style.format_bibliography(bib_data)

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
        print(contents)

if __name__=='__main__':
    B = BibTeXerClass()
    B.loadLibrary('../../../Web/fluid_properties/Incompressibles.bib')
    print(B.entry2rst('Cesar2013'))
    print(B.entries2All(keys=['Cesar2013']))
#     B = BibTeXerClass(fName='../../../Web/CoolPropBibTeXLibrary.bib')
#     print(B.entry2rst('Mulero-JPCRD-2012'))
#     print(B.entries2All(keys=['Mulero-JPCRD-2012']))