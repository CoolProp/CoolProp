#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import division, absolute_import, print_function
from pybtex.database.input import bibtex



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


    def loadLibrary(self, path):
        parser = bibtex.Parser()
        self.library = parser.parse_file(path)
        self.entries = self.library.entries


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

        entry = self.entries[key]

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

        entry = self.entries[key]

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

    def findentry(self, key):
        for entry in self.entries:
            if entry['key'] == key:
                return entry

if __name__=='__main__':
    B = BibTeXerClass()
    print(B.entry2rst('Mulero-JPCRD-2012'))