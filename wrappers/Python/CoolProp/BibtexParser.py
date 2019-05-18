#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import division, absolute_import, print_function
from __future__ import generators
import pybtex.plugin, pybtex.database.input.bibtex, pybtex.errors
import io
import codecs, latexcodec
import os
import six

# Here we are going to simply hack the formatting of the article class such that it preserves the
# desired spacing and capitalization.  Probably this problem could be better solved with new styles,
# but this seems to work...
import pybtex.style.formatting as formatting
from pybtex.style.template import node, join
# Create our new function


@node
def toplevel(children, data):
    return join(sep=' ')[children].format_data(data)


# And now, we over-write the function with our desired function.  Et voila! It works!
formatting.toplevel = toplevel


class BibTeXerClass(object):
    """
    A base class that defines all the variables needed
    to print a very basic bibliography.
    """

    def __init__(self, fName=u'../../../CoolPropBibTeXLibrary.bib'):
        self.loadLibrary(fName)

    def loadLibrary(self, path, keys=[], encoding="latex"):
        """Open a BibTex file and do some initial parsing.
        The path to the file has to be provided, otherwise an exception is raised.

        Optionally, you can provide an array of keys to reduce the amount of
        processed entries immediately.

        By default, the Bibtex file gets scanned with the artificial \"latex\"
        encoding, which translates Latex commands to their unicode equivalents.
        If you need Latex output, you can skip this step by passing another
        codec to the \"encoding\" parameter, for example \"ascii\" or \"utf-8\".
        """

        if path is None:
            raise ValueError("You have to provide a path to a bibtex file.")

        # Create the parser object
        if len(keys) > 0:
            bib_parser = pybtex.database.input.bibtex.Parser(wanted_entries=keys)
        else:
            bib_parser = pybtex.database.input.bibtex.Parser()

        # TODO: not needed anymore?
        oldLatexCodec = False
        if oldLatexCodec:
            # Do not print that many warnings
            pybtex.errors.set_strict_mode(enable=False)

            # TODO: Remove empty lines to keep Pybtex from choking
            with open(path, "r") as f:
                lines = f.readlines()
                cleaned = [l.strip() for l in lines if l.strip()]

            path = path + ".filtered.bib"
            with open(path, "w") as f:
                f.writelines('\n'.join(cleaned))

        # Open the file and convert it according to the encoding
        with codecs.open(path, encoding=encoding) as stream:
            self.library = bib_parser.parse_stream(stream)

        # os.remove(path)

        # Do some post-processing if encoding was latex
        if encoding == "latex":
            for tag in self.library.entries:
                entry = self.library.entries[tag]
                for key, value in six.iteritems(entry.fields):
                    entry.fields[key] = self.stripCurls(value)
                    if key == 'Title':
                        entry.fields[key] = u'{' + entry.fields[key] + '}'
                for key in entry.persons.keys():
                    for i in range(len(entry.persons[key])):
                        entry.persons[key][i]._first = self.stripCurls(entry.persons[key][i].first())
                        entry.persons[key][i]._middle = self.stripCurls(entry.persons[key][i].middle())
                        entry.persons[key][i]._prelast = self.stripCurls(entry.persons[key][i].prelast())
                        entry.persons[key][i]._last = self.stripCurls(entry.persons[key][i].last())
                        entry.persons[key][i]._lineage = self.stripCurls(entry.persons[key][i].lineage())

    def stripCurls(self, text):
        """Remove curly brackets from processed Latex code

        A function that always returns unicode. It also tries to recurse into
        lists and dicts, but be careful this is not thoroughly tested.
        """
        table = {ord(u'{'): None, ord(u'}'): None}

        if isinstance(text, str):  # ordinary string
            pass
        else:
            try:
                for k in text.keys():
                    text[k] = self.stripCurls(text[k])
            except AttributeError:
                pass
            # Check for list
            if isinstance(text, type([])):
                for i in range(len(text)):
                    text[i] = self.stripCurls(text[i])
            # return the plain object
            return text

        stripped = text.translate(table)
        return stripped

    def getBibliography(self, keys=None, fmt="plaintext", style="unsrtalpha", enc=None, objects=False):
        """This function creates a formatted bibliography according to the
        defined parameters using Pybtex.
        Specify your desired output format using the \"fmt\" parameter. Supported
        formats are defined by Pybtex and the parameter is only passed on to the
        plugin search function: latex, html, plaintext and markdown.
        The \"style\" parameter is handled similarly. At the moment, Pybtex partly
        supports the styles: alpha, plain, unsrt and unsrtalpha.
        Use the objects parameter to obtain the formatted bibliography and the
        backend renderer instead of the unicode object.
        """

        if self.library is None:
            raise ValueError("No library has been loaded, yet.")

        style_cls = pybtex.plugin.find_plugin('pybtex.style.formatting', style)
        style = style_cls()
        data = self.library
        biblio = style.format_bibliography(data, citations=keys)

        backend_cls = pybtex.plugin.find_plugin('pybtex.backends', fmt)
        backend = backend_cls(encoding=enc)

        if objects:
            return biblio, backend

        stream = io.StringIO()
        backend.write_to_stream(biblio, stream)
        contents = stream.getvalue()
        stream.close()
        return contents

    def getEntry(self, key, label=False, fmt="plaintext", style="unsrtalpha", enc=None):
        """If you only want a single entry from your bibliography and not the
        whole thing, this function is for you.
        It uses Pybtex, but removes the prologue and the epilogue and optionally
        deletes the label.
        For the other parameters, have a look at :func:`getBibliography`.
        """
        # TODO: What is this for?
        if key.startswith('__'):
            return None

        biblio, backend = self.getBibliography(keys=[key], fmt=fmt, style=style, enc=enc, objects=True)

        # Disable prologue and the epilogue
        backend.write_prologue = lambda: None
        backend.write_epilogue = lambda: None

        stream = io.StringIO()
        backend.write_to_stream(biblio, stream)
        contents = stream.getvalue()
        stream.close()

        if label:
            return contents

        label_table = {
            'latex': '}',
            'html': '</dt>',
            'markdown': ']',
            'plaintext': ']',
            }  # How do we find the end of the label?

        end_of_label = contents.index(label_table[fmt])
        contents = contents[end_of_label + len(label_table[fmt]):].strip()

        # {} are needed to preserve capitalization, but we don't want them in the final output, so remove them
        contents = contents.replace('{', '').replace('}', '')

        if fmt == "latex":
            contents = contents.replace(u"\\newblock ", "")
            contents = codecs.encode(contents, "latex")
        elif fmt == "html":
            contents = contents.replace(u"<dd>", "")
            contents = contents.replace(u"</dd>", "")

        contents = contents.replace(u"\n", "")

        return contents


# getEntry(self, key, label=False, fmt="markdown", style="unsrtalpha", enc=None):


if __name__ == '__main__':
    B = BibTeXerClass('../../../CoolPropBibTeXLibrary.bib')
    print("\nLatex:")
    print(B.getEntry(key='Mulero-JPCRD-2012', fmt='latex'))
    print("\nHTML:")
    print(B.getEntry(key='Mulero-JPCRD-2012', fmt='html'))
    print("\nMarkdown:")
    print(B.getEntry(key='Mulero-JPCRD-2012', fmt='markdown'))
    print("\nText:")
    print(B.getEntry(key='Mulero-JPCRD-2012', fmt='plaintext'))
