# coding: utf-8

from pybtex.database.input import bibtex

def accent_substitutions(name):

    mapping = [('{\\~n}','\xf1'), # ñ
               ('{\\`e}','\xe8'), # è
               ("{\\'e}",'\xe9'), # é
               ("{\\'a}",'\xe1'), # á
               ("{\\`a}",'\xe0'), # à
               ("{\\'i}",'\xed'), # í
               ("{\\'i}",'\xed'), # í
               ('{\\\"o}','\xf6'), # ö
               ('{\\\"u}','\xfc'), # ü
               ('{\\v s}','\x161'), # š
               ]
    for old, new in mapping:
        name = name.replace(old, new)
    return name

def count_substr(s, ss):
    c = 0
    for e in s:
        if e == ss:
            c += 1
    return c

def DE(s):
    try:
        return s.decode('ascii').encode('utf-8')
    except UnicodeEncodeError:
        print 'Decoding error for',s
    
class BibTeXerClass:
    
    def __init__(self, fName = '../CoolProp/CoolPropBibTeXLibrary.bib'):
        parser = bibtex.Parser()
        bib_data = parser.parse_file('../CoolProp/CoolPropBibTeXLibrary.bib')
        self.entries = bib_data.entries
                
    def entry2rst(self, key):
        
        if key.startswith('__'):
            return ''
        
        entry = self.entries[key]
        
        if entry is None:
            return ''
            
        try:
            authors = '; '.join([accent_substitutions(unicode(author).decode('ascii').encode('utf-8')) for author in entry.persons['author']])
        except UnicodeEncodeError:
            print 'Decoding error for',[author for author in entry.persons['author']]
            
        if authors.find('{') > -1 or authors.find('}') > -1:
            print authors
            raise ValueError("authors [{authors:s}] may not have '{{' or '}}' character".format(authors = authors))
            
        fields = entry.fields
            
        # Strip off the opening and closing brackets
        fields['title'] = fields['title'].strip()
        if fields['title'].startswith('{') and fields['title'].endswith('}'):
            fields['title'] = fields['title'][1:len(entry.fields['title'])-1]
        
        f = fields
        for key in f:
            f[key] = DE(f[key])
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
            print entry
            
    def entry2HTML(self, key):
        
        if key.startswith('__'):
            return ''
        
        entry = self.entries[key]
        
        if entry is None:
            return ''
            
        try:
            authors = '; '.join([accent_substitutions(unicode(author).decode('ascii').encode('utf-8')) for author in entry.persons['author']])
        except UnicodeEncodeError:
            print 'Decoding error for',[author for author in entry.persons['author']]
            
        if authors.find('{') > -1 or authors.find('}') > -1:
            print authors
            raise ValueError("authors [{authors:s}] may not have '{{' or '}}' character".format(authors = authors))
            
        fields = entry.fields
            
        # Strip off the opening and closing brackets
        fields['title'] = fields['title'].strip()
        if fields['title'].startswith('{') and fields['title'].endswith('}'):
            fields['title'] = fields['title'][1:len(entry.fields['title'])-1]
        
        f = fields
        for key in f:
            f[key] = DE(f[key])
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
            print entry
        
    def findentry(self, key):
        for entry in self.entries:
            if entry['key'] == key:
                return entry

if __name__=='__main__':
    B = BibTeXerClass()
    print B.entry2rst('Mulero-JPCRD-2012')