
#!/usr/bin/env python3
#-*- coding: utf-8 -*-

### Adapted from gist of Filip Dominec:
### https://gist.github.com/FilipDominec/9ff081952dbc4aae1df657a56c3db4ea

### Usage: python journal_abbreviation.py filename.bib


import sys, os
import re

try:    bibtexdb = open(sys.argv[1]).read()
except: print("Error: specify the file to be processed!")

if not os.path.isfile('journalList.txt'):
    import urllib.request
    urllib.request.urlretrieve("https://gist.githubusercontent.com/FilipDominec/6df14b3424e335c4a47a96640f7f0df9/raw/74876d2d5df9ed60492ef3a14dc3599a6a6a9cfc/journalList.txt", 
            filename="journalList.txt")
rulesfile = open('journalList.txt')

for rule in rulesfile.readlines()[::-1]:           ## reversed alphabetical order matches extended journal names first
    pattern1, pattern2 = rule.strip().split(" = ")
    if pattern1 != pattern1.upper() and (' ' in pattern1):        ## avoid mere abbreviations
    #bibtexdb = bibtexdb.replace(pattern1.strip(), pattern2.strip())    ## problem - this is case sensitive
        repl = re.compile(re.escape(pattern1), re.IGNORECASE)               ## this is more robust, although ca. 10x slower
        (bibtexdb, num_subs) = repl.subn(pattern2, bibtexdb)
        if num_subs > 0:
            print("Replacing '%s' FOR '%s'" % (pattern1, pattern2))

with open('abbreviated.bib', 'w') as outfile:
    outfile.write(bibtexdb)
print("Bibtex database with abbreviated files saved into 'abbreviated.bib'")
