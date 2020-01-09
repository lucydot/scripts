#!/usr/bin/python

# usage: python compile_refs.py main_file_name

import subprocess, sys

commands = [
    ['pdflatex', sys.argv[1] + '.tex','-interaction=nonstopmode'],
    ['bibtex', sys.argv[1] + '.aux'],
    ['pdflatex', sys.argv[1] + '.tex','-interaction=nonstopmode'],
    ['pdflatex', sys.argv[1] + '.tex','-interaction=nonstopmode']
]

for c in commands:
    subprocess.call(c)
