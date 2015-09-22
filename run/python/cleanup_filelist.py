#!/bin/env python

# loop over a filelist, and identify the files with empty trees
#
# davide.gerbaudo@gmail.com
# September 2015

import argparse
import sys
import ROOT as R
R.gROOT.SetBatch(1)

usage = """
cat filelist.txt | %(prog)s [options] 2>/dev/null | tee output.txt

Note that you need to throw away the stderr produced by root.
"""
def main():
    parser = argparse.ArgumentParser(description='comment out files with empty trees.')
    parser.add_argument('-i', '--input', type = argparse.FileType('r'), default = '-')
    parser.add_argument('-t', '--tree-name', default='susyNt')
    parser.add_argument('-c', '--comment-char', default='#')
    parser.add_argument('-d', '--drop-lines', action='store_true', help='drop empty files rather than commenting them out')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    number_of_affected_lines = 0
    comment_char = args.comment_char
    drop_lines = args.drop_lines
    for line in args.input:
        empty_tree = has_empty_tree(line.strip(), 'susyNt')
        if empty_tree and drop_lines:
            pass
        else:
            sys.stdout.write((comment_char+' 0 entries ' if empty_tree else '')+line)
        number_of_affected_lines += 1 if empty_tree else 0
    if args.verbose:
        print comment_char+" Number of affected lines: %d" % number_of_affected_lines

def has_empty_tree(filename='', treename=''):
    empty_tree = False
    input_file = R.TFile.Open(filename) if filename else None
    if input_file:
        input_tree = input_file.Get(treename)
        if input_tree:
            empty_tree = input_tree.GetEntries()==0
        input_file.Close()
        input_file.Delete()
    return empty_tree

if __name__ == "__main__":
    main()
