#!/usr/bin/env python

import sys, itertools

def record_generator(filename, sep="\t", skip_n_lines=0, skip_until_regex_line=""):
    """Given a file with field headers on the first line and records on subsequent lines,
generates a dictionary for each line keyed by the header fields"""
    fin = open(filename)

    if skip_n_lines > 0:
        for i in range(skip_n_lines): # Skip a number of lines
            fin.readline()

    found_regex = False
    if skip_until_regex_line != "":
        import re
        regex_line = re.compile(skip_until_regex_line)
        for line in fin:
            match = regex_line.search(line)
            if match:
                found_regex = line
                break
        if not found_regex:
            print "Warning: Regex "+skip_until_regex_line+" not found in FlatFileTable:record_generator"

    if found_regex:
        header = found_regex.rstrip().split(sep) # Parse header
    else:
        header = fin.readline().rstrip().split(sep) # Pull off header
    
    for line in fin: # 
        fields = line.rstrip().split(sep)
        record = dict(itertools.izip(header, fields))
        yield record

def record_matches_values(record, match_field_values):
    for match_field, match_values in match_field_values:
        if record[match_field] not in match_values:
            return False
    return True
