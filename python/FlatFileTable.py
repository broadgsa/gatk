#!/usr/bin/env python

import sys, itertools

def record_generator(filename, sep="\t"):
    """Given a file with field headers on the first line and records on subsequent lines,
generates a dictionary for each line keyed by the header fields"""
    fin = open(filename)
    header = fin.readline().rstrip().split() # pull off header
    for line in fin:
        fields = line.rstrip().split(sep)
        record = dict(itertools.izip(header, fields))
        yield record

def record_matches_values(record, match_field_values):
    for match_field, match_values in match_field_values:
        if record[match_field] not in match_values:
            return False
    return True
