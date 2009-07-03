#!/usr/bin/env python

import sys

def chopped_line_generator(filename):
    fin = open(filename)
    for line in fin:
        line = line.rstrip()
        yield line

def subset_list_by_indices(indices, list):
    subset = []
    for index in indices:
        subset.append(list[index])
    return subset

def chunk_generator(line_gen, key_fields):
    """Input:
  line_gen: generator that produces lines with linefeeds chopped off
  key_fields: field numbers in each record used to determine chunk membership
Output:
  locus_chunk: list of consecutive lines that have the same key_fields
"""
  
    locus_chunk = []
    last_key = ""
    for line in line_gen:
        fields = line.split()
        key = subset_list_by_indices(key_fields, fields)
        if key == last_key:
            locus_chunk.append(line)
        else:
            last_key =key
            if locus_chunk != []:
                yield locus_chunk
                locus_chunk = []
    yield locus_chunk

def chunk_stats(chunk):
    records = 0
    correct_genotype = 0
    for record in chunk:
        fields = record.split()
        if fields[0] == fields[8]:
            correct_genotype += 1
        records += 1
    return float(correct_genotype) / records

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: CoverageEval.py geli_file")
    filename = sys.argv[1]

    fin = open(filename)
    locus_gen = chunk_generator(chopped_line_generator(filename), (3,4))
    print "Fraction correct genotype\tCoverage sampled\tLocus\tReference base\tHapmap chip genotype (Max. coverage genotype call for reference calls)"
    for locus in locus_gen:
        #print "NEW LOCUS"
        covs = dict()
        coverage_chunk_gen = chunk_generator(locus, (1,3,4))
        for cov_chunk in coverage_chunk_gen:
            #print "NEW COVERAGE"
            #print "\n".join(cov_chunk)
            fields = cov_chunk[0].split()
            coverage = fields[1]
            print "\t".join(map(str,("%.2f"%chunk_stats(cov_chunk), coverage, fields[3]+":"+fields[4],fields[5],fields[0])))
            
            #covs[coverage] = cov_chunk
            
            
