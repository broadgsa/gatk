#!/usr/bin/env python

import sys, itertools, FlatFileTable

def subset_list_by_indices(indices, list):
    subset = []
    for index in indices:
        subset.append(list[index])
    return subset

def chunk_generator(record_gen, key_fields):
    """Input:
  line_gen: generator that produces lines with linefeeds chopped off
  key_fields: field numbers in each record used to determine chunk membership
Output:
  locus_chunk: list of consecutive lines that have the same key_fields"""
  
    locus_chunk = []
    last_key = ""
    first_record = True
    for record in record_gen:
        key = [record[f] for f in key_fields]
        if key == last_key or first_record:
            locus_chunk.append(record)
            first_record = False
        else:
            if locus_chunk != []:
                yield locus_chunk
                locus_chunk = [record]
        last_key = key
    yield locus_chunk

def chunk_stats(chunk):
    records = 0
    conf_calls = 0
    correct_genotype = 0
    for record in chunk:
        if abs(float(record["BtnbLod"])) >= 5:
            conf_calls += 1
            if record["HapmapChipGenotype"] == record["BestGenotype"]:
                correct_genotype += 1
        records += 1

    return float(correct_genotype) / max(conf_calls,1), float(conf_calls) / max(records,1)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: CoverageEval.py geli_file")
    filename = sys.argv[1]

    locus_gen = chunk_generator(FlatFileTable.record_generator(filename, None), ("Sequence","Position"))
    print "Fraction correct genotype\tCoverage sampled\tLocus\tReference base\tHapmap chip genotype (Max. coverage genotype call for reference calls)"
    for locus in locus_gen:
        #print "NEW LOCUS"
        covs = dict()
        coverage_chunk_gen = chunk_generator(locus, ("DownsampledCoverage", "Sequence", "Position"))
        for cov_chunk in coverage_chunk_gen:
            #print "NEW COVERAGE"
            #print "\n".join(cov_chunk)
            record = cov_chunk[0]
            print "\t".join(map(str,("%.4f\t%.4f"%chunk_stats(cov_chunk), record["DownsampledCoverage"], record["Sequence"]+":"+record["Position"],record["ReferenceBase"],record["HapmapChipGenotype"])))
            
            
