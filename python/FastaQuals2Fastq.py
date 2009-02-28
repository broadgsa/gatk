#!/usr/bin/env python

# Utility script to convert FASTA + QUALS file into FASTQ required for MAQ
# MAQ needs to convert fastq into bfq file using fastq2bfq in order to use
# this data

import sys, fasta, os
from itertools import *

def usage():
    print "Usage: fastq.py fasta_in quals_in fastq_out"
    print "  fasta_in: *.fasta file"
    print "  quals_in: *.quals.txt OR *.quala"
    sys.exit()

def output_fastq_record(file_handle, fasta, quals):
    #print 'fasta, quals', fasta, quals
    #for fasta, qual in zip(fastas, quals):
    
    qualsVector = map(int, quals.seq.split())
    
    print >> fqout, "@"+fasta.id
    print >> fqout, fasta.seq
    print >> fqout, "+"
    print >> fqout, "".join( map(lambda q: chr(q+33), qualsVector ) )

def qualStr2fasta( line ):
    elts = line.split()
    id = elts[0]
    quals = ' '.join(elts[1:])
    return fasta.fasta_record(id, quals)

if __name__ == "__main__":
    if len(sys.argv[1:]) != 3:
        usage()
    else:
        fasta_in_file, quals_in_file, fastq_out_file = sys.argv[1:]

    print "Quals input file: "+quals_in_file
    print "Fasta input file: "+fasta_in_file
    
    if os.path.splitext(quals_in_file)[1] == ".txt":
        qual_factory = imap( qualStr2fasta, file(quals_in_file) )
    else:
        if os.path.splitext(quals_in_file)[1] == ".quala":
            # Create fasta factory from the quala file and treat it 
            # as if it was a fasta file and pull the numbers
            qual_factory = fasta.fasta_file(quals_in_file, cleanup=False)
        else:
            print "Qual input file should be a *.quals.txt or *.quala file"

    print "Fastq ouptut file: "+fastq_out_file
    fqout = file(fastq_out_file, "w")

    for qual_record, fasta_record in izip_longest(qual_factory, fasta.fasta_file(fasta_in_file)):  
        if qual_record == None or fasta_record == None:
            # We ran out of one record type, so warn the user!
            print "WARNING: Different number of quals and fastas in input files!"
            sys.exit(1)
        #print 'qual_record', qual_record
        #print 'fasta_record', fasta_record
        output_fastq_record(fqout, fasta_record, qual_record)


    
