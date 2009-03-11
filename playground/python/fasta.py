#!/usr/bin/env python

import string

class fasta_record:
    "Record containing one FASTA sequence"
    def __init__(self, id, seq):
        self.id = id
        self.seq = seq
    
    def __str__(self):
        return '['+self.id+" "+self.seq+']'

class fasta_file:
    "Iterable object based on FASTA file format"
    def __init__(self, filename, cleanup=True):
        "cleanup removes spaces from fasta text (default: True)"
        self.filename = filename
        self.fin = open(self.filename)
        self._cleanup = cleanup

    def RecordGenerator(self):
        line = self.fin.readline().rstrip()
        assert line[0] == ">"
        id = line[1:]
        seq = ""
        for line in self.fin:
            line = line.rstrip()
            if line[0] == ">":
                yield fasta_record(id, seq)
                id = line[1:]
                seq = ""
            else:
                if self._cleanup:
                    seq += line.replace(" ","")
                else:
                    seq += line

        yield fasta_record(id, seq) # Yield last seq
        raise StopIteration # No more lines
                        
    def __iter__(self):
        return self.RecordGenerator()
          

if __name__ == "__main__":
    print "Testing fast.py on file 5seqs.fa..."
    for fasta_rec in fasta_file("5seqs.fa"):
         print fasta_rec

