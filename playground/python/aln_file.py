#!/usr/bin/env python

import string, sys

class aln_record:
    "Stores one record of data from a MAQ .aln.txt file"
    field_names = (
        "read name",
        "chromosome",
        "position",
        "strand",
        "insert size from the outer coorniates of a pair",
        "paired flag",
        "mapping quality",
        "single-end mapping quality",
        "alternative mapping quality",
        "number of mismatches of the best hit",
        "sum of qualities of mismatched bases of the best hit",
        "number of 0-mismatch hits of the first 24bp",
        "number of 1-mismatch hits of the first 24bp on the reference",
        "length of the read",
        "read sequence",
        "sequence quality"
        )
    max_field_name_len = max(map(len, field_names))

    def __init__(self, obj, parse = True):
        self.fields = []
        if type(obj) == str:
            self.setValuesFromString( obj, parse );
        else:
            raise TypeError("aln_record did not recognize type: "+str(type(obj)))
    
    def setValuesFromString( self, line, parse = True ):
        if parse:
            formats = [str, str, int, str, int, int, int, int, int, int, int, int, int, int, str, str]
            self.fields = map( lambda f, v: f(v), formats, line.split() )
        else:
            self.fields = line.split()

    def __str__(self):
        s = ""
        for n,v in zip(aln_record.field_names, self.fields):
            s += ("%"+str(aln_record.max_field_name_len)+"s : %s\n") % (n, str(v))
        return s

        #return string.join( map( str, self.fields ), ' ')

    def id(self): return self.fields[0]
    def contig(self): return self.fields[1]

    def offset(self): return self.fields[2]-1
    def pos(self): return self.fields[2]

    # Quality of read mapping (only maq gives this field)
    def map_qual(self): return self.fields[6]

    #def offset_end(self): return self.fields[8]
    #def pos_end(self): return self.fields[8]+1

    #def linear_start(self): return self.fields[9]


class aln_file:
    def __init__(self, filename, parse=True):
        self.filename = filename
        self.parse = parse
        self.faln = open(self.filename)

    def RecordGenerator(self):
        for line in self.faln:
            yield aln_record(line, self.parse)
        raise StopIteration
                        
    def __iter__(self):
        return self.RecordGenerator()
          
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "To test aln_file class:\naln_file.py ALN_FILE"
    else:
        count = 0
        for aln in aln_file(sys.argv[1]):
            print aln
            count += 1
            #if count > 5:
            #    break
