#!/usr/bin/env python

import string

class align_block:
    def __init__(self, indel_count, length, mismatches):
        self.indel_count = int(indel_count)
        self.length = int(length)
        self.mismatches = int(mismatches)
    def __str__(self):
        return string.join( map( str, [self.indel_count, self.length, self.mismatches]), ' ')

class qltout_record:

    def __init__(self, obj):
        self.fields = []
        self.align_blocks = []
        if type(obj) == str:
            self.setValuesFromString( obj );
        else:
            raise TypeError("qltout_record did not recognize type: "+str(type(obj)))
    
    def setValuesFromString( self, line ):
        formats = [str, int, int, int, int, int, int, int, int, int, int]
        
        split_line = line.split()
        self.fields = map( lambda f, v: f(v), formats, split_line[:11] )

        next_align_block = 11
        while next_align_block < len(split_line):
            self.align_blocks.append( align_block(*split_line[next_align_block:next_align_block+3]) )
            next_align_block += 3
            #print self.__str__()

    def indelled_bases(self):
        return sum([ab.indel_count for ab in self.align_blocks])

    def __str__(self):
        return string.join( map( str, self.fields ), ' ')+' '+string.join( map( str, self.align_blocks ), ' ')

    def id(self): return self.fields[1]
    def contig(self): return self.fields[6]

    def offset(self): return self.fields[7]
    def pos(self): return self.fields[7]+1

    def offset_end(self): return self.fields[8]
    def pos_end(self): return self.fields[8]+1

    def linear_start(self): return self.fields[9]
    def map_qual(self): return 100 # This only exists in maq files
    # for ILT and Merlin, we assume the same confidence for all mappings


class qltout_file:
    def __init__(self, filename):
        self.filename = filename
        self.fqlt = open(self.filename)

    def RecordGenerator(self):
        
        for line in self.fqlt:
            #fields = line.rstrip("\n").split("\t")
            if not line.startswith("QUERY"):
                continue
            yield qltout_record(line) #fields)
        raise StopIteration
                        
    def __iter__(self):
        return self.RecordGenerator()
          
            
