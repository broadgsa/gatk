import string
import operator
from pushback_file import pushback_file

#   Type   TAG   Req?   Default     Description
SAMHeaderFields = \
    [['HD', 'VN', True,  '0.1.0',   'File format version'], \
     ['HD', 'SO', False, None,      'Sort order. Valid values are: unsorted, queryname or coordinate'], \
     ['SQ', 'SN', True,  None,      'Sequence name. Unique among all sequence records in the file. The value of this field is used in alignment records.'], \
     ['SQ', 'LN', True,  None,      'Sequence length'], \
     ['SQ', 'AS', False, None,      'Genome assembly identifier. Refers to the reference genome assembly in an unambiguous form. Example: HG18.'], \
     ['SQ', 'M5', False, None,      'MD5 checksum of the sequence in the uppercase (gaps and space are removed)'], \
     ['SQ', 'UR', False, None,      'URI of the sequence'],  \
     ['SQ', 'SP', False, None,      'Species'], \
     ['RG', 'ID', True, 'ReadGroup1', 'Unique read group identifier. The value of the ID field is used in the RG tags of alignment records.'], \
     ['RG', 'SM', True, 'SampleA',  'Sample (use pool name where a pool is being sequenced)'], \
     ['RG', 'LB', False, None,      'Library'], \
     ['RG', 'DS', False, None,      'Description'], \
     ['RG', 'PU', False, None,      'Platform unit (e.g. lane for Illumina or slide for SOLiD)'], \
     ['RG', 'PI', False, None,      'Predicted median insert size (maybe different from the actual median insert size)'], \
     ['RG', 'CN', False, None,      'Name of sequencing center producing the read'], \
     ['RG', 'DT', False, None,      'Date the run was produced (ISO 8601 date or date/time)'], \
     ['RG', 'PL', False, None,      'Platform/technology used to produce the read'], \
     ['PG', 'ID', True,  'ProgramA','Program name'], \
     ['PG', 'VN', False, None,      'Program version'], \
     ['PG', 'CL', False, None,      'Command line']] 

def SAMHeaderEncode( type, tag ):
    return type + '.' + tag

SAMHeaderDict = dict()
for record in SAMHeaderFields:
    type, tag, req, default, desc = record
    SAMHeaderDict[SAMHeaderEncode(type, tag)] = [req, default, desc]

# -----------------------------------------------------------------------------------------------
#
# A SAM header is a potentially complex object, so we just punt on it.  The only really required
# outputs that *must* be user specified are the SQ: SN and SQ: LN fields.  
# Everything else is just split up and stored in our hash table by reading the SAMHeaderFields
# data structure (defined above).  All required fields are initialized to their default values
#
# -----------------------------------------------------------------------------------------------
class SAMHeader(dict):
    def __init__(self, seqName, seqLen, **keys):
        # setup initial header
        self.fields = dict()
        for record in SAMHeaderFields:
            type, tag, req, default, desc = record
            self.setField( type, tag, default )
        
        self.setField('SQ', 'SN', seqName )
        self.setField('SQ', 'LN', seqLen )
                        
        # add keyword arguments
        for key, value in keys.iteritems():
            type, tag = key.split('.')
            self.setField( type, tag, value )
        
    def isReq( self, type, tag ):
        return SAMHeaderDict[SAMHeaderEncode(type,tag)][0]
            
    def setField( self, type, tag, value ):
        #print 'Setting', type, tag, value
        self.fields[SAMHeaderEncode(type, tag)] = value
        
    def getField( self, type, tag ):
        return self.fields[SAMHeaderEncode(type, tag)]
        
    def __str__( self ):
        types = ['HD', 'SQ', 'RG']
        
        def formatType( type ):
            s = '@' + type + '\t'
            for record in SAMHeaderFields:
                qtype, tag, req, default, desc = record
                if type == qtype and self.isReq(type, tag):
                    value = self.getField(type, tag)
                    if value == None:
                        raise Error('Unset required SAM header ' + type + ' ' + tag)
                    s += tag + ':' + str(value) + '\t'
            return s
        return string.join(map( formatType, types ),'\n')

SAM_SEQPAIRED    = 0x0001 # the read is paired in sequencing, no matter whether it is mapped in a pair
SAM_MAPPAIRED    = 0x0002 # the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1 
SAM_UNMAPPED     = 0x0004 # the query sequence itself is unmapped 
SAM_MATEUNMAPPED = 0x0008 # the mate is unmapped 1 
SAM_QUERYSTRAND  = 0x0010 # strand of the query (0 for forward; 1 for reverse strand) 
SAM_MATESTRAND   = 0x0020 # strand of the mate 1 
SAM_ISFIRSTREAD  = 0x0040 # the read is the first read in a pair 1,2 
SAM_ISSECONDREAD = 0x0080 # the read is the second read in a pair 1,2 
SAM_NOTPRIMARY   = 0x0100 # the alignment is not primary (a read having split hits may have multiple primary alignment records) 

SAM_FLAGS = {
    SAM_SEQPAIRED       : 'the read is paired in sequencing, no matter whether it is mapped in a pair',
    SAM_MAPPAIRED       : 'the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1', 
    SAM_UNMAPPED        : 'the query sequence itself is unmapped',
    SAM_MATEUNMAPPED    : 'the mate is unmapped 1', 
    SAM_QUERYSTRAND     : 'strand of the query (0 for forward; 1 for reverse strand)' ,
    SAM_MATESTRAND      : 'strand of the mate 1' ,
    SAM_ISFIRSTREAD     : 'the read is the first read in a pair 1,2' ,
    SAM_ISSECONDREAD    : 'the read is the second read in a pair 1,2',
    SAM_NOTPRIMARY      : 'the alignment is not primary (a read having split hits may have multiple primary alignment records)'
    }
   
def SAMRecordFromArgs( qname, flags, rname, pos, mapq, cigar, seq, quals, pairContig = '*', pairPos = 0, insertSize = 0 ):
    r = SAMRecord()
    r.setValuesFromArgs( qname, flags, rname, pos, mapq, cigar, seq, quals, pairContig, pairPos, insertSize )
    return r

def SAMRecordFromString( str ):
    r = SAMRecord()
    r.setValuesFromString( str )
    return r
    
def SAMFlagValue( flags, testFlag ):
    return testFlag & flags

def SAMFlagIsSet( flags, testFlag ):
    return SAMFlagValue(flags, testFlag) <> 0

def SAMFlagsDescs( flags ):
    def keepMe(p):
        flagKey, flagDesc = p
        return [flagKey, SAMFlagIsSet(flags, flagKey), flagDesc]
    return sorted(map( keepMe, SAM_FLAGS.iteritems() ))            

# -----------------------------------------------------------------------------------------------
#
# This is really the meat of the SAM I/O system.  
#
# setValuesFromArgs takes the required arguments for a SAM record and stores them in a list
# (in SAM ordering) so that writing the values to a string is trivial.  Note the conversion
# of int quals to ASCII encoded quals.  
#
# -----------------------------------------------------------------------------------------------
class SAMRecord:
    def __init__( self, vals = []):
        self.vals = []
    
    def setValuesFromArgs( self, qname, flags, rname, pos, mapq, cigar, seq, quals, pairContig = '*', pairPos = 0, insertSize = 0 ):
        self.vals = [qname, self.formatFlags(flags), rname, pos, mapq, \
                     cigar, pairContig, pairPos, insertSize, seq.lower(), self.toASCII33(quals) ]

    def setValuesFromString( self, line ):
        #print 'line is', line
        formats = [str, int, str, int, int, str, str, int, int, str, str]
        self.vals = map( lambda f, v: f(v), formats, line.split()[0:len(formats)] )

    def formatFlags( self, flags ):
        b = reduce( operator.__or__, flags, 0 )
        #print 'FormatFlags', flags, b
        return b
        
    def getQName(self): return self.vals[0]
    def getFlags(self): return self.vals[1]
    def getRname(self): return self.vals[2] # Reference sequence name (can be chr1, chr2, etc...)
    def getPos(self): return self.vals[3]
    def getMapq(self): return self.vals[4]
    def getCigar(self): return self.vals[5] # Returns CIGAR which gives match, insertion, and other info as "75M2I"
    def getSeq(self): return self.vals[9]
    def getQuals(self): return self.fromASCII33(self.vals[10])

    def toASCII33( self, quals ):
        return string.join( map( lambda q: chr(q+33), quals ), '' )

    def fromASCII33( self, qualStr ):
        return map( lambda q: ord(q)-33, qualStr) 

    def __str__( self ):
        return string.join( map( str, self.vals ), '\t')

# -----------------------------------------------------------------------------------------------
#
# Wrapper class for reading and writing records to a SAM file
#
# -----------------------------------------------------------------------------------------------
class SAMIO:
    def __init__( self, fileName, header = None, debugging = False, func = None ):
        self.header = header
        self.debugging = debugging
        
        if func == None:
            func = lambda x: x
        self.func = func
        
        if self.header == None:
            self.fileHandle = pushback_file(fileName, "r")
            self.readHeader()
        else:
            self.fileHandle = file(fileName, "w")
            self.writeHeader()
            
    def close(self):
        self.fileHandle.close()
        
    def readHeader(self):
        self.header = ""
        while True:
            line = self.fileHandle.readline()
            if line.startswith("@"):
                self.header += line
            else:
                self.fileHandle.pushback(line)
                if self.debugging:
                    pass #print str(self.header)
                return True

        return self.header

    def RecordGenerator( self ):
        for line in self.fileHandle:
            #print 'Reading line', line
            line = line.rstrip("\n")
            yield self.func(SAMRecordFromString(line))
        raise StopIteration
        return
        
    def __iter__(self):
        return self.RecordGenerator()
        
    def writeRecord( self, recordIterator ):
        for record in recordIterator:
            self.writeRecord(record)
            
    def writeRecord( self, record ):
        if self.debugging:
            print record
        print >> self.fileHandle, record
        return True
        
    def writeHeader( self ):
        if self.debugging:
            print str(self.header)
        print >> self.fileHandle, str(self.header)
        return True
