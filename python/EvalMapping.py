#!/usr/bin/env python

# Evaluates mapping results produced my RunMapper.py by comparig
# the result to the original simulated reads in a sam file

import getopt, sys, os, re
from operator import attrgetter, itemgetter
from subprocess import Popen, STDOUT, PIPE
from SAM import SAMIO
from pushback_file import pushback_file
from qltout import qltout_file
from aln_file import aln_file


#!/usr/bin/env python

import string, sys

class SAMRecordWrapper:
    def __init__(self, rec):
        #print 'Rec', rec
        self.rec = rec
    
    def id(self): return self.rec.getQName()
    def contig(self): return self.rec.getRname()
    def pos(self): return self.rec.getPos()
    def offset(self): return self.pos()-1
    def map_qual(self): return self.rec.getMapq()
    
    def __str__(self): return str(self.rec)

def SAMIOWrapper( file ):
    print 'SAMIOWrapper', file
    return SAMIO( file, debugging=True, func=SAMRecordWrapper)

class mapper_info:
    def __init__(self, name, ext, file_class):
        self.name = name.lower()
        self.ext = ext
        self.file_class = file_class

mapper_defs = [ mapper_info("ilt", "ilt.qltout", qltout_file),
            mapper_info("merlin", "qltout", qltout_file),
            mapper_info("swmerlin", "qltout", qltout_file),
            mapper_info("bwa", "sam", SAMIOWrapper),
            mapper_info("bwa32", "sam", SAMIOWrapper),
            mapper_info("maq", "aln.txt", aln_file) ]
mapper_defs = dict( [(m.name, m) for m in mapper_defs] )

def Comp_sam_map(comp_head, mapper, file_output, sam_indel):
    
    # sam_indel is a +/- number corresponding to the size of insertions / deletions
    # in generated SAM file that are fixed (1,2,4,8,16...), this doesn't always
    # translate to the SAM's cigar, but the full value should be used

    file_output_head = comp_head+"."+mapper.name
    if file_output:
        # if non-false, output to filename in file_output_head
        # if "" or False, leave output going to stdout
        #saveout = sys.stdout
        feval = open(file_output_head+".eval", "w")
        #sys.stdout = feval

        # create a file with a one-line tabular output of stats 
        # (for merging into Excel / R input)
        fevaltab = open(file_output_head+".tab", "w")

    # if maq file, convert aln.map to aln.txt
    aln_txt_filename = file_output_head+"."+mapper.ext
    if mapper.name == "maq" and not os.path.exists(aln_txt_filename):
        maq_exe = "/seq/dirseq/maq-0.7.1/maq"
        
        # SHOULD BE THIS
        #cmd_str = maq_exe+" mapview "+file_output_head+".out.aln.map"

        cmd_str = maq_exe+" mapview "+file_output_head+".out.aln.map"

        print >> sys.stderr, "Executing "+cmd_str
        fout = open(aln_txt_filename, "w")
        p = Popen( cmd_str , shell=True, stdout=fout)

    map_filename = comp_head+"."+mapper.name+"."+mapper.ext
    print 'map_filename', map_filename
    map_file = mapper.file_class(map_filename)

    SAM_filename = comp_head+".sam"
    sam_file = SAMIO(SAM_filename, debugging=True)

    #sams = [s for s in sam_file]
    print "Reading map file..."
    maps = [q for q in map_file]

    print "Reading SAM file..."
    sams = {}
    if mapper.name in ["maq", 'bwa', 'bwa32']:
        for s in sam_file:
            sams[s.getQName()] = s
    else:
        for i, s in enumerate(sam_file):
            sams[i] = s
            
    #maps = {}
    #for m in map_file:
    #    maps[m.id()] = m
    #sams = [s for s in sam_file]

    class smq_t:
        def __init__(self, sam, map, qual):
            self.sam = sam # SAM record (exists for all reads)
            self.map = map # Mapping record, may be None
            self.qual = qual # Quality of mapping

    print "Making SMQ records..."
    SMQ = []
    for maprec in maps:
        samrec = sams.get(maprec.id()) # get mapping otherwise None
        SMQ.append( smq_t(samrec, maprec, maprec.map_qual()) )

    print "Sorting..."
    SMQ.sort(key=attrgetter('qual'), reverse=True)

    print "Evaluating placements..."
    placed = 0
    incorrectly_placed = 0
    correctly_placed = 0
    fevaltab.write( "\t".join( ["last_qual", "total_reads", "placed", "correctly_placed", "incorrectly_placed", "not_placed"] ) + "\n" )
    if len(SMQ):
        last_qual = SMQ[0].qual # grab 1st qual
        for smq in SMQ:
            if smq.sam == None:
                continue
        
            if smq.qual != last_qual:
                total_reads = len(sams)#placed + not_placed
                not_placed = len(sams) - placed
                fevaltab.write( "\t".join( map(str, [last_qual, total_reads, placed, correctly_placed, incorrectly_placed, not_placed]) ) + "\n" )
                #print "Wrote all reads with qual >= "+str(last_qual)
                last_qual = smq.qual

            placed += 1

            # Convert the CIGAR string - e.g. 75M1I, 74M2D - into an int corresponding to bases
            # inserted / deleted - e.g. +1, -2
            #cigar = smq.sam.getCigar()
            #indelled_bases = 0
            #match = re.search(r'(\d+)I', cigar)
            #if match:
            #    indelled_bases += int(match.groups()[0])
            #match = re.search(r'(\d+)D', cigar)
            #if match:
            #    indelled_bases -= int(match.groups()[0])

            # We consider a read properly aligned if 
            # - both start positions are the same OR
            # - the end position of the mapping = end position of the samq
            #if smq.sam.getPos() != smq.map.pos() and smq.sam.getPos() - indelled_bases != smq.map.pos():
                
            try:
                map_indel_bases = smq.map.indelled_bases() # total number of indels in mapped seq
                # are generally at the beginning of an alignment spanning the anchor base
            except AttributeError:
                map_indel_bases = 0 # Since MAQ doesn't do local alignment, we don't care 
                # if its sam.py interface doesn't have access to indel data since there are 
                # no indels in its alignmnets

            if smq.sam.getPos() != smq.map.pos() and smq.map.pos() + sam_indel + map_indel_bases != smq.sam.getPos():
                #print >> feval, "Indelled SAM bases: "+str(indelled_bases)
                try:
                    print >> feval, "Indelled map bases: "+str(map_indel_bases)
                except AttributeError:
                    pass
                print >> feval, "SAM ref: %5s %s" % (smq.sam.getRname(), smq.sam.getPos())
                print >> feval, "Mapping: %5s %s" % (smq.map.contig(), smq.map.pos())
                print >> feval, "SAM record:\n"+str(smq.sam)
                print >> feval, "Mapping record:\n"+str(smq.map)+"\n"
                incorrectly_placed += 1
            else:
                correctly_placed +=1
            
            
    #indexed_maps = [None] * len(sams)
    #for q in maps:
    #    indexed_maps[q.id()] = q
            
    #def f(s,q): 
    #    if q == None: 
    #        return s
    #    else:
    #        return False

    #missing = filter( None, map(f, sams, indexed_maps))#, range(len(indexed_maps))) )
    

    #for miss in missing:
    #    print miss
    #not_placed = len(missing)
    #not_placed = len(sams) - len(maps)
    total_reads = len(sams)  #placed + not_placed
    not_placed = len(sams) - placed

    print >> feval, "Total reads  : "+str(total_reads)
    print >> feval, "+-Placed     : "+str(placed)
    print >> feval, "| +-Correct  : "+str(correctly_placed)
    print >> feval, "| +-Incorrect: "+str(incorrectly_placed)
    print >> feval, "+-Not placed : "+str(not_placed)

    fevaltab.write( "\t".join( map(str, [0, total_reads, placed, correctly_placed, incorrectly_placed, not_placed]) ) + "\n" )
    print "Wrote all reads with all quals"

    #if file_output_head:
    #    sys.stdout = saveout

def remove_escape_codes(s):
    import re
    return re.sub(r'\x1b\[(\d\d)?\w', "", s)

def usage():
    print "EvalMapping.py\n"
    print "Required arguments:\n"
    print "  -h HEAD    All input and output files start with HEAD"
    print "   --OR--"
    print "  <list of *.sam files on STDIN to use as filename templates>"
    print "    i.e. ls *.sam | EvalMapping.py -m MAPPER"
    print
    print "  -m MAPPER  Mapping program output to compare (e.g. merlin, ILT)"
    print
    print "Optional arguments:\n"
    print 
    print "  -o         Output to screen instead of HEAD.MAPPER.eval"
    print 
    print "Input files:"
    print "  HEAD.sam"
    print 
    print "Result files expected for Merlin and ILT:"
    print "  HEAD.MAPPER.qltout"
    print 
    print "Result files expected for MAQ:"
    print
    sys.exit(2)

if __name__ == "__main__":
    opts = None
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:m:o", ["head","mapper","file-output"])
    except getopt.GetoptError:
        usage()

    comp_heads = []
    mapper_str = ""
    file_output = True
    for opt, arg in opts:
        if opt in ("-h", "--head"):
            comp_heads = [arg]
        if opt in ("-m", "--mapper"):
            mapper_str = arg
            possible_mappers = mapper_defs.keys()
            # Check it later now...
            #if mapper_str not in possible_mappers:
            #    print "--mapper argument was \""+mapper_str+"\"; expected one of the following: "+", ".join(possible_mappers)
            #    sys.exit(2)
        if opt in ("-o", "--file-output"):
            file_output = False

    # Check for required args
    if (mapper_str == ""):
        usage()
    elif (mapper_str == "all"):
        mappers = mapper_defs.values()
        print mappers
    else:
        try:
            mapper_str_items = mapper_str.split(",")
            mappers = [mapper_defs.get(mapper.lower()) for mapper in mapper_str_items]
        except KeyError:
            print "Don't know this mapper given in -m option: "+mapper
            print "Try: ilt, merlin, maq, bwa"
            sys.exit(1)

    # if we didn't already get a single file from -h 
    if len(comp_heads) == 0:
        # and there's a list on STDIN,
        if not sys.stdin.isatty(): # redirected from file or pipe
            comp_heads = [ os.path.splitext( remove_escape_codes(file).rstrip() )[0] for file in sys.stdin.readlines() ]
            comp_heads = [ f for f in comp_heads if len(f) > 0 ]
        else:
            # no single file AND no STDIN file list
            usage()
    
    print "\nFile heads to evaluate:"
    print "\n".join(comp_heads)+"\n"

    comp_heads.reverse()
    for comp_head in comp_heads: #[:1]
        for mapper in mappers:
            # Guess the original size of simulated indel from filename
            m = re.search(r'_den\.(\d+)_', comp_head)
            if m:
                if ".DELETION_" in comp_head:
                    indel_size = -1 * int(m.groups()[0])
                elif ".INSERTION_" in comp_head:
                    indel_size =  int(m.groups()[0])
                else:
                    indel_size = 0
            else:
                sys.exit("Couldn't guess simulated read indel size from filename")
    
            print "Evaluating "+mapper.name+" results with indel size: "+str(indel_size)+" and file header: "+comp_head
            Comp_sam_map(comp_head, mapper, file_output, indel_size)

