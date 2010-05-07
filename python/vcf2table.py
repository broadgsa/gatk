import os.path
import sys
from optparse import OptionParser
from vcfReader import *

if __name__ == "__main__":
    usage = "usage: %prog [file.vcf | if absent stdin] [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--f", dest="fields",
                        type='string', default=None,
                        help="Comma separated list of fields to exact")
    parser.add_option("-e", "--filter", dest="filter",
                        action='store_true', default=False,
                        help="If true, only includes records that aren't filtered in the output")
    parser.add_option("-s", "--s", dest="skip",
                        type='int', default=0,
                        help="Only print out every 1 / skip records")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) > 1:
        parser.error("incorrect number of arguments")

    counter = OPTIONS.skip
    src = sys.stdin
    if len(args) > 0:
        src = open(args[0])

    if OPTIONS.fields == None:
        sys.exit("Fields argument must be provided")

    fields = OPTIONS.fields.split(',')
    for header, vcf, count in lines2VCF(src, extendedOutput = True):
        #print vcf, count
        if count == 1 and vcf.hasHeader():
            print '\t'.join(fields)

        if counter > 0:
            counter -= 1                      
        else:
            counter = OPTIONS.skip
            if OPTIONS.filter and vcf.failsFilters():
                pass
            else:
                print '\t'.join([ str(vcf.getField(field, '0')) for field in fields])
