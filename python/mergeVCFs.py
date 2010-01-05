import os.path
import sys
from optparse import OptionParser
from vcfReader import *
from itertools import *
import faiReader

def main():
    global OPTIONS
    usage = "usage: %prog [options] file1 ... fileN"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--f", dest="fai",
                        type='string', default=None,
                        help="FAI file defining the sort order of the VCF")
    parser.add_option("-a", "--assumeSorted", dest="assumeSorted",
                        action='store_true', default=False,
                        help="If provided, this assumes the input VCF files are themselves sorted, enabling a simple efficent merge")         
         
    (OPTIONS, args) = parser.parse_args()
    if len(args) == 0:
        parser.error("Requires at least 1 VCF to merge")

    order = None
    if OPTIONS.fai <> None: order = faiReader.readFAIContigOrdering(OPTIONS.fai)
    #print 'Order', order

    if OPTIONS.assumeSorted:
        mergeSort(args, order)
    else:
        memSort(args, order)

def cmpVCFRecords(order, r1, r2):
    if order <> None:
        c1 = order[str(r1.getChrom())]
        c2 = order[str(r2.getChrom())]
        orderCmp = cmp(c1, c2)
        if orderCmp <> 0:
            return orderCmp
    return cmp(r1.getPos(), r2.getPos())

def mergeSort(args, order):
    #print 'MergeSort', args, order
    header = None
    
    orderMap = []
    for file in args:
        #print file
        openedFile = open(file)
        for header, record, counter in lines2VCF(openedFile, extendedOutput = True, decodeAll = False):
            orderMap.append([record, file])
            break
        openedFile.close()

    #print orderMap
    sortedOrderMap = sorted(orderMap, key=lambda x: x[0], cmp = lambda r1, r2: cmpVCFRecords(order, r1, r2))
    #print sortedOrderMap

    for headerLine in header: print headerLine
    for file in map( lambda x: x[1], sortedOrderMap):
        #print file
        for record in lines2VCF(open(file), extendedOutput = False, decodeAll = False):
            print record.format()
        
def memSort(args, order):
    header = None
    records = []
    for file in args:
        #print file
        for header, record, counter in lines2VCF(open(file), extendedOutput = True, decodeAll = False):
            records.append(record)

    records.sort(lambda r1, r2: cmpVCFRecords(order, r1, r2)) 
    for line in formatVCF(header, records):
        #pass
        print line

PROFILE = False
if __name__ == "__main__":
    if PROFILE:
        import cProfile
        cProfile.run('main()', 'fooprof')
        import pstats
        p = pstats.Stats('fooprof')
        p.sort_stats('cumulative').print_stats(10)
        p.sort_stats('time').print_stats(10)
        p.sort_stats('time', 'cum').print_stats(.5, 'init')
    else:
        main()