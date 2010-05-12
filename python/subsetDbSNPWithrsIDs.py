import sys
from optparse import OptionParser

def getRsIDSet(dbSNPIds, idIndex):
    s = set()
    for line in open(dbSNPIds):
        #print line
        rsID = line.split()[idIndex]
        s.add(rsID)

    return frozenset(s)


def main():
    global OPTIONS
    usage = "usage: %prog [options] dbSNP.in.rsids dbSNP.to.match dbSNP.out"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="")
                        
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")

    dbSNPIds = args[0]
    dbSNPMatch = args[1]
    dbSNPOut = args[2]
    idIndex = 4

    rsSet = getRsIDSet(dbSNPIds, idIndex)
    print 'rsID set has %d elements' % len(rsSet)
    
    # 
    count = 0
    matched = 0 
    out = open(dbSNPOut, 'w')
    for line in open(dbSNPMatch):
        count += 1
        rsID = line.split()[idIndex]
        if rsID in rsSet:
            #sys.stdout.write(line)
            matched += 1
            out.write(line)
    print 'Processed %d lines, matching %d elements, excluding %d' % ( count, matched, count - matched )
    out.close()
            
if __name__ == "__main__":
    main()
