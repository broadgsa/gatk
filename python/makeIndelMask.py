from optparse import OptionParser

def main():
    global OPTIONS
    usage = "usage: %prog [options] indelCalls.bed maskSize indelsMask.bed"
    parser = OptionParser(usage=usage)
    #parser.add_option("", "--dry", dest="dry",
    #                    action='store_true', default=False,
    #                    help="If provided, nothing actually gets run, just a dry run")
                       
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")

    indelCalls, maskSize, indelsMask = args
    maskSize = int(maskSize)

    out = open(indelsMask, 'w')
    for line in open(indelCalls):
        # chr1    71996   72005   -AAAAAAAAA:4/6
        chr, indelStart, indelStop, notes = line.split()
        maskStart = int(indelStart) - maskSize
        maskStop = int(indelStop) + maskSize
        maskNotes = notes + ":+/-" + str(maskSize)
        print >> out, '\t'.join([chr, str(maskStart), str(maskStop), maskNotes])
    out.close()

if __name__ == "__main__":
    main()
