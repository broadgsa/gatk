from optparse import OptionParser
from math import log10
from math import floor

def foreach( appl, col ):
    for z in col:
        appl(z)

def parseWiggle(line):
    return int(line.strip())

def parseTDF(line):
    return int(line.split("Q: ")[1].strip())

def recalibrateJointly(files,JECDF):
    print("Not implemented.")

def recalibrate(fname,hist,out):
    print(fname)
    output = open(out + "." + fname.rsplit(".",1)[0].rsplit(".",1)[1], 'w')
    denom = sum(hist.values())
    cumQuals = dict()
    for key1 in hist:
        sumMore = 0
        for key2 in hist:
            if ( key2 >= key1 ):
                sumMore += hist[key2]
        cumQuals[key1] = min(int(-10*log10((0.0+sumMore)/denom)),150)
    use = parseTDF
    if ( fname.endswith(".wig") ):
        use = parseWiggle
    inFile = open(fname)

    if ( fname.endswith(".wig") ):
        output.write(inFile.readline())

    def rewrite(line,val):
        if ( line.find("Q: ") == -1 ):
            return str(val)+"\n"
        else:
            spline = line.split("Q: ")
            bef = spline[0]
            af = spline[1]
            afSp = af.split("\t")
            afSp[0] = str(val)
            return bef + "\t".join(afSp)

    foreach( lambda u: output.write(rewrite(u,cumQuals[use(u)])), inFile.readlines() )

def run(opt,arg):
    print(opt)
    print(arg)
    files = map(lambda u: open(u), opt.in_files)
    fn = dict(zip(files,opt.in_files))
    foreach( lambda u: u[1].readline(), filter( lambda u: u[0].endswith('.wig'), zip(opt.in_files,files)))
    isWig = map(lambda u: u.endswith('.wig'), opt.in_files)
    isWig = dict(zip(files,isWig))
    marginalCumulativeDists = dict(map( lambda u: [u,dict()], opt.in_files))

    def addVal(f,v):
        d = marginalCumulativeDists[fn[f]]
        if ( v not in d ):
            d[v] = 0
        d[v] += 1

    def parseLine(fl,line):
        if ( isWig[fl] ):
            return parseWiggle(line)
        else:
            return parseTDF(line)

    foreach( lambda f: foreach( lambda v: addVal(f,v), map(lambda u: parseLine(f,u), f.readlines() )), files )
    foreach( lambda f: f.close(), files )
    #print(marginalCumulativeDists)
    if ( opt.joint ):
        recalibrateJointly(opt.in_files,marginalCumulativeDists)
    else:
        foreach( lambda s: recalibrate(s,marginalCumulativeDists[s],opt.out), opt.in_files )

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("-I","--in",dest="in_files",help="file(s) to recalibrate",action="append")
    parser.add_option("-J","--recalibrateJointly",dest="joint",action="store_true",help="Recalibrate quality scores jointly across files, rather than for each independently. Assumes lines match up exactly.")
    parser.add_option("-O","--output",dest="out",action="store",help="The base name for the recalibrated output files")
    (options,args) = parser.parse_args()
    run(options,args)

if __name__ == "__main__":
    main()
