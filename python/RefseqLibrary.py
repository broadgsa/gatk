import math
def chr2int(chr):
    try:
        chr = chr.split("chr")[1]
        if ( chr == "M" ):
            return 0
        if ( chr == "X" ):
            return 23
        if ( chr == "Y" ):
            return 24
        return int(chr)
    except IndexError:
        print("Index error: "+chr)
        return -1

def intervalCompare(int1,int2):
    chr1 = chr2int(int1.chromosome)
    chr2 = chr2int(int2.chromosome)
    if ( chr1 < chr2 ):
        return -1
    elif ( chr1 > chr2 ):
        return 1
    else:
        start1 = int1.start
        start2 = int2.start
        return start1 - start2

class Interval:
    def __init__(self,chrom,start,stop):
        self.chromosome = chrom
        self.start = start
        self.stop = stop
        if ( chrom == "NONE" ):
            self.isEmpty = True
        else:
            self.isEmpty = False
        self.basesCovered = None

    def size(self):
        return self.stop - self.start

    def overlaps(self,other):
        if ( self.chromosome == other.chromosome ):
            if ( other.stop < self.stop and not other.start < self.start ):
                return True
            if ( other.stop > self.start and not other.start > self.stop ):
                return True
        return False

    def intersect(self,other):
        if ( self.overlaps(other) ):
            return Interval(self.chromosome, max(self.start,other.start), min(self.stop,other.stop))
        else:
            return Interval("NONE",-1,-1)

    def isBefore(self,other):
        if ( chr2int(self.chromosome) < chr2int(other.chromosome) ):
            return True
        elif ( chr2int(self.chromosome) > chr2int(other.chromosome) ):
            return False
        else:
            if ( other.start > self.stop ):
                return True
            return False

    def __str__(self):
        return self.chromosome + ":" + str(self.start) + "-" + str(self.stop)

    def __cmp__(self,other):
        return intervalCompare(self,other)

class CoveredInterval(Interval):
    def __init__(self,chrom,start,stop):
        Interval.__init__(self,chrom,start,stop)
        self.overlappingSubIntervals = list()

    def updateCoverage(self,other):
        if ( other.overlaps(self) ):
            self.overlappingSubIntervals.append(other)

    def getBaseCoverage(self):
        if ( self.basesCovered is None ):
            basesCovered = 0
            intersects = list()
            self.overlappingSubIntervals.sort(intervalCompare)
            for overlap in self.overlappingSubIntervals:
                ival = self.intersect(overlap)
                intersects.append(ival)
            for i in range(len(intersects)):
                basesCovered = basesCovered + intersects[i].size()
                for j in range(i+1,len(intersects)):
                    basesCovered = basesCovered - (intersects[i].intersect(intersects[j])).size()

            self.basesCovered = basesCovered
        return self.basesCovered

    def getOverlappingIntervals(self):
        return self.overlappingSubIntervals

class Exon:
    def __init__(self,geneName,exonid,chrom,start,stop):
        #print("Adding exon for "+geneName+" with "+chrom+":"+str(start)+"-"+str(stop))
        self.interval = CoveredInterval(chrom,start,stop)
        self.gene = geneName
        self.id = exonid

    def getOverlappingIntervals(self):
        return self.interval.getOverlappingIntervals()

    def updateCoverage(self,target):
        self.interval.updateCoverage(target)

    def overlaps(self,target):
        return self.interval.overlaps(target)

    def isBefore(self,target):
        return self.interval.isBefore(target)

    def getBaseCoverage(self):
        return self.interval.getBaseCoverage()

    def size(self):
        return self.interval.size()

    def getBedEntry(self):
        return "\t".join([self.interval.chromosome,str(self.interval.start),str(self.interval.stop),self.gene,self.id,str(self.getCoverageProportion())])

    def getCoverageProportion(self):
        if ( self.size() > 0 ):
            return float(self.getBaseCoverage())/float(self.size())
        else:
            return -1

    def __str__(self):
        return self.gene+"("+self.id+") "+str(self.interval)

    def getInterval(self):
        return self.interval

class Gene:
    def __init__(self,name):
        self.name = name
        self.exons = list()

    def addExon(self,exon):
        self.exons.append(exon)

    def size(self):
        size = 0
        for exon in self.exons:
            size = size + exon.size()
        return size

    def getBaseCoverage(self):
        coverage = 0
        for exon in self.exons:
            coverage = coverage + exon.getBaseCoverage()
        return coverage

    def __str__(self):
        exonString = list()
        for exon in exons:
            exonString.append(str(exon))
        return name+"\t"+"\t".join(exonString)

class ExonRecord(Exon):
    def __init__(self,geneName,exonid,chrom,start,stop,prop):
        Exon.__init__(self,geneName,exonid,chrom,start,stop)
        self.coverageProportion = prop
        self.baseCoverage = math.ceil(prop*self.size())
        self.records = list()

    def getBaseCoverage(self):
        return self.baseCoverage

    def getCoverageProportion(self):
        return self.coverageProportion

    def addRecord(self,record):
        self.records.append(record)

    def getData(self):
        toRet = ""
        for record in self.records:
            toRet += record.dataString()+"\t"
        return toRet
class CoverageRecord:
    def __init__(self,chrom,start,stop,sampleName,mean,median,q1,q3):
        self.interval = Interval(chrom,start,stop)  #indexing issues
        self.sampleName = sampleName
        self.mean = mean
        self.median = median
        self.q1 = q1
        self.q3 = q3

    def dataString(self):
        return "\t".join([self.sampleName,str(self.mean),str(self.median),str(self.q1),str(self.q3)])

    def getInterval(self):
        return self.interval
