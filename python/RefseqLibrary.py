import math
def chrFormat(chr):
    if ( chr == "chr0"):
        return "chrM"
    if ( chr == "chr23" ):
        return "chrX"
    if ( chr == "chr24" ):
        return "chrY"
    else:
        return chr

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

    def bedFormat(self):
        return self.chromosome+"\t"+str(self.start)+"\t"+str(self.stop)+"\t+\ttarget_x"

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
        return "\t".join([chrFormat(self.interval.chromosome),str(self.interval.start),str(self.interval.stop),self.gene,self.id,str(self.getCoverageProportion())])

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

    def getExonIntervals(self):
        intervals = list()
        for exon in self.exons:
            intervals.append(exon.getInterval())
        return intervals

    def getBaseCoverage(self):
        coverage = 0
        for exon in self.exons:
            coverage = coverage + exon.getBaseCoverage()
        return coverage

    def __str__(self):
        exonString = list()
        for exon in self.exons:
            exonString.append(str(exon))
        return self.name+"\t"+"\t".join(exonString)

    def getGeneName(self):
        return self.name

    def setGeneName(self,newName):
        self.name = newName

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

def getRefseqGenes(names):
    names = list(names)
    refGene = open("/humgen/gsa-hpprojects/GATK/data/refGene.sorted.txt")
    refSeq = open("/humgen/gsa-hpprojects/GATK/data/refseq/hg18.ref_gene.cds.bed")
    refSeqGeneNames = list()
    refNamesToAltNames = dict()
    for name in names:
        if ( name.startswith("NM_") ):
            refSeqGeneNames.append(name)
            refNamesToAltNames[name]=name
    
    if ( len(names) > 0 ):
        for line in refGene.readlines():
            spline = line.strip().split("\t")
            altName = spline[len(spline)-4]
            if ( altName in names ):
                if ( not ( altName in refNamesToAltNames.values() ) ): 
                    refSeqGeneNames.append(spline[1])
                    refNamesToAltNames[spline[1]]=altName
                else:
                    print("WARNING: multiple transcripts found for gene "+altName+" using first available transcript from refseq export")

    if ( len(names) > len(refSeqGeneNames) ):
        for g in refSeqGeneNames:
            if ( refNamesToAltNames[g] in names ):
                names.remove(refNamesToAltNames[g])

        raise ValueError("No entry found for genes: "+str(names))

    # build up the gene list
    genes = dict()
    for geneName in refSeqGeneNames:
        genes[geneName] = Gene(geneName)

    for line in refSeq.readlines():
        spline = line.strip().split("\t")
        geneName = spline[3].split("_cds")[0]
        if ( geneName in refSeqGeneNames ):
            chrom = spline[0]
            start = int(spline[1])
            stop = int(spline[2])
            id = "cds_"+spline[3].split("_cds_")[1].split("_")[0]
            genes[geneName].addExon(Exon(geneName,id,chrom,start,stop))

    toReturn = list()
    for gene in genes.values():
        gene.setGeneName(refNamesToAltNames[gene.getGeneName()])
        toReturn.append(gene)
    return toReturn


def getIntervalHeaderLines():
    whole_exome_file = open("/humgen/gsa-hpprojects/GATK/data/whole_exome_agilent_1.1_refseq_plus_3_boosters.targets.hg18.interval_list")
    header = list()
    line = whole_exome_file.readline()
    while ( line.startswith("@") ):
        header.append(line)
        line = whole_exome_file.readline()
    whole_exome_file.close()
    return header

def parseDesignFile(file):
    designFile = open(file)
    genes = dict()
    for line in designFile.readlines():
        if ( line.startswith("TARGET") ):
            spline = line.strip().split()
            if ( spline[1].startswith("chr") ):
                chrom = spline[1]
            else:
                chrom = "chr"+spline[1]
            start = 1 + int(spline[2])
            stop = 1 + int(spline[3])
            gene_name = spline[4].split("#")[1].split("_")[0]
            try:
                exon_id = spline[4].split(gene_name)[1]
            except IndexError:
                exon_id = "_".join(spline[5:len(spline)-1])
            exon = Exon(gene_name,exon_id,chrom,start,stop)
            if ( gene_name in genes.keys() ):
                genes[gene_name].addExon(exon)
            else:
                genes[gene_name] = Gene(gene_name)
                genes[gene_name].addExon(exon)
    return genes.values()
