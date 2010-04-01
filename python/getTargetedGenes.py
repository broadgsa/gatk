import RefseqLibrary as exLib

def chr2int(chr):
    chr = chr.split("chr")[1]
    if ( chr == "M" ):
        return 0
    if ( chr == "X" ):
        return 23
    if ( chr == "Y" ):
        return 24
    return int(chr)

def intervalCompareExon(exon1,exon2):
    int1 = exon1.interval
    int2 = exon2.interval
    return intervalCompare(int1,int2)

def intervalCompareExonTarget(exon1,target):
    int1 = exon1.interval
    return intervalCompare(int1,target)

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

test_first = exLib.Interval("chr2",100000,100001)
test_second  = exLib.Interval("chr2",100001,100004)
shouldBeNegative = intervalCompare(test_first,test_second)
if ( not shouldBeNegative < 0 ):
    raise ValueError("Interval compare is not performing properly.")

test_first = exLib.CoveredInterval("chr2",10400,10800)
test_first.updateCoverage(exLib.Interval("chr2",10387,10420))
test_first.updateCoverage(exLib.Interval("chr2",10410,10560))
test_first.updateCoverage(exLib.Interval("chr2",10600,10820))
test_first.updateCoverage(exLib.Interval("chr2",10890,10990))
if ( test_first.getBaseCoverage() != 360 ):
    print("Base coverage of test was "+str(test_first.getBaseCoverage())+" (400 expected)")
    print("Testing intersection....")
    g = exLib.Interval("chr2",10387,10420).intersect(exLib.Interval("chr2",10410,10560))
    if ( g.chromosome != "chr2" or g.start != 10410 or g.stop != 10420 ):
        print("Bad intersection! "+str(g))
    if ( g.size() != 10 ):
        print("Size not performing correctly")
    raise ValueError("Covered interval is not performing properly")

refseq_exons = open("/humgen/gsa-hpprojects/exome/gene_interval_lists/refseq_exons/hg18.ref_gene.cds.bed")
refseq_utr3 = open("/humgen/gsa-hpprojects/exome/gene_interval_lists/refseq_exons/hg18.ref_gene.utr3.bed")
refseq_utr5 = open("/humgen/gsa-hpprojects/exome/gene_interval_lists/refseq_exons/hg18.ref_gene.utr5.bed")
cancer_6000 = open("/seq/references/HybSelOligos/tcga_6k_genes.design")

counter = 0
exons = list()
def isNotValid(rline):
    try:
        n = chr2int(line.split()[0])
    except ValueError:
        return True
    return False

def parseLine(rline,token):
    spline = line.strip().split()
    chrom = spline[0]
    start = int(spline[1])
    stop = int(spline[2])
    longname = spline[3]
    geneName = longname.split("_"+token)[0]
    exonID = longname.split("_"+token+"_")[1].split("_chr")[0]+"_"+"_".join([token,longname[len(longname)-1]])
    return exLib.Exon(geneName,exonID,chrom,start,stop)

for line in refseq_exons.readlines():
    if ( isNotValid(line) ):
        continue
    counter = 1 + counter
    exons.append(parseLine(line,"cds"))
    if ( counter % 25000 == 0 ):
        print(str(counter)+" exons added")

for line in refseq_utr3:
    if ( isNotValid(line) ):
        continue
    counter = 1 + counter
    exons.append(parseLine(line,"utr3"))
    if ( counter % 25000 == 0 ):
        print(str(counter)+" exons added")

for line in refseq_utr5:
    if ( isNotValid(line) ):
        continue
    counter = 1 + counter
    exons.append(parseLine(line,"utr5"))
    if ( counter % 25000 == 0 ):
        print(str(counter)+" exons added")

for line in cancer_6000:
    if ( not line.startswith("TARGET") ):
        continue
    counter = 1 + counter
    spline = line.strip().split()
    chrom = "chr"+spline[1]
    start = int(spline[2])
    stop = int(spline[3])
    longname = spline[4]
    gene_name = longname.split("_")[0].split("#")[1]
    exon_id = "tcga_"+longname.split("_")[1]
    exons.append(exLib.Exon(gene_name,exon_id,chrom,start,stop))
    if ( counter % 25000 == 0 ):
        print(str(counter)+" exons added")

exons.sort(intervalCompareExon)
print(str(len(exons))+" exons added")
start_index = 0
counter = 0
import sys
interval_list = sys.argv[1]
headerlines = list()
for line in open(interval_list).readlines():
    if ( not line.startswith("@") and not line.startswith("#") ):
        counter = 1 + counter
        s = line.split()
        target = exLib.Interval(s[0],int(s[1]),int(s[2]))
        while ( exons[start_index].isBefore(target) and start_index < len(exons)-1):
            start_index = 1 + start_index
        index = start_index
        while(exons[index].overlaps(target) and index < len(exons) - 1):
            exons[index].updateCoverage(target)
            index = 1 + index
        if ( counter % 25000 == 0 ):
            print("Read "+str(counter)+" lines from interval list")
        #if ( counter % 5000 == 0 ):
        #    break
    else:
        if ( line.startswith("@") ):
            headerlines.append(line)

print("Done reading interval file. Creating genes.")
exons.sort(key = lambda x: x.gene)
genes = list()
counter = 0
prevName = "@T#h12_3"
nGenes = 0
for exon in exons:
    counter = 1 + counter
    if ( exon.gene != prevName ):
        genes.append(exLib.Gene(exon.gene))
        nGenes = 1 + nGenes
        genes[nGenes - 1].addExon(exon)
        prevName = exon.gene
    else:
        genes[nGenes - 1].addExon(exon)
    if ( counter % 20000 == 0 ):
        print("Processed "+str(counter)+" exons. "+str(len(genes))+" gene records created.")

def byProportion(gene1,gene2):
    gene1targets = gene1.size()
    gene1cvg = gene1.getBaseCoverage()
    gene1prop = float(gene1cvg)/float(gene1targets)
    gene2targets = gene2.size()
    gene2cvg = gene2.getBaseCoverage()
    gene2prop = float(gene2cvg)/float(gene2targets)
    if ( gene1prop > gene2prop ):
        return -1
    elif ( gene1prop < gene2prop ):
        return 1
    return 0

genes.sort(byProportion)

print("Writing coverage table...")
coverage = open("geneCoverage.txt",'w')
counter = 0
exonsInWellCoveredGenes = list()
for gene in genes:
    try:
        coverage.write(gene.name+"\t"+str(float(gene.getBaseCoverage())/float(gene.size()))+"\t"+str(gene.size())+"\t"+str(gene.getBaseCoverage()))
        if ( float(gene.getBaseCoverage())/float(gene.size()) > 0.8 ):
            for exon in gene.exons:
                exonsInWellCoveredGenes.append(exon)
        coverage.write("\n")
    except ZeroDivisionError:
        continue
    counter = 1 + counter
    if ( counter % 5000 == 0 ):
        print("Written: "+str(counter)+" genes.")
coverage.close()
print("Sorting exons by start location...")
exons.sort(intervalCompareExon)
exonsInWellCoveredGenes.sort(intervalCompareExon)
print("Writing all exon targets...")
exon_targets = open("exonTargets.interval_list",'w')
exon_targets.write("".join(headerlines))
for exon in exons:
        exon_targets.write(exon.getBedEntry()+"\n")
exon_targets.close()
print("Writing exons for well-covered genes...")
good_exon_targets = open("exonTargets_well_covered_genes.interval_list",'w')
good_exon_targets.write("".join(headerlines))
for exon in exonsInWellCoveredGenes:
    good_exon_targets.write(exon.getBedEntry()+"\n")
good_exon_targets.close()
print("Writing all exons and their overlapping targets...")
debug = open("overlapping_targets.txt",'w')
for exon in exons:
    debug.write(exon.id+"\t"+str(exon.getCoverageProportion())+"\t"+str(exon.interval)+"\t")
    intervals = list()
    for inter in exon.getOverlappingIntervals():
        intervals.append(str(inter))
    debug.write("\t".join(intervals)+"\n")
