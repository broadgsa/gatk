print("opening...")
ea_vcf = open("/humgen/gsa-hpprojects/analysis/privateMutations/eomi+autism/resources/callsets/eomi+autism/eomi+autism_batch.merged.vcf")
print("reading past header...")
ln = ea_vcf.readline()
import random
while (ln.startswith("#")):
    ln = ea_vcf.readline()

numcalc = 0
nread = -1

def getAC(e):
    if ( e.startswith("0/0") ):
        return 0
    elif ( e.startswith("0/1") or e.startswith("1/0") ):
        return 1
    elif ( e.startswith("1/1") ):
        return 2
    else:
        print("Warning: "+e)
        return 0

def calcTrans(line):
    spline = line.strip().split("\t")
    gtypes = filter(lambda y: y.find("./.") == -1, spline[9:len(spline)])
    if ( len(gtypes) < 1800 ):
        return (-1,-1)
    random.shuffle(gtypes)
    firstAC = reduce(lambda x,y: x + y , map(lambda u: getAC(u),gtypes[0:900]))
    if ( firstAC > 5 ):
        return (-1,-1)
    secondAC = reduce(lambda x,y: x + y, map(lambda u: getAC(u),gtypes[900:1800]))
    return (firstAC,secondAC)

print("Calculating...")
counts = filter(lambda u: u[0] > -1, map(lambda z: calcTrans(z) ,ea_vcf.readlines()))
print("Lines actually processed: %d" % len(counts))

cdict = dict()
for c in counts:
    if ( not c in cdict ):
        cdict[c] = 0
    cdict[c] += 1

out = open("posterior_counts.txt",'w')
for c in cdict:
    out.write("%d\t%d\t%d\n" % (c[0],c[1],cdict[c]))
