import os
base1 = "/humgen/gsa-hphome1/chartl/projects/exome/coverage/whole_exome_broad/broad/"
base2 = "/humgen/gsa-hphome1/chartl/projects/exome/coverage/whole_exome_broad/intersect/"

def getField(fname,header):
    f = open(fname)
    h = f.readline().split("\t")
    idx = h.index(header)
    return map(lambda x: x.strip().split("\t")[idx],filter(lambda u: not u.startswith("Tot"),f.readlines()))

def getFiles(base):
    num_jobs = len(os.listdir(base+"GATKDispatcher/"))
    sum_file = base+"GATKDispatcher/dispatch%d/job%d.sample_summary"
    int_file = base+"GATKDispatcher/dispatch%d/job%d.sample_interval_summary"
    files = map( lambda x: (sum_file % (x,x), int_file % (x,x)), range(num_jobs))
    return files

def getBases(base):
    files = getFiles(base)
    lines = map(lambda x: getField(x[1],"Target"),files)
    sizes = map( lambda iList: reduce ( lambda u,v: u + v , map(lambda ival: int(ival.split(":")[1].split("-")[1])-int(ival.split(":")[1].split("-")[0]), iList ) ) , lines )
    return(sizes)

def getPctAbove(base):
    files = getFiles(base)
    pct_above = map(lambda y: map(lambda z: float(z),y),map(lambda x: getField(x[0],"%_bases_above_20"),files))
    return pct_above

def getTotals(base):
    files = getFiles(base)
    totals = map(lambda y: map(lambda z: int(z),y),map(lambda x: getField(x[0],"total"),files))
    return totals

bases = getBases(base1) + getBases(base2)
totals = getTotals(base1) + getTotals(base2)
pct_above = getPctAbove(base1) + getPctAbove(base2)
files = getFiles(base1) + getFiles(base2)

total_bases = reduce(lambda u,v: u + v , bases)

total_cvg = map(lambda u: 0, range(len(pct_above[0])))
total_covered_20 = map(lambda u: 0, range(len(pct_above[0])))
for f_idx in range(len(files)):
    for sam_idx in range(len(pct_above[0])):
        total_covered_20[sam_idx] += pct_above[f_idx][sam_idx]*bases[f_idx]
        total_cvg[sam_idx] += totals[f_idx][sam_idx]

pct_above_20 = map( lambda x: x/total_bases, total_covered_20)
mean_cvgs = map( lambda x: (x+0.0)/total_bases, total_cvg)
#print(total_bases)
#print(total_cvg)
print("Mean\t%_above_20")
mean_str = map(lambda u: "%.1f" % u, mean_cvgs)
above_str = map(lambda u: "%.1f" % u, pct_above_20)
print(reduce( lambda u,v: u + "\n" + v, map(lambda idx: mean_str[idx]+"\t"+above_str[idx],range(len(mean_str)))))
