#!/usr/bin/env python

import os
import sys

input_verbose_beds = sys.argv[1].split(",")
output_vcf_name = sys.argv[2]

vcf_header = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
format = "GT:FT"
class Indel:
    def __init__(self,chrom,start,stop,bases,sample,isHet,isDeletion,isFiltered,isCoding,
                 consF,consR,refF,refR,consMM,refMM):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.bases = bases
        self.isDeletion = isDeletion
        self.isFiltered = isFiltered
        self.isHet = isHet
        self.sample = sample
        self.isCoding = isCoding
        self.consForward = int(consF)
        self.consReverse = int(consR)
        self.refForward = int(refF)
        self.refReverse = int(refR)
        self.consMM = float(consMM)
        self.refMM = float(refMM)

    def getAnnotations__DO_NOT_USE(self):
        ## note: log(1/2) ~~ -0.3
        forwardLod = -0.3*(self.consForward+self.refForward) + 2*self.consForward
        reverseLod = -0.3*(self.consReverse+self.refReverse) + 2*self.consReverse
        totalLod = -0.3*(self.consForward+self.consReverse+self.refForward+self.refReverse) + 2*(self.consForward+self.consReverse)
        strand_score = max(forwardLod-totalLod,reverseLod-totalLod)
        SB = "SB="+str(strand_score)
        mismatch = "MMR="+str(self.consMM)
        ref_mismatch = "refMMR="+str(self.refMM)
        return [SB,mismatch,ref_mismatch]

    def getTotalLod(self):
        return -0.3*(self.consForward+self.consReverse+self.refForward+self.refReverse) + 2*(self.consForward+self.consReverse)

    def getFwdLod(self):
        return -0.3*(self.consForward+self.refForward) + 2*self.consForward

    def getRevLod(self):
        return -0.3*(self.consReverse+self.refReverse) + 2*self.consReverse

def getAnnotations(indelList):
    total_lod = 0.0
    total_fwd_lod = 0.0
    total_rev_lod = 0.0
    avg_mmr = 0.0
    avg_ref_mmr = 0.0
    for indel in indelList:
        total_lod += indel.getTotalLod()
        total_fwd_lod += indel.getFwdLod()
        total_rev_lod += indel.getRevLod()
        avg_mmr += indel.consMM
        avg_ref_mmr += indel.refMM
    avg_mmr = avg_mmr/len(indelList)
    avg_ref_mmr = avg_ref_mmr/len(indelList)
    strand_score = max(total_fwd_lod-total_lod,total_rev_lod-total_lod)
    SB = "SB="+str(strand_score)
    mismatch = "MMR="+str(avg_mmr)
    ref_mismatch = "refMMR="+str(avg_ref_mmr)
    return [SB,mismatch,ref_mismatch]

def getGenotypes(samples,indels,alts):
    s2i = dict()
    alts = alts.split(",")
    genotypes = list()
    for indel in indels:
        s2i[indel.sample] = indel
    for sample in samples:
        if ( sample in s2i.keys() ):
            if ( s2i[sample].isDeletion ):
                bases = "D"+str(s2i[sample].bases)
            else:
                bases = "I"+s2i[sample].bases
            gt_index = str(1+alts.index(bases))
            if ( s2i[sample].isHet ):
                gtstr = "0/"+gt_index
            else:
                gtstr = gt_index+"/"+gt_index
            if ( s2i[sample].isFiltered ):
                gtstr += ":1"
            else:
                gtstr += ":0"
        else:
            gtstr = "0/0:0"
        genotypes.append(gtstr)
    return genotypes            

def getAlts(indel_list):
    alts = list()
    for indel in indel_list:
        if ( indel.isDeletion ):
            if ( not "D"+indel.bases in alts ):
                alts.append("D"+indel.bases)
        else:
            if ( not "I"+indel.bases in alts ):
                alts.append("I"+indel.bases)
                
    alt = ",".join(alts)

    return alt

def fixAlts(alts_in):
    alts = alts_in.split(",")
    fixed_alts = list()
    for alt in alts:
        if ( alt.startswith("D") ):
            fixed_alts.append("D"+str(len(alt)-1))
        else:
            fixed_alts.append(alt)
    return ",".join(fixed_alts)

def fixChrom(chrom):
    if ( chrom == 0 ):
        return "chrM")
    if ( chrom == "23" ):
        return "chrX"
    if ( chrom == "24" ):
        return "chrY"
    return "chr"+chrom

def writeVCFLine(out_stream,indel_list,sample_list):
    alts = getAlts(indel_list)
    ID = "." ## ignore dbsnp annotation for now
    chrom = str(indel_list[0].chrom)
    start = str(indel_list[0].start)
    if ( indel_list[0].isCoding ):
        info = "type=Codon;"
    else:
        info = "type=Intron;"
    info += "count="+str(len(indel_list))
    # format is global
    def inSampleOrdering(indel1,indel2):
        return sample_list.index(indel1.sample)-sample_list.index(indel2.sample)
    indel_list.sort(inSampleOrdering)
    fixed_alts = fixAlts(alts)
    if ( True or not fixed_alts.find(",") > -1 ):
        info += ";"+";".join(getAnnotations(indel_list))
    entries = [fixChrom(chrom),start,ID,"A",fixed_alts,"50","0",info,format]
    for e in getGenotypes(sample_list,indel_list,alts):
        entries.append(e)
    out_stream.write("\t".join(entries)+"\n")

def startCompare(ind1,ind2):
    if ( ind1.chrom != ind2.chrom ):
        return ind1.chrom - ind2.chrom
    else:
        return ind1.start - ind2.start

def output(indels,popname,samples):
    outfile = open(output_vcf_name+"_BAD_REFERENCE.vcf",'w')
    outfile.write("##VCFv3.2\n##Created by gatherIndelsToVCF.py\n")
    outfile.write("\t".join(vcf_header)+"\t"+"\t".join(samples)+"\n")
    indels_for_line = list()
    for indel in indels:
        if ( len(indels_for_line) == 0 or startCompare(indel,indels_for_line[len(indels_for_line)-1]) == 0 ):
            indels_for_line.append(indel)
        else:
            writeVCFLine(outfile,indels_for_line,samples)
            indels_for_line = list()
    outfile.close()

def parseSample(filename):
    return filename.split(".",1)[0]

def parseIndels(filePath,sampleName):
    f = open(filePath)
    indels = list()
    for line in f.readlines():
        if ( not line.startswith("[") ):
            spline = line.split("\t")
            spline[0] = spline[0].split("chr")[1]
            if ( spline[0] != "X" and spline[0] != "Y" and spline[0] != "M"):
                chrom = int(spline[0])
            else:
                if ( spline[0] == "M" ):
                    chrom = 0
                elif ( spline[0] == "X"):
                    chrom = 23
                else:
                    chrom = 24
            start = int(spline[1])
            stop = int(spline[2])
            rawbase = spline[3]
            isDeletion = rawbase.startswith("-")
            if ( isDeletion ):
                bases = spline[3].split("-")[1]
            else:
                bases = spline[3].split("+")[1]
            sample = sampleName
            isFiltered = False
            if ( line.find("AUTOFILTER") > -1 ):
                if ( line.find("CONS_AV_MM") > -1 ):
                    mm = float(line.split("AV_MM[C/R]:")[1].split("/")[0])
                    if ( mm >= 3.5 ):
                        isFiltered = True
                else:
                    isFiltered = True
            isCoding = line.find("CODING") > -1
            isHet = True ## haven't seen a hom yet---todo --- fix this
            # STRAND_COUNTS[C/C/R/R]:
            strand_counts_field = line.split("STRAND_COUNTS[C/C/R/R]:")[1].split()[0]
            strand_counts = strand_counts_field.split("/")
            consF = strand_counts[0]
            consR = strand_counts[1]
            refF = strand_counts[2]
            refR = strand_counts[3]
            # AV_MM[C/R]:
            mismatch_field = line.split("AV_MM[C/R]:")[1].split()[0]
            consMM = mismatch_field.split("/")[0]
            refMM = mismatch_field.split("/")[1]
            # recall order: self,chrom,start,stop,bases,sample,isHet,isDeletion,isFiltered,isCoding
            indels.append(Indel(chrom,start,stop,bases,sample,isHet,isDeletion,isFiltered,isCoding,consF,consR,refF,refR,consMM,refMM))
    return indels

def writeVCF(verbose_bed_list):
    for file in verbose_bed_list:
        samples.append(parseSample(file))
        indels.extend(parseIndels(indel_dir+file,parseSample(file)))
    indels.sort(startCompare)
    output(indels,popname,samples)

writeVCF(input_verbose_beds)

