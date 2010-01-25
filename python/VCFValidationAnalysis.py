# VCF Validation Analysis library -- classes and functions to help VCF validation

header = ["Name","Filtered","Called","False Positive","2+ Alleles","Discordant","True Positives","Now not called singleton"]

vcf_header = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]

def COUNT_ALLELES(vcf_fields):
    al = 0
    for i in range(len(vcf_fields) ):
        if ( i > 8 ):
            gt = vcf_fields[i].split(":")[0]
            #print(gt)
            if ( gt.find("1") != -1 ):
              al = 1 + al

    return al

class VCFGenotypeValidator:
    def __init__(self,name):
        self.name = name
        self.chipHets = set()
        self.chipHoms = set()
        self.chipRefs = set()
        self.chipNoCalls = set()
        self.concordantGenotypes = set()
        self.falseNegativeGenotypes = set()
        self.falsePositiveGenotypes = set()
        self.homsCalledHets = set()
        self.hetsCalledHoms = set()
        self.header = "INITIALIZED"

    def getFalsePositives(self):
        return self.falsePositiveGenotypes

    def setHeader(self, header):
        try:
            self.offset = header.index(self.name)
            self.header = header
        except ValueError:
            print("Index error -- name is "+self.name)
            self.offset = -1

    def importData(self,splitVCFLine):
        chrpos = splitVCFLine[0]+":"+splitVCFLine[1]
        genotype = splitVCFLine[self.offset].split(":")[0]
        if ( genotype  == "1/1" ):
            self.chipHoms.add(chrpos)
        elif ( genotype == "0/0" ):
            self.chipRefs.add(chrpos)
        elif ( genotype.startswith(".") ):
            self.chipNoCalls.add(chrpos)
        else:
            self.chipHets.add(chrpos)

    def checkOffset(self):
        if ( self.header[self.offset] != self.name ):
            raise ValueError("Header[offset] is not name")

    def checkOffset(self,newHeader):
        if ( newHeader[self.offset] != self.name ):
            raise ValueError("The internal offset is not appropriate for the input header.")

    def checkData(self, splitVCFLine):
        chrpos = splitVCFLine[0]+":"+splitVCFLine[1]
        genotype = splitVCFLine[self.offset].split(":")[0]
        if ( chrpos in self.chipRefs ):
            self.checkFalsePositive(chrpos,genotype)
            self.chipRefs.remove(chrpos)
        elif ( chrpos in self.chipHets ):
            self.checkChipHet(chrpos,genotype)
            self.chipHets.remove(chrpos)
        elif ( chrpos in self.chipHoms ):
            self.checkChipHom(chrpos,genotype)
            self.chipHoms.remove(chrpos)
        # else ignore the site

    def checkFalsePositive(self,chrpos,genotype):
        if ( genotype.count("1") != 0 ):
            self.falsePositiveGenotypes.add(chrpos)

    def checkChipHet(self,chrpos,genotype):
        if ( genotype.count("1") == 2 ):
            self.hetsCalledHoms.add(chrpos)
        elif ( genotype.count("1") == 0 ):
            self.falseNegativeGenotypes.add(chrpos)
        else:
            self.concordantGenotypes.add(chrpos)

    def checkChipHom(self,chrpos,genotype):
        if ( genotype.count("1") == 2 ):
            self.concordantGenotypes.add(chrpos)
        elif ( genotype.count("1") == 0 ):
            self.falseNegativeGenotypes.add(chrpos)
        else:
            self.homsCalledHets.add(chrpos)

    def __str__(self):
        false_negs = str(len(self.falseNegativeGenotypes))
        false_pos = str(len(self.falsePositiveGenotypes))
        hetsCalledHoms = str(len(self.hetsCalledHoms))
        homsCalledHets = str(len(self.homsCalledHets))
        concordant = str(len(self.concordantGenotypes))
        variants_called_only_in_chip = str(len(self.chipHets)+len(self.chipHoms))
        out_fields = [self.name, concordant, homsCalledHets, hetsCalledHoms, false_pos, false_negs, variants_called_only_in_chip]
        return "\t".join(out_fields)
        

class VCFSingletonValidator:

    def __init__(self, infokeys, name):
        self.name = name
        self.infokeys = infokeys
        self.numCalls = 0
        self.concordantTPCalls = 0
        self.numFalsePositives = 0
        self.numSingletonsActuallyDoublePlus = 0
        self.filteredCalls = 0
        self.falsePositives = list()
        self.wrongAF = list()
        self.allSites = list()

    def update(self, vcfLine, args="NONE"):
        if ( self.useLine(vcfLine) ):
            self.numCalls = 1 + self.numCalls
            fields = vcfLine.strip().split("\t")
            self.allSites.append(fields[0]+":"+fields[1])
            info = fields[7]
            info = info.split(";")
            infodict = dict()
            for pair in info:
                keyval = pair.split("=")
                infodict[keyval[0]]=keyval[1]
            ref = fields[3]
            alt = fields[4]
            numAlleles = int(infodict["nonrefAlleles"])
            snpName = infodict["snpID"]
            call = snpName.split("_g")[1].split("_")[0].upper()
            if ( numAlleles == 0 ):
                self.numFalsePositives = 1 + self.numFalsePositives
                self.falsePositives.append(fields[0]+":"+fields[1])
            else:
                if ( numAlleles > 1 ):
                    if ( args == "CHECK_NEW_CALL" ):
                        contig = fields[0]+":"+fields[1]
                        newcallAlleles = self.calledCounts.get(contig)
                        if ( int(newcallAlleles > 1) ):
                            self.originalSingletonsNowCalledHigher = 1 + self.originalSingletonsNowCalledHigher
                        else:
                            self.numSingletonsActuallyDoublePlus = 1 + self.numSingletonsActuallyDoublePlus
                            self.wrongAF.append(fields[0]+":"+fields[1])
                    else:
                        self.numSingletonsActuallyDoublePlus = 1 + self.numSingletonsActuallyDoublePlus
                        self.wrongAF.append(fields[0]+":"+fields[1])
                if ( call.find(alt) != -1 ):
                    self.concordantTPCalls = 1 + self.concordantTPCalls
                elif ( alt.find(",") == -1 ):
                    print("Discordant site at "+fields[0]+":"+fields[1]+" with call "+call+"and alt "+alt)
                else:
                    print("Tri allelic site at "+fields[0]+":"+fields[1])

    def useLine(self, vcfLine):
        for key in self.infokeys:
            if ( vcfLine.find(key) == -1 ):
                return False
        
        if ( vcfLine.split("\t")[6] != "." ):
            self.filteredCalls = 1 + self.filteredCalls
            return False
        
        return True

    def __str__(self):
        entries = [self.name, str(self.filteredCalls), str(self.numCalls), str(self.numFalsePositives),
                   str(self.numSingletonsActuallyDoublePlus),
                   str(self.numCalls-self.concordantTPCalls-self.numFalsePositives),str(self.concordantTPCalls)]
        return ("\t".join(entries))

    def falsePositiveSites(self):
        return "\n".join(self.falsePositives)

    def wrongAFSites(self):
        return "\n".join(self.wrongAF)

    def printAllSites(self):
        return "\n".join(self.allSites)

class ContigFilteredValidator(VCFSingletonValidator):
    # only need a couple of things
    def __init__(self,keyset,name,alleledict):
        VCFSingletonValidator.__init__(self,keyset,name)
        self.calledCounts = alleledict
        self.originalSingletonsNowCalledHigher = 0

    def useLine(self, vcfLine):
        spline = vcfLine.strip().split("\t")
        contig = spline[0]+":"+spline[1]
        if ( contig in self.infokeys):
            if ( vcfLine.split("\t")[6] != "." ):
                self.filteredCalls = 1 + self.filteredCalls
                return False
            else:
                return True
        else:
            return False

    def update(self, vcfLine):
        VCFSingletonValidator.update(self,vcfLine,"CHECK_NEW_CALL")

    def __str__(self):
        fields = VCFSingletonValidator.__str__(self).split("\t")
        fields.append(str(self.originalSingletonsNowCalledHigher))
        return "\t".join(fields)

class SingletonExclusionValidator(VCFSingletonValidator):
    # override the useLine function -- must have the first
    # key but no others
    def __init__(self,keyset,name):
        VCFSingletonValidator.__init__(self,keyset,name)

    def useLine(self,vcfLine):
        if ( vcfLine.find(self.infokeys[0]) == -1 ):
            return False
        else:
            for key in self.infokeys:
                if ( key != self.infokeys[0] and vcfLine.find(key) != -1 ):
                    return False

        if ( vcfLine.split("\t")[6]!="." ):
            self.filteredCalls = 1 + self.filteredCalls
            return False
        else:
            return True

    def update(self, vcfLine):
        VCFSingletonValidator.update(self,vcfLine)
