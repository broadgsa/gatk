class VCFRecord:
    """Simple support for accessing a VCF record"""
    def __init__(self, basicBindings, header=None):
        self.header = header
        self.bindings = basicBindings
        self.info = parseInfo(basicBindings["INFO"])
    
    def hasHeader(self): return self.header <> None
    def getHeader(self): return self.header
    
    def get(self, key): return self.bindings[key]
    
    def getChrom(self): return self.get("CHROM")
    def getPos(self): return self.get("POS")

    def getID(self): return self.get("ID")
    def isNovel(self): return self.getID() == "."
    def isKnown(self): return not self.isNovel()
    
    def getRef(self): return self.get("REF")
    def getAlt(self): return self.get("ALT")
    def getQual(self): return self.get("QUAL")
    
    def getFilter(self): return self.get("FILTER")
    def failsFilters(self): return not self.passesFilters()
    def passesFilters(self):
        #print self.getFilter(), ">>>", self
        return self.getFilter() == "." or self.getFilter() == "0"
    
    def getField(self, field, default = None):
        if field in self.bindings:
            return self.get(field)
        elif field in self.getInfoDict():
            return self.getInfoKey(field)
        else:
            return default
    
    def getInfo(self): return self.get("INFO")
    def getInfoDict(self): return self.info
    
    def getInfoKey(self, name, default = None): 
        info = self.getInfoDict()
        if name in info:
            return info[name]
        else:
            return default
            
    def infoHasKeys(self, keys):
        return all(map(lambda key: key in self.getInfo(), keys))
       
    def __str__(self):
        return str(self.bindings) + " INFO: " + str(self.info) 

def parseInfo(s):
    def handleBoolean(key_val):
        if len(key_val) == 1:
            return [key_val[0], True]
        else:
            return key_val

    key_val = map( lambda x: handleBoolean(x.split("=")), s.split(";"))
    return dict(key_val)

def string2VCF(line, header=None):
    if line[0] != "#":
        s = line.split()
        keys = "CHROM  POS     ID        REF   ALT    QUAL  FILTER  INFO".split()
        bindings = dict(zip(keys, s[0:8]))
        return VCFRecord(bindings, header)
    else:
        return None

def lines2VCF(lines):
    header = None
    counter = 0
    for line in lines:
        if line[0] != "#":
            counter += 1
            vcf = string2VCF(line, header=header)
            if vcf <> None:
                yield vcf, counter
        else:
            header = line[1:].split()
    raise StopIteration()
