import itertools

VCF_KEYS = "CHROM  POS     ID        REF   ALT    QUAL  FILTER  INFO".split()

TRANSITIONS = dict()
for p in ["AG", "CT"]:
    TRANSITIONS[p] = True
    TRANSITIONS[''.join(reversed(p))] = True

def convertToType(d, onlyKeys = None):
    out = dict()
    types = [int, float, str]
    for key, value in d.items():
        if onlyKeys == None or key in onlyKeys:
            for type in types:
                try:
                    out[key] = type(value)
                    break
                except:
                    pass
        else:
            out[key] = value
    return out
    
class VCFRecord:
    """Simple support for accessing a VCF record"""
    def __init__(self, basicBindings, header=None, rest=[], decodeAll = True):
        self.header = header
        self.info = parseInfo(basicBindings["INFO"])
        if decodeAll: self.info = convertToType(self.info)
        self.bindings = convertToType(basicBindings, onlyKeys = ['POS', 'QUAL'])
        self.rest = rest
    
    def hasHeader(self): return self.header <> None
    def getHeader(self): return self.header
    
    def get(self, key): return self.bindings[key]
    
    def getChrom(self): return self.get("CHROM")
    def getPos(self): return self.get("POS")
    def getLoc(self): return str(self.getChrom()) + ':' + str(self.getPos())

    def getID(self): return self.get("ID")
    def isNovel(self): return self.getID() == "."
    def isKnown(self): return not self.isNovel()
    
    def getRef(self): return self.get("REF")
    def getAlt(self): return self.get("ALT")
    def getQual(self): return self.get("QUAL")
    
    def getVariation(self): return self.getRef() + self.getAlt()

    def isTransition(self): 
        #print self.getVariation(), TRANSITIONS
        return self.getVariation() in TRANSITIONS 
    def isTransversion(self): 
        return not self.isTransition() 
    
    def getFilter(self): return self.get("FILTER")
    def failsFilters(self): return not self.passesFilters()
    def passesFilters(self):
        #print self.getFilter(), ">>>", self
        return self.getFilter() == "." or self.getFilter() == "0"

    def hasField(self, field):
        return field in self.bindings or field in self.info 

    def setField(self, field, value):
        assert value <> None
        
        #print 'setting field', field, value
        #print 'getInfo', self.getInfo()
        if field in self.bindings:
            self.bindings[field] = value
        else: 
            self.info[field] = value
            self.setField("INFO", self.getInfo())
        #print 'getInfo', self.getInfo()
    
    def getField(self, field, default = None):
        if field in self.bindings:
            return self.get(field)
        elif field in self.getInfoDict():
            return self.getInfoKey(field)
        else:
            return default
    
    #def getInfo(self): return self.get("INFO")
    def getInfo(self): 
        def info2str(x,y):
            if type(y) == bool:
                return str(x)
            else:
                return str(x) + '=' + str(y)
        return ';'.join(map(lambda x: info2str(*x), self.info.iteritems()))

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
        #return str(self.bindings) + " INFO: " + str(self.info) 
        return ' '.join(['%s=%s' % (x,y) for x,y in self.bindings.iteritems()])
        
    def format(self):
        return '\t'.join([str(self.getField(key)) for key in VCF_KEYS] + self.rest)

def parseInfo(s):
    d = dict()
    for elt in s.split(";"):
        if '=' in elt:
            key, val = elt.split('=')
        else:
            key, val = elt, 1
        d[key] = val
    return d

# def parseInfo(s):
#     def handleBoolean(key_val):
#         if len(key_val) == 1:
#             return [key_val[0], 1]
#         else:
#             return key_val
# 
#     key_val = map( lambda x: handleBoolean(x.split("=")), s.split(";"))
#     return dict(key_val)

def string2VCF(line, header=None, decodeAll = True):
    if line[0] != "#":
        s = line.split()
        bindings = dict(zip(VCF_KEYS, s[0:8]))
        return VCFRecord(bindings, header, rest=s[8:], decodeAll = decodeAll)
    else:
        return None

def readVCFHeader(lines):
    header = []
    columnNames = None
    for line in lines:
        if line[0] == "#":
            header.append(line.strip())
        else:
            if header <> []:
                columnNames = header[-1]
            return header, columnNames, itertools.chain([line], lines)
            
    # we reach this point for empty files    
    return header, columnNames, []

def lines2VCF(lines, extendedOutput = False, decodeAll = True):
    header, columnNames, lines = readVCFHeader(lines)
    counter = 0
    
    for line in lines:
        if line[0] != "#":
            counter += 1
            vcf = string2VCF(line, header=columnNames, decodeAll = decodeAll)
            if vcf <> None:
                if extendedOutput:
                    yield header, vcf, counter
                else:
                    yield vcf
    raise StopIteration()


def formatVCF(header, records):
    #print records
    #print records[0]
    return itertools.chain(header, map(VCFRecord.format, records)) 
 
