TRANSITIONS = dict()
for p in ["AG", "CT"]:
    TRANSITIONS[p] = True
    TRANSITIONS[''.join(reversed(p))] = True

def convertToType(d):
    out = dict()
    types = [int, float, str]
    for key, value in d.items():
        for type in types:
            try:
                #print 'Parsing', key, value, type
                out[key] = type(value)
                #print '  Parsed as', key, value, type
                break
            except:
                pass
    return out
    
class VCFRecord:
    """Simple support for accessing a VCF record"""
    def __init__(self, basicBindings, header=None):
        self.header = header
        self.info = convertToType(parseInfo(basicBindings["INFO"]))
        self.bindings = convertToType(basicBindings)
    
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

def parseInfo(s):
    def handleBoolean(key_val):
        if len(key_val) == 1:
            return [key_val[0], 1]
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
