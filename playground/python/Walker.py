class Walker:
    """Represents a functional walking object to be applied to SAM records.
    
    Actually useful programs will inherit from this walker object to implement
    one or more of the blank operators provided by this module.  Effectly, this
    engine applies a walker object in a standard order to each record in a SAM/BAM
    file:
        
        walker = MyWalker(args)
        samRecords = parseReadsInFiles( files ):
        results = []
        
        foreach record in files:
            if walker.filterfunc( record ): 
                x = walker.mapfunc( record )
                results.append(x)
        reduced = reduce( walker.reducefunc, results, walker.reduceDefault )


        where:
    
        filterfunc( data ) -> returns true or false if the data should be
            included or excluded from mapping
        
        mapfunc( data ) -> applied to each record, locus, or pair of records, depending 
            on the walker type
            
        reducefunc( x, sum ) -> combines partial results x with the accumulating
            results sum across all the data
    
    """
    def __init__(self, walkerType, options, name = None, desc = None ):
        self.options = options
        self.setOption( 'walkerType', walkerType )
        self.setOption( 'name', name )
        self.setOption( 'desc', desc )
        
        print self.options

    reduceDefault = 0
    
    #
    # walker types
    #
    def isRecordWalker( self ):
        return self.getOption('walkerType') == 'byRecord'
    def isPairedRecordWalker( self ):
        return self.getOption('walkerType') == 'byRecordPairs'
    def isLociWalker( self ):
        return self.getOption('walkerType') == 'byLoci'

    #
    # Option processing
    #
    def getOption( self, flag ):
        if self.hasOption( flag ):
            return self.options[flag]
        else:
            return None

    def setOption( self, flag, value ):
        self.options[flag] = value

    def hasOption(self, flag):
        return flag in self.options

    def optionEnabled(self, flag):
        return flag in self.options and self.options[flag] <> False

    def optionDisabled(self, flag):
        return not self.optionEnabled(flag)

    #
    # Default functions
    #
    def filterfunc( self, datum ):
        return True
        
    def mapfunc( self, datum ):
        return datum
     
    def reducefunc( self, mapResult, sum ):
        return sum + 1
    
if __name__ == "__main__":
    main()


