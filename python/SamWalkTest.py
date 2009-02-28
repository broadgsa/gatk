#!/util/bin/python

#from SamWalk import *
from SimpleSAM import *
import Walker

class SamWalkTest(object, Walker.Walker):
    def __init__(self):
        Walker.Walker.__init__( self, 'byRecord', {}, name = 'SamWalkTest' )

    def filterfunc( self, read ):
        """Return true if you want to keep the read for mapping"""
        return True
        #return datum.getPos() % 2 == 0

    def mapfunc( self, read ):
        """Do something to the read, and return a result for reduction"""
        #print datum
        return 1
        #return len(datum.getSeq())
        #return datum.getSeq() + '/'
        #return "\n>%s\n%s" % (datum.getQName(), datum.getSeq())
     
    #reduceDefault = ''
    def reducefunc( self, x, sum ):
        """Take the result of mapping, and previous reduce results in sum, and combine them"""
        #print x
        return x + sum



