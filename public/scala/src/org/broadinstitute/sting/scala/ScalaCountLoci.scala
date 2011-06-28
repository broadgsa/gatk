package org.broadinstitute.sting.scala

import org.broadinstitute.sting.gatk.walkers.LocusWalker
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker
import org.broadinstitute.sting.gatk.contexts.ReferenceContext
import org.broadinstitute.sting.gatk.contexts.AlignmentContext

class ScalaCountLoci extends LocusWalker[Int,Int] {
    override def map(tracker: RefMetaDataTracker, ref: ReferenceContext, context: AlignmentContext): Int = {
        return 1
    }

    def reduceInit(): Int = {
        return 0;
    }

    def reduce(m: Int, r: Int): Int = {
        return m + r;
    }
    
    def main(args: Array[String]) {
       println("Hello, world!")
    }    
}
