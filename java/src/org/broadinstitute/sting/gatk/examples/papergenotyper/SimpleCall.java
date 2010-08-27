package org.broadinstitute.sting.gatk.examples.papergenotyper;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Nov 19, 2009
 * Time: 2:07:25 AM
 *
 * This is a simple call class that stores the data for the per-locus calls of the GATKPaperGenotyper.
 *
 */
class SimpleCall {
    public String genotype;
    public double LOD;
    public GenomeLoc loc;
    public char ref;
    SimpleCall(GenomeLoc location, String gt, double lod, char reference) {
        genotype = gt;
        LOD = lod;
        loc = location;
        this.ref = reference;
    }

    public String toString() {
        return String.format("%s\t%s\t%.4f\t%c", loc, genotype, LOD,ref);
    }
}
