package org.broadinstitute.sting.playground.gatk.walkers.papergenotyper;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Nov 19, 2009
 * Time: 2:07:25 AM
 *
 * This simple call class stores the data for the per-locus calls of the GATKPaperGenotyper.
 *
 */
class SimpleCall {
    public String genotype;
    public double LOD;
    public GenomeLoc loc;

    SimpleCall(GenomeLoc location, String gt, double lod) {
        genotype = gt;
        LOD = lod;
        loc = location;
    }

    public String toString() {
        return String.format("Location %s : %s with LOD %.2f", loc, genotype, LOD);
    }
}
