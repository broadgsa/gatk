package org.broadinstitute.sting.gatk.walkers.compression.reducereads;

import org.broadinstitute.sting.utils.GenomeLoc;

/**
 * GenomeLocs are very useful objects to keep track of genomic locations and perform set operations
 * with them.
 *
 * However, GenomeLocs are bound to strict validation through the GenomeLocParser and cannot
 * be created easily for small tasks that do not require the rigors of the GenomeLocParser validation
 *
 * SimpleGenomeLoc is a simple utility to create GenomeLocs without going through the parser. Should
 * only be used outside of the engine.
 *
 * User: carneiro
 * Date: 10/16/12
 * Time: 2:07 PM
 */
public class SimpleGenomeLoc extends GenomeLoc {
    private boolean finished;

    public SimpleGenomeLoc(String contigName, int contigIndex, int start, int stop, boolean finished) {
        super(contigName,  contigIndex, start, stop);
        this.finished = finished;
    }

    public boolean isFinished() {
        return finished;
    }
}
