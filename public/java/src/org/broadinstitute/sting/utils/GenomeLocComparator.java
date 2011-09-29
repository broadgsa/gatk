package org.broadinstitute.sting.utils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;

import java.util.Comparator;

/**
 *
 * @author Mauricio Carneiro
 * @since 9/28/11
 */
public class GenomeLocComparator implements Comparator<GenomeLoc> {
    /**
     * compares genomeLoc's contigs
     *
     * @param gl1 the genome loc to compare contigs
     * @param gl2 the genome loc to compare contigs
     * @return 0 if equal, -1 if gl2.contig is greater, 1 if gl1.contig is greater
     */
    @Requires("gl2 != null")
    @Ensures("result == 0 || result == 1 || result == -1")
    public final int compareContigs( GenomeLoc gl1, GenomeLoc gl2 ) {
        if (gl1.contigIndex == gl2.contigIndex)
            return 0;
        else if (gl1.contigIndex > gl2.contigIndex)
            return 1;
        return -1;
    }

    @Requires("gl2 != null")
    @Ensures("result == 0 || result == 1 || result == -1")
    public int compare ( GenomeLoc gl1, GenomeLoc gl2 ) {
        int result = 0;

        if ( gl1 == gl2 ) {
            result = 0;
        }
        else if(GenomeLoc.isUnmapped(gl1))
            result = 1;
        else if(GenomeLoc.isUnmapped(gl2))
            result = -1;
        else {
            final int cmpContig = compareContigs(gl1, gl2);

            if ( cmpContig != 0 ) {
                result = cmpContig;
            } else {
                if ( gl1.getStart() < gl2.getStart() ) result = -1;
                if ( gl1.getStart() > gl2.getStart() ) result = 1;
            }
        }

        return result;
    }
}
