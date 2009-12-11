/*
 * Copyright (c) 2009 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.Utils;

import java.util.ArrayList;

/**
 * samtools pileup [-f in.ref.fasta] [-t in.ref_list] [-l in.site_list] [-iscg] [-T theta] [-N nHap] [-r pairDiffRate] <in.alignment>
 *
 * Print the alignment in the pileup format. In the pileup format, each line represents a genomic position,
 * consisting of chromosome name, coordinate, reference base, read bases, read qualities and alignment mapping
 * qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all
 * encoded at the read base column. At this column, a dot stands for a match to the reference base on the forward strand,
 * a comma for a match on the reverse strand, 'ACGTN' for a mismatch on the forward strand and 'acgtn' for a mismatch on the
 * reverse strand.
 *
 * A pattern '\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next
 * reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence.
 * Similarly, a pattern '-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference.
 * Also at the read base column, a symbol '^' marks the start of a read segment which is a contiguous subsequence on the read
 * separated by 'N/S/H' CIGAR operations. The ASCII of the character following '^' minus 33 gives the mapping quality.
 * A symbol '$' marks the end of a read segment.
 *
 * @help.description Prints the alignment in the pileup format. In the pileup format, each line represents a genomic position,
 * consisting of chromosome name, coordinate, reference base, read bases, read qualities and alignment mapping
 * qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all
 * encoded at the read base column. At this column, a dot stands for a match to the reference base on the forward strand,
 * a comma for a match on the reverse strand, 'ACGTN' for a mismatch on the forward strand and 'acgtn' for a mismatch on the
 * reverse strand.
 */
public class PileupWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    @Argument(fullName="alwaysShowSecondBase",doc="If true, prints dummy bases for the second bases in the BAM file where they are missing",required=false)
    public boolean alwaysShowSecondBase = false;

    @Argument(fullName="qualsAsInts",doc="If true, prints out qualities in the pileup as comma-separated integers",required=false)
    public boolean qualsAsInts = false;

    @Argument(fullName="ignore_uncovered_bases",shortName="skip_uncov",doc="Output nothing when a base is uncovered")
    public boolean IGNORE_UNCOVERED_BASES = false;
    
    public void initialize() {
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ReadBackedPileup pileup = context.getPileup();
        
        String secondBasePileup = "";
        if(shouldShowSecondaryBasePileup(pileup))
            secondBasePileup = getSecondBasePileup(pileup);
        String rods = getReferenceOrderedData( tracker );

        out.printf("%s%s %s%n", pileup.getPileupString(ref.getBase()), secondBasePileup, rods);

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return treeReduce(sum,value);
    }
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    /**
     * Should the secondary base be shown under all circumstances?
     * @param pileup The ReadBackedPileup at the current locus.
     * @return True, if a secondary base pileup should always be shown.
     */
    private boolean shouldShowSecondaryBasePileup( ReadBackedPileup pileup ) {
        return ( pileup.hasSecondaryBases() || alwaysShowSecondBase );
    }

    /**
     * Gets second base information for the pileup, if requested.
     * @param pileup Pileup from which to extract secondary base info.
     * @return String representation of the secondary base.
     */
    private String getSecondBasePileup( ReadBackedPileup pileup ) {
        if( pileup.hasSecondaryBases() )
            return " " + new String(pileup.getSecondaryBases());
        else
            return " " + Utils.dupString('N', pileup.size());
    }

    /**
     * Get a string representation the reference-ordered data.
     * @param tracker Container for the reference-ordered data.
     * @return String representation of the reference-ordered data.
     */
    private String getReferenceOrderedData( RefMetaDataTracker tracker ) {
        ArrayList<String> rodStrings = new ArrayList<String>();
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) {
            if ( datum != null && ! (datum instanceof rodDbSNP)) {
                rodStrings.add(datum.toSimpleString());
            }
        }
        String rodString = Utils.join(", ", rodStrings);

        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbSNP", null);
        if ( dbsnp != null )
            rodString += dbsnp.toMediumString();

        if ( !rodString.equals("") )
            rodString = "[ROD: " + rodString + "]";

        return rodString;
    }
}
