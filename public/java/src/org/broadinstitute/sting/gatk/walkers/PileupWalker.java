/*
 * Copyright (c) 2010 The Broad Institute
 *
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Prints the alignment in the pileup format. In the pileup format, each line represents a genomic position,
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
 * Associated command:
 * samtools pileup [-f in.ref.fasta] [-t in.ref_list] [-l in.site_list] [-iscg] [-T theta] [-N nHap] [-r pairDiffRate] <in.alignment>
 */
public class PileupWalker extends LocusWalker<Integer, Integer> implements TreeReducible<Integer> {
    @Output
    PrintStream out;

    @Argument(fullName="showIndelPileups",shortName="show_indels",doc="In addition to base pileups, generate pileups of extended indel events")
    public boolean SHOW_INDEL_PILEUPS = false;

    @Input(fullName="metadata",shortName="metadata",doc="Add these ROD bindings to the output Pileup", required=false)
    public List<RodBinding<Feature>> rods = Collections.emptyList();

    public void initialize() {
    }

    public boolean generateExtendedEvents() { return SHOW_INDEL_PILEUPS; }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        String rods = getReferenceOrderedData( tracker );

        if ( context.hasBasePileup() ) {
            ReadBackedPileup basePileup = context.getBasePileup();
            out.printf("%s %s%n", basePileup.getPileupString(ref.getBaseAsChar()), rods);
        }

        if ( context.hasExtendedEventPileup() ) {
            ReadBackedExtendedEventPileup indelPileup = context.getExtendedEventPileup();
            List<Pair<String,Integer>> eventCounts = indelPileup.getEventStringsWithCounts(ref.getBases());

            out.printf("%s %s ", indelPileup.getShortPileupString(), rods);
            int i = 0;
            for ( ; i < eventCounts.size() - 1 ; i++ ) {
                out.printf("%s:%d,",eventCounts.get(i).first,eventCounts.get(i).second);
            }
            out.printf("%s:%d%n",eventCounts.get(i).first,eventCounts.get(i).second);
        }
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
     * Get a string representation the reference-ordered data.
     * @param tracker Container for the reference-ordered data.
     * @return String representation of the reference-ordered data.
     */
    private String getReferenceOrderedData( RefMetaDataTracker tracker ) {
        ArrayList<String> rodStrings = new ArrayList<String>();
        for ( Feature datum : tracker.getValues(rods) ) {
            rodStrings.add(datum.toString());
        }
        String rodString = Utils.join(", ", rodStrings);

        if ( !rodString.equals("") )
            rodString = "[ROD: " + rodString + "]";

        return rodString;
    }

    @Override
    public void onTraversalDone(Integer result) {
        out.println("[REDUCE RESULT] Traversal result is: " + result);
    }    
}
