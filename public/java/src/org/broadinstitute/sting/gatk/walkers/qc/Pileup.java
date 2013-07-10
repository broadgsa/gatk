/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.qc;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.NanoSchedulable;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Emulates the samtools pileup command to print aligned reads
 *
 * <p>Prints the alignment in something similar to the samtools pileup format.  Each line represents a genomic position,
 * consisting of chromosome name, coordinate, reference base, read bases, and read qualities.
 *
 * Emulated command:
 * samtools pileup [-f in.ref.fasta] [-t in.ref_list] [-l in.site_list] [-iscg] [-T theta] [-N nHap] [-r pairDiffRate] <in.alignment>
 *
 * <h3>Input</h3>
 * <p>
 * A BAM file and the interval to print.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Formatted pileup-style alignment of reads.
 * </p>
 *
 * <h3>Example</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -T Pileup \
 *   -R ref.fasta \
 *   -I aligned_reads.bam \
 *   -o output.txt
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class Pileup extends LocusWalker<String, Integer> implements TreeReducible<Integer>, NanoSchedulable {

    private static final String verboseDelimiter = "@"; // it's ugly to use "@" but it's literally the only usable character not allowed in read names

    @Output
    PrintStream out;

    /**
     * In addition to the standard pileup output, adds 'verbose' output too.  The verbose output contains the number of spanning deletions,
     * and for each read in the pileup it has the read name, offset in the base string, read length, and read mapping quality.  These per
     * read items are delimited with an '@' character.
     */
    @Argument(fullName="showVerbose",shortName="verbose",doc="Add an extra verbose section to the pileup output", required=false)
    public boolean SHOW_VERBOSE = false;

    @Input(fullName="metadata",shortName="metadata",doc="Add these ROD bindings to the output Pileup", required=false)
    public List<RodBinding<Feature>> rods = Collections.emptyList();

    @Hidden
    @Argument(fullName="outputInsertLength",shortName = "outputInsertLength",doc="Add a column which contains the length of the insert each base comes from.",required=false)
    public boolean outputInsertLength=false;

    @Override
    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        final String rods = getReferenceOrderedData( tracker );

        ReadBackedPileup basePileup = context.getBasePileup();

        final StringBuilder s = new StringBuilder();
        s.append(String.format("%s %s", basePileup.getPileupString((char)ref.getBase()), rods));
        if ( outputInsertLength )
            s.append(" ").append(insertLengthOutput(basePileup));
        if ( SHOW_VERBOSE )
            s.append(" ").append(createVerboseOutput(basePileup));
        s.append("\n");

        return s.toString();
    }

    // Given result of map function
    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(String value, Integer sum) {
        out.print(value);
        return sum + 1;
    }

    @Override
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
    private static String insertLengthOutput(final ReadBackedPileup pileup) {

        Integer[] insertSizes=new Integer[pileup.depthOfCoverage()];

        int i=0;
        for ( PileupElement p : pileup ) {
            insertSizes[i]=p.getRead().getInferredInsertSize();
            ++i;
        }
        return Utils.join(",",insertSizes);
    }


    private static String createVerboseOutput(final ReadBackedPileup pileup) {
        final StringBuilder sb = new StringBuilder();
        boolean isFirst = true;

        sb.append(pileup.getNumberOfDeletions());
        sb.append(" ");

        for ( PileupElement p : pileup ) {
            if ( isFirst )
                isFirst = false;
            else
                sb.append(",");
            sb.append(p.getRead().getReadName());
            sb.append(verboseDelimiter);
            sb.append(p.getOffset());
            sb.append(verboseDelimiter);
            sb.append(p.getRead().getReadLength());
            sb.append(verboseDelimiter);
            sb.append(p.getRead().getMappingQuality());
        }
        return sb.toString();
    }

    @Override
    public void onTraversalDone(Integer result) {
        out.println("[REDUCE RESULT] Traversal result is: " + result);
    }    
}
