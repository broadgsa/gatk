/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.qc;

import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.NanoSchedulable;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Print read alignments in Pileup-style format
 *
 * <p>This tool emulates the 'samtools pileup' command. It prints the alignment in a format that is very similar to
 * the Samtools pileup format (see the
 * <a href="http://samtools.sourceforge.net/pileup.shtml">Pileup format documentation</a> for more details about
 * the original format). There is one line per genomic position, listing the chromosome name, coordinate, reference
 * base, read bases, and read qualities. In addition to these default fields, additional information can be added to
 * the output as extra columns; see options detailed below.</p>
 *
 * <h4>Emulated command:</h4>
 * <pre>
 *  samtools pileup -f in.ref.fasta -l in.site_list input.bam
 * </pre>
 *
 * <h3>Input</h3>
 * <p>
 * A BAM file and the interval to print.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 *  Alignment of reads formatted in the Pileup style.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T Pileup \
 *   -R reference.fasta \
 *   -I my_reads.bam \
 *   -L chr1:257-267
 *   -o output.txt
 * </pre>
 * <h4>Expected output</h4>
 * <pre>
 *     chr1 257 A CAA '&=
 *     chr1 258 C TCC A:=
 *     chr1 259 C CCC )A=
 *     chr1 260 C ACC (=<
 *     chr1 261 T TCT '44
 *     chr1 262 A AAA '?:
 *     chr1 263 A AGA 1'6
 *     chr1 264 C TCC 987
 *     chr1 265 C CCC (@(
 *     chr1 266 C GCC ''=
 *     chr1 267 T AAT 7%>
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class Pileup extends LocusWalker<String, Integer> implements TreeReducible<Integer>, NanoSchedulable {

    private static final String verboseDelimiter = "@"; // it's ugly to use "@" but it's literally the only usable character not allowed in read names

    @Output
    PrintStream out;

    /**
     * In addition to the standard pileup output, adds 'verbose' output too. The verbose output contains the number of spanning deletions,
     * and for each read in the pileup it has the read name, offset in the base string, read length, and read mapping quality.  These per
     * read items are delimited with an '@' character.
     */
    @Argument(fullName="showVerbose",shortName="verbose",doc="Add an extra verbose section to the pileup output", required=false)
    public boolean SHOW_VERBOSE = false;
    /**
     * This enables annotating the pileup to show overlaps with metadata from a ROD file.
     * For example, if you provide a VCF and there is a SNP at a given location covered by the pileup, the pileup
     * output at that position will be annotated with the corresponding source ROD identifier.
     */
    @Input(fullName="metadata",shortName="metadata",doc="ROD file containing metadata", required=false)
    public List<RodBinding<Feature>> rods = Collections.emptyList();
    /**
     * Adds the length of the insert each base comes from to the output pileup. Here, "insert" refers to the DNA insert
     * produced during library generation before sequencing.
     */
    @Hidden
    @Argument(fullName="outputInsertLength",shortName = "outputInsertLength",doc="Output insert length",required=false)
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
