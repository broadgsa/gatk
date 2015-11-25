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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RefWalker;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.exceptions.GATKException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.PrintStream;

/**
 * Quality control for the reference fasta
 *
 *
 * <h3>Input</h3>
 * <p>
 * One reference file only.  And optionally -L intervals
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 *     If the reference is fully valid, the run will complete successfully. If not, an error message will be produced
 *     at the site where the program encountered a problem.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T QCRef \
 *   -R reference.fasta
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
public class QCRef extends RefWalker<Integer, Integer> {
    @Output
    public PrintStream out;

    String contigName = "";
    int contigStart, contigEnd;
    IndexedFastaSequenceFile uncachedRef;
    byte[] uncachedBases;

    @Override
    public void initialize() {
        super.initialize();    //To change body of overridden methods use File | Settings | File Templates.
        uncachedRef = getToolkit().getReferenceDataSource().getReference();
    }

    private final void throwError(ReferenceContext ref, String message) {
        throw new GATKException(String.format("Site %s failed: %s", ref.getLocus(), message));
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        final String locusContigName = ref.getLocus().getContig();
        if ( ! locusContigName.equals(contigName) ) {
            contigName = locusContigName;
            ReferenceSequence refSeq = uncachedRef.getSequence(contigName);
            contigStart = 1;
            contigEnd = contigStart + refSeq.length() - 1;
            uncachedBases = uncachedRef.getSubsequenceAt(contigName, contigStart, contigEnd).getBases();
            logger.info(String.format("Loading contig %s (%d-%d)", contigName, contigStart, contigEnd));
        }

        final byte refBase = ref.getBase();
        if (! ( BaseUtils.isRegularBase(refBase) || isExtendFastaBase(refBase) ) )
            throwError(ref, String.format("Refbase isn't a regular base (%d %c)", refBase, (char)refBase));

        // check bases are equal
        final int pos = (int)context.getPosition() - contigStart;
        if ( pos > contigEnd )
            throwError(ref, String.format("off contig (len=%d)", contigEnd));
        final byte uncachedBase = uncachedBases[pos];

        if ( uncachedBase != refBase )
            throwError(ref, String.format("Provided refBase (%d %c) not equal to uncached one (%d %c)",
                    refBase, (char)refBase, uncachedBase, (char)uncachedBase));

        return 1;
    }

    private static final boolean isExtendFastaBase(final byte b) {
        switch ( b ) {
            case 'U':
            case 'R':
            case 'Y':
            case 'K':
            case 'M':
            case 'S':
            case 'W':
            case 'B':
            case 'D':
            case 'H':
            case 'V':
            case 'N':
            case 'X':
            case '-':
                return true;
            default:
                return false;
        }
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer one, Integer sum) {
        return one + sum;
    }
}