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

package org.broadinstitute.sting.gatk.walkers.qc;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.*;

/**
 * Emits intervals present in either the original or reduced bam but not the other.
 *
 * <h2>Input</h2>
 * <p>
 * The original and reduced BAM files.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A list of intervals present in one bam but not the other.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -I:original original.bam \
 *   -I:reduced reduced.bam \
 *   -R ref.fasta \
 *   -T AssessReducedCoverage \
 *   -o output.intervals
 * </pre>
 *
 * @author ebanks
 */
@DocumentedGATKFeature( groupName = "Quality Control and Simple Analysis Tools", extraDocs = {CommandLineGATK.class} )
@ReadFilters({UnmappedReadFilter.class, NotPrimaryAlignmentFilter.class, DuplicateReadFilter.class, FailsVendorQualityCheckFilter.class, BadCigarFilter.class})
@Hidden
public class AssessReducedCoverage extends LocusWalker<GenomeLoc, GenomeLoc> implements TreeReducible<GenomeLoc> {

    private static final String original = "original";
    private static final String reduced = "reduced";

    @Output
    protected PrintStream out;

    @Override
    public boolean includeReadsWithDeletionAtLoci() { return true; }

    @Argument(fullName = "output_reduced_only_coverage", shortName = "output_reduced_only_coverage", doc = "Output an interval if the reduced bam has coverage where the original does not", required = false)
    public boolean OUTPUT_REDUCED_ONLY_INTERVALS = false;

    public void initialize() {}

    public GenomeLoc map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        if ( tracker == null )
            return null;

        Set<String> tags = getAllTags(context.getBasePileup());
        return (tags.contains(original) && !tags.contains(reduced)) ||
                (OUTPUT_REDUCED_ONLY_INTERVALS && tags.contains(reduced) && !tags.contains(original)) ? ref.getLocus() : null;
    }

    private Set<String> getAllTags(final ReadBackedPileup pileup) {

        final Set<String> tags = new HashSet<String>(10);

        for ( final PileupElement p : pileup ) {
            if ( (int)p.getQual() > 2 && p.getMappingQual() > 0 && !p.isDeletion() )
                tags.addAll(getToolkit().getReaderIDForRead(p.getRead()).getTags().getPositionalTags());
        }

        return tags;
    }

    public void onTraversalDone(GenomeLoc sum) {
        if ( sum != null )
            out.println(sum);
    }

    public GenomeLoc reduceInit() {
        return null;
    }

    public GenomeLoc treeReduce(GenomeLoc lhs, GenomeLoc rhs) {
        if ( lhs == null )
            return rhs;

        if ( rhs == null )
            return lhs;

        // if contiguous, just merge them
        if ( lhs.contiguousP(rhs) )
            return getToolkit().getGenomeLocParser().createGenomeLoc(lhs.getContig(), lhs.getStart(), rhs.getStop());

        // otherwise, print the lhs and start over with the rhs
        out.println(lhs);
        return rhs;
    }

    public GenomeLoc reduce(GenomeLoc value, GenomeLoc sum) {
        if ( value == null )
            return sum;

        if ( sum == null )
            return value;

        // if contiguous, just merge them
        if ( sum.contiguousP(value) )
            return getToolkit().getGenomeLocParser().createGenomeLoc(sum.getContig(), sum.getStart(), value.getStop());

        // otherwise, print the sum and start over with the value
        out.println(sum);
        return value;
    }
}