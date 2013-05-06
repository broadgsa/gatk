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

package org.broadinstitute.sting.gatk.walkers.readutils;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.downsampling.DownsamplingUtils;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.NanoSchedulable;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;

/**
 */
@Requires({DataSource.READS, DataSource.REFERENCE})
public class DownsampleReadsQC extends ReadWalker<GATKSAMRecord, Collection<GATKSAMRecord>> implements NanoSchedulable {
    @Output(doc="Write output to this BAM filename instead of STDOUT", required = true)
    StingSAMFileWriter out;

    @Argument(fullName = "minReadsPerAlignmentStart", shortName = "minReadsPerAlignmentStart", doc ="", required = false)
    private int minReadsPerAlignmentStart = 5;

    @Argument(fullName = "downsampleTo", shortName = "downsampleTo", doc ="", required = false)
    private int downsampleTo = 1000;

    /**
     * The initialize function.
     */
    public void initialize() {
//        final boolean preSorted = true;
//        if (getToolkit() != null && getToolkit().getArguments().BQSR_RECAL_FILE != null && !NO_PG_TAG ) {
//            Utils.setupWriter(out, getToolkit(), getToolkit().getSAMFileHeader(), !preSorted, keep_records, this, PROGRAM_RECORD_NAME);
//        }
    }

    /**
     * The reads map function.
     *
     * @param ref  the reference bases that correspond to our read, if a reference was provided
     * @param readIn the read itself, as a GATKSAMRecord
     * @return the read itself
     */
    public GATKSAMRecord map( ReferenceContext ref, GATKSAMRecord readIn, RefMetaDataTracker metaDataTracker ) {
        return readIn;
    }

    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     *
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    public Collection<GATKSAMRecord> reduceInit() {
        return new LinkedList<GATKSAMRecord>();
    }

    /**
     * given a read and a output location, reduce by emitting the read
     *
     * @param read   the read itself
     * @param output the output source
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public Collection<GATKSAMRecord> reduce( GATKSAMRecord read, Collection<GATKSAMRecord> output ) {
        output.add(read);
        return output;
    }

    @Override
    public void onTraversalDone(Collection<GATKSAMRecord> result) {
        for ( final GATKSAMRecord read : DownsamplingUtils.levelCoverageByPosition(new ArrayList<GATKSAMRecord>(result), downsampleTo, minReadsPerAlignmentStart) )
            out.addAlignment(read);
    }
}
