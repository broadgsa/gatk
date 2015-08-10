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

package org.broadinstitute.gatk.engine.walkers;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.engine.io.NWaySAMFileWriter;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMFileWriter;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

public class TestPrintReadsWalker extends ReadWalker<GATKSAMRecord, SAMFileWriter> implements NanoSchedulable {
    @Output
    private GATKSAMFileWriter out;

    @Argument(fullName = "no_pg_tag", shortName = "npt", doc ="", required = false)
    public boolean NO_PG_TAG = false;

    @Override
    public void initialize() {
        // All for the no_pg_tag. Should this be in the engine and not in the walker?
        final GenomeAnalysisEngine toolkit = getToolkit();
        final SAMFileHeader outputHeader = toolkit.getSAMFileHeader().clone();
        final String PROGRAM_RECORD_NAME = "GATK PrintReads";
        final boolean preSorted = true;
        if (toolkit.getArguments().BQSR_RECAL_FILE != null && !NO_PG_TAG ) {
            NWaySAMFileWriter.setupWriter(out, toolkit, outputHeader, preSorted, this, PROGRAM_RECORD_NAME);
        } else {
            out.writeHeader(outputHeader);
            out.setPresorted(preSorted);
        }
    }

    @Override
    public GATKSAMRecord map(final ReferenceContext ref, final GATKSAMRecord read, final RefMetaDataTracker metaDataTracker) {
        return read;
    }

    @Override
    public SAMFileWriter reduceInit() {
        return out;
    }

    @Override
    public SAMFileWriter reduce(final GATKSAMRecord read, final SAMFileWriter output) {
        output.addAlignment(read);
        return output;
    }
}
