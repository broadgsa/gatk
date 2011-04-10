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

package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: April 7, 2011
 */
public class ReduceReadsWalker extends ReadWalker<SAMRecord, Collection<SAMRecord>> {
    @Output
    protected StingSAMFileWriter out;

    @Output(fullName="bedOut", shortName = "bedOut", doc="BED output", required = false)
    protected PrintStream bedOut = null;

    @Argument(fullName = "contextSize", shortName = "CS", doc = "", required = false)
    protected int contextSize = 10;

    @Argument(fullName = "INCLUDE_RAW_READS", shortName = "IRR", doc = "", required = false)
    protected boolean INCLUDE_RAW_READS = false;

    @Argument(fullName = "useRead", shortName = "UR", doc = "", required = false)
    protected Set<String> readNamesToUse;

    protected int totalReads = 0;

    @Override
    public void initialize() {
        super.initialize();
        out.setPresorted(false);
    }

    @Override
    public SAMRecord map( ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        totalReads++;
        return read; // all the work is done in the reduce step for this walker
    }


    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     * @return SAMFileWriter, set to the BAM output file if the command line option was set, null otherwise
     */
    @Override
    public Collection<SAMRecord> reduceInit() {
        return new ArrayList<SAMRecord>();
    }

    /**
     * given a read and a output location, reduce by emitting the read
     * @param read the read itself
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public Collection<SAMRecord> reduce( SAMRecord read, Collection<SAMRecord> reads ) {
        if ( readNamesToUse == null || readNamesToUse.contains(read.getReadName()) )
            reads.add(read);
        return reads;
    }

    @Override
    public void onTraversalDone( Collection<SAMRecord> reads ) {
        String contig = reads.iterator().next().getReferenceName();
        ConsensusReadCompressor compressor =
                new MultiSampleConsensusReadCompressor(getToolkit().getSAMFileHeader(),
                        contextSize, getToolkit().getGenomeLocParser(), contig);

        // add all of the reads to the compressor
        for ( SAMRecord read : reads ) {
            if ( INCLUDE_RAW_READS )
                out.addAlignment(read);
            compressor.addAlignment(read);
        }

        // compress the reads
        //compressor.writeConsensusBed(bedOut);
        int nCompressedReads = 0;
        for ( SAMRecord read : compressor ) {
            out.addAlignment(read);
            nCompressedReads++;
        }

        logger.info("Compressed reads : " + nCompressedReads);
        logger.info("Total reads      : " + totalReads);
    }
    
}
