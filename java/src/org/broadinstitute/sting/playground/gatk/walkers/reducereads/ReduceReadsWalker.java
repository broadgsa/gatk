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

package org.broadinstitute.sting.playground.gatk.walkers.reducereads;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: April 7, 2011
 */
public class ReduceReadsWalker extends ReadWalker<SAMRecord, ConsensusReadCompressor> {
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

    @Argument(fullName = "minBpForRunningConsensus", shortName = "mbrc", doc = "", required = false)
    protected int minBpForRunningConsensus = 1000;

    @Argument(fullName = "maxReadsAtVariableSites", shortName = "mravs", doc = "", required = false)
    protected int maxReadsAtVariableSites = 500;

    protected int totalReads = 0;
    int nCompressedReads = 0;

    MultiSampleConsensusReadCompressor compressor;

    @Override
    public void initialize() {
        super.initialize();

        compressor = new MultiSampleConsensusReadCompressor(getToolkit().getSAMFileHeader(),
                contextSize, getToolkit().getGenomeLocParser(),
                minBpForRunningConsensus, maxReadsAtVariableSites);

        out.setPresorted(false);

        for ( SAMReadGroupRecord rg : compressor.getReducedReadGroups())
            out.getFileHeader().addReadGroup(rg);
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
    public ConsensusReadCompressor reduceInit() {
        return compressor;
    }

    /**
     * given a read and a output location, reduce by emitting the read
     * @param read the read itself
     * @return the SAMFileWriter, so that the next reduce can emit to the same source
     */
    public ConsensusReadCompressor reduce( SAMRecord read, ConsensusReadCompressor comp ) {
        if ( readNamesToUse == null || readNamesToUse.contains(read.getReadName()) ) {
            if ( INCLUDE_RAW_READS )
                out.addAlignment(read);

            // write out compressed reads as they become available
            for ( SAMRecord consensusRead : comp.addAlignment(read) ) {
                out.addAlignment(consensusRead);
                nCompressedReads++;
            }
        }

        return comp;
    }

    @Override
    public void onTraversalDone( ConsensusReadCompressor compressor ) {
        //compressor.writeConsensusBed(bedOut);
        // write out any remaining reads
        for ( SAMRecord consensusRead : compressor.close() ) {
            out.addAlignment(consensusRead);
            nCompressedReads++;
        }

        double percent = (100.0 * nCompressedReads) / totalReads;
        logger.info("Compressed reads : " + nCompressedReads + String.format(" (%.2f%%)", percent));
        logger.info("Total reads      : " + totalReads);
    }
    
}
