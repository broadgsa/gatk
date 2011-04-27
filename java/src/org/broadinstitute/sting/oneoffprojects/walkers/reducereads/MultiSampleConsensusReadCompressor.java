package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

//import org.broadinstitute.sting.utils.SimpleTimer;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 *
 * @author depristo
 * @version 0.1
 */
public class MultiSampleConsensusReadCompressor implements ConsensusReadCompressor {
    protected static final Logger logger = Logger.getLogger(MultiSampleConsensusReadCompressor.class);

    Map<String, SingleSampleConsensusReadCompressor> compressorsPerSample = new HashMap<String, SingleSampleConsensusReadCompressor>();

    public MultiSampleConsensusReadCompressor(SAMFileHeader header,
                                              final int readContextSize,
                                              final GenomeLocParser glParser,
                                              final String contig,
                                              final int minBpForRunningConsensus,
                                              final int maxReadsAtVariableSites) {
        for ( String name : SampleUtils.getSAMFileSamples(header) ) {
            compressorsPerSample.put(name,
                    new SingleSampleConsensusReadCompressor(readContextSize,
                            glParser, contig, minBpForRunningConsensus, maxReadsAtVariableSites));
            // todo -- argument for minConsensusSize
        }
    }

    @Override
    public Iterable<SAMRecord> addAlignment(SAMRecord read) {
        String sample = read.getReadGroup().getSample();
        SingleSampleConsensusReadCompressor compressor = compressorsPerSample.get(sample);
        if ( compressor == null )
            throw new ReviewedStingException("No compressor for sample " + sample);
        return compressor.addAlignment(read);
    }

    @Override
    public Iterable<SAMRecord> close() {
        List<SAMRecord> reads = new LinkedList<SAMRecord>();
        for ( SingleSampleConsensusReadCompressor comp : compressorsPerSample.values() )
            for ( SAMRecord read : comp.close() )
                reads.add(read);
        return reads;
    }
}
