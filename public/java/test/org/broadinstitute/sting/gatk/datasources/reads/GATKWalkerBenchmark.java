/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.reads;

import com.google.caliper.Param;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.qc.CountLociWalker;
import org.broadinstitute.sting.gatk.walkers.qc.CountReadsWalker;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.File;
import java.util.Collections;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Feb 25, 2011
 * Time: 10:16:54 AM
 * To change this template use File | Settings | File Templates.
 */
public class GATKWalkerBenchmark extends ReadProcessingBenchmark {
    @Param
    private String bamFile;

    @Param
    private Integer maxReads;

    @Param
    private String referenceFile;

    @Param
    private WalkerType walkerType;

    @Override
    public String getBAMFile() { return bamFile; }

    @Override
    public Integer getMaxReads() { return maxReads; }    

    @Override
    public void setUp() {
        super.setUp();
    }

    public void timeWalkerPerformance(final int reps) {
        for(int i = 0; i < reps; i++) {
            GenomeAnalysisEngine engine = new GenomeAnalysisEngine();

            // Establish the argument collection
            GATKArgumentCollection argCollection = new GATKArgumentCollection();
            argCollection.referenceFile = new File(referenceFile);
            argCollection.samFiles = Collections.singletonList(inputFile.getAbsolutePath());

            engine.setArguments(argCollection);
            // Bugs in the engine mean that this has to be set twice.
            engine.setSAMFileIDs(Collections.singletonList(new SAMReaderID(inputFile,new Tags())));
            engine.setFilters(Collections.<ReadFilter>singletonList(new UnmappedReadFilter()));
            engine.setReferenceMetaDataFiles(Collections.<RMDTriplet>emptyList());

            // Create the walker
            engine.setWalker(walkerType.create());

            engine.execute();
        }
    }

    private enum WalkerType {
        COUNT_READS {
            @Override
            Walker create() { return new CountReadsWalker(); }
        },
        COUNT_BASES_IN_READ {
            @Override
            Walker create() { return new CountBasesInReadPerformanceWalker(); }
        },
        COUNT_LOCI {
            @Override
            Walker create() {
                CountLociWalker walker = new CountLociWalker();
                JVMUtils.setFieldValue(JVMUtils.findField(CountLociWalker.class,"out"),walker,System.out);
                return walker;
            }
        };
        abstract Walker create();
    }
}

class CountBasesInReadPerformanceWalker extends ReadWalker<Integer,Long> {
    private long As;
    private long Cs;
    private long Gs;
    private long Ts;

    public Integer map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker tracker) {
        for(byte base: read.getReadBases()) {
            switch(base) {
                case 'A': As++; break;
                case 'C': Cs++; break;
                case 'G': Gs++; break;
                case 'T': Ts++; break;
            }
        }
        return 1;
    }

    public Long reduceInit() { return 0L; }
    public Long reduce(Integer value, Long accum) { return value + accum; }
}
