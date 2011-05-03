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

package org.broadinstitute.sting.oneoffprojects.walkers.qc;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.instrumentation.Sizeof;

import java.io.PrintStream;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;

/**
 * Analyzes the memory consumption required by the input data set as a percentage of the total heap consumed.
 * Uses the Sizeof operator, which means that some supplemental data must be added to the command line.
 *
 * add -javaagent:$STING_HOME/dist/StingUtils.jar as a command-line
 * JVM argument.
 *
 * For up-to-the-minute documentation, see the org.broadinstitute.sting.utils.instrumentation.Sizeof class.
 */
public class AnalyzeMemoryConsumption extends LocusWalker<LocusContext,Long> {
    @Output(doc="Write output to this file, or /dev/stdout if unspecified.")
    private PrintStream out;

    @Argument(doc="How frequently should we emit heap usage data",required=false)
    private int frequency = 1000;

    private MemoryMXBean monitor;

    public void initialize() {
        monitor = ManagementFactory.getMemoryMXBean();
        out.println("contig\tlocus\tref (bytes)\tref (% of max)\treads (count)\treads (bytes)\treads (% of max)\tRODs (bytes)\tRODs (% of max heap)\tHeap Used (bytes)\tHeap Used (% of max)\tMax Heap");
    }

    public LocusContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext reads) {
        return new LocusContext(ref,reads,tracker);
    }

    /**
     *
     */
    public Long reduceInit() {
        return 0L;
    }

    /**
     *
     */
    public Long reduce(LocusContext locusContext, Long sum) {
        sum++;

        if(sum % frequency == 0) {
            long refSizeInBytes = Sizeof.getObjectGraphSize(locusContext.reference);
            long numReads = locusContext.alignedReads.size();
            long readsSizeInBytes = Sizeof.getObjectGraphSize(locusContext.alignedReads);
            long trackerSizeInBytes = Sizeof.getObjectGraphSize(locusContext.referenceOrderedData);

            MemoryUsage memoryUsage = monitor.getHeapMemoryUsage();
            long memoryUsed = memoryUsage.getUsed();
            long maxMemory = memoryUsage.getMax();

            out.printf("%s\t%s\t%d\t%.3f\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d%n",
                    locusContext.reference.getLocus().getContig(),
                    locusContext.reference.getLocus().getStart(),
                    refSizeInBytes,refSizeInBytes*100.0/maxMemory,
                    numReads,readsSizeInBytes,readsSizeInBytes*100.0/maxMemory,
                    trackerSizeInBytes,trackerSizeInBytes*100.0/maxMemory,
                    memoryUsed,memoryUsed*100.0/maxMemory,
                    maxMemory);
        }

        return sum;
    }

}

/**
 * Allows the user to easily pass data specific to a locus from map to reduce.
 */
class LocusContext {
    public final ReferenceContext reference;
    public final AlignmentContext alignedReads;
    public final RefMetaDataTracker referenceOrderedData;

    public LocusContext(final ReferenceContext reference, final AlignmentContext alignedReads, final RefMetaDataTracker referenceOrderedData) {
        this.reference = reference;
        this.alignedReads = alignedReads;
        this.referenceOrderedData = referenceOrderedData;
    }
}
