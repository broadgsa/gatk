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

package org.broadinstitute.sting.gatk.datasources.reads.performance;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.DownsamplingMethod;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.utils.SimpleTimer;

import java.io.File;

/**
 * Basic suite for testing idealized and actual performance of read processing.
 */
public class BAMProcessingPerformanceMeter extends CommandLineProgram {
    @Input(fullName = "input_file", shortName = "I", doc = "SAM or BAM file(s)", required = true)
    File samFile;

    @Input(fullName = "reference_file", shortName="R", doc = "Associated FASTA sequence", required = true)
    File referenceFile;

    @Argument(fullName="test_repetitions", shortName = "test_reps", doc="Number of times to repeat each test", required = false)
    int testRepetitions = 5;

    @Argument(fullName="print_frequency", shortName = "pf", doc="Print cumulative time after x # reads", required = false)
    int printFrequency = 100000;

    private void testBAMFileProcessingThroughput(ReadProcessor readProcessor) {
        readProcessor.execute(samFile,referenceFile);
    }

    public int execute() {
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(new NoAdditionalProcessing(this));
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(new IterateOverEachBase(this));
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(new IterateOverCigarString(this));
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(new ExtractTag(this,"OQ"));
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(new InvokeSamLocusIterator(this));
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(new InvokeLocusIteratorByState(this, GATKArgumentCollection.getDefaultDownsamplingMethod()));
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(new InvokeLocusIteratorByState(this, DownsamplingMethod.NONE));
        GATKWalkerInvoker countReadsInvoker = new GATKWalkerInvoker(this);
        CountReadsPerformanceWalker countReadsWalker = new CountReadsPerformanceWalker(countReadsInvoker);
        countReadsInvoker.setWalker(countReadsWalker);
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(countReadsInvoker);

        GATKWalkerInvoker countBasesInReadInvoker = new GATKWalkerInvoker(this);
        CountBasesInReadPerformanceWalker countBasesInReadWalker = new CountBasesInReadPerformanceWalker(countBasesInReadInvoker);
        countBasesInReadInvoker.setWalker(countBasesInReadWalker);
        for(int i = 0; i < testRepetitions; i++) testBAMFileProcessingThroughput(countBasesInReadInvoker);

        return 0;
    }

    /**
     * Required main method implementation.
     * @param argv Command-line argument text.
     * @throws Exception on error.
     */
    public static void main(String[] argv) throws Exception {
        int returnCode = 0;
        try {
            BAMProcessingPerformanceMeter instance = new BAMProcessingPerformanceMeter();
            start(instance, argv);
            returnCode = 0;
        }
        catch(Exception ex) {
            returnCode = 1;
            ex.printStackTrace();
            throw ex;
        }
        finally {
            System.exit(returnCode);
        }
    }
}

abstract class ReadProcessor {
    private final SimpleTimer timer;
    private final int printFrequency;
    protected int iterations = 0;

    public ReadProcessor(BAMProcessingPerformanceMeter performanceMeter) {
        timer = new SimpleTimer("timer");
        this.printFrequency = performanceMeter.printFrequency;
    }

    public abstract String getTestName();
    public String getIterationType() { return "loci"; }

    public void processRead(final SAMRecord read) { }
    public void execute(File bamFile,File fastaFile) {
        SAMFileReader reader = new SAMFileReader(bamFile);
        startTest();
        for(SAMRecord read: reader) {
            processRead(read);
            updateIterationCount();
        }
        stopTest();
        reader.close();
    }

    protected void startTest() {
        timer.start();
    }

    protected void stopTest() {
        timer.stop();
        printStatus("TEST COMPLETE");
    }

    protected void updateIterationCount() {
        if(++iterations % printFrequency == 0) printStatus("ONGOING");
    }

    private void printStatus(String prefix) {
        System.out.printf("%s: %s printed %d %s in %f seconds.%n",prefix,getTestName(),iterations,getIterationType(),timer.getElapsedTime());
    }
}