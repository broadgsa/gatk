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

import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.DownsamplingMethod;
import org.broadinstitute.sting.gatk.ReadProperties;
import org.broadinstitute.sting.gatk.arguments.ValidationExclusion;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.datasources.sample.SampleDataSource;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.iterators.LocusIteratorByState;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.baq.BAQ;

import java.io.File;
import java.util.Collections;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Feb 25, 2011
 * Time: 10:16:54 AM
 * To change this template use File | Settings | File Templates.
 */
class InvokeLocusIteratorByState extends ReadProcessor {
    private final DownsamplingMethod downsamplingMethod;

    public InvokeLocusIteratorByState(final BAMProcessingPerformanceMeter performanceMeter,DownsamplingMethod downsamplingMethod) {
        super(performanceMeter);
        this.downsamplingMethod = downsamplingMethod;
    }

    @Override
    public String getTestName() {
        if(downsamplingMethod != DownsamplingMethod.NONE)
            return String.format("invoke locus iterator by state; downsampling by sample to coverage = %d; ",downsamplingMethod.toCoverage);
        else
            return String.format("invoke locus iterator by state; no downsampling; ");
    }

    @Override
    public String getIterationType() { return "loci"; }

    @Override
    public void execute(File samFile, File fastaFile) {
        SAMFileReader reader = new SAMFileReader(samFile);
        ReadProperties readProperties = new ReadProperties(Collections.<SAMReaderID>singletonList(new SAMReaderID(samFile,new Tags())),
                                                           reader.getFileHeader(),
                                                           false,
                                                           SAMFileReader.ValidationStringency.SILENT,
                                                           0,
                                                           downsamplingMethod,
                                                           new ValidationExclusion(Collections.singletonList(ValidationExclusion.TYPE.ALL)),
                                                           Collections.<SamRecordFilter>emptyList(),
                                                           false,
                                                           false,
                                                           BAQ.CalculationMode.OFF,
                                                           BAQ.QualityMode.DONT_MODIFY,
                                                           null,
                                                           (byte)0);

        GenomeLocParser genomeLocParser = new GenomeLocParser(reader.getFileHeader().getSequenceDictionary());
        SampleDataSource sampleDataSource = new SampleDataSource();
        sampleDataSource.addSamplesFromSAMHeader(reader.getFileHeader());

        // Filter unmapped reads.  TODO: is this always strictly necessary?  Who in the GATK normally filters these out?
        Iterator<SAMRecord> readIterator = new FilteringIterator(reader.iterator(),new UnmappedReadFilter());
        LocusIteratorByState locusIteratorByState = new LocusIteratorByState(readIterator,readProperties,genomeLocParser,sampleDataSource);
        startTest();
        while(locusIteratorByState.hasNext()) {
            locusIteratorByState.next();
            updateIterationCount();
        }
        stopTest();

        reader.close();
    }
}
