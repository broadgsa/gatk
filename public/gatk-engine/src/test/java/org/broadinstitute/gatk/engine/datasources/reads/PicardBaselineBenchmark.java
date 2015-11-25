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

package org.broadinstitute.gatk.engine.datasources.reads;

import com.google.caliper.Param;
import com.google.caliper.SimpleBenchmark;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 22, 2011
 * Time: 3:51:06 PM
 * To change this template use File | Settings | File Templates.
 */
public class PicardBaselineBenchmark extends ReadProcessingBenchmark {
    @Param
    private String bamFile;

    @Param
    private Integer maxReads;

    @Override
    public String getBAMFile() { return bamFile; }

    @Override
    public Integer getMaxReads() { return maxReads; }
    
    public void timeDecompressBamFile(int reps) {
        for(int i = 0; i < reps; i++) {
            SAMFileReader reader = new SAMFileReader(inputFile);
            CloseableIterator<SAMRecord> iterator = reader.iterator();
            while(iterator.hasNext())
                iterator.next();
            iterator.close();
            reader.close();
        }
    }

    public void timeExtractTag(int reps) {
        for(int i = 0; i < reps; i++) {
            SAMFileReader reader = new SAMFileReader(inputFile);
            CloseableIterator<SAMRecord> iterator = reader.iterator();
            while(iterator.hasNext()) {
                SAMRecord read = iterator.next();
                read.getAttribute("OQ");
            }
            iterator.close();
            reader.close();
        }
    }

    public void timeSamLocusIterator(int reps) {
        for(int i = 0; i < reps; i++) {
            SAMFileReader reader = new SAMFileReader(inputFile);
            long loci = 0;

            SamLocusIterator samLocusIterator = new SamLocusIterator(reader);
            samLocusIterator.setEmitUncoveredLoci(false);
            Iterator<SamLocusIterator.LocusInfo> workhorseIterator = samLocusIterator.iterator();

            while(workhorseIterator.hasNext()) {
                SamLocusIterator.LocusInfo locusInfo = workhorseIterator.next();
                // Use the value of locusInfo to avoid optimization.
                if(locusInfo != null) loci++;
            }
            System.out.printf("Total loci = %d%n",loci);

            reader.close();
        }
    }
}
