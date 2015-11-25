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
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.io.File;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Apr 22, 2011
 * Time: 4:04:38 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReadProcessingBenchmark extends SimpleBenchmark {
    protected abstract String getBAMFile();
    protected abstract Integer getMaxReads();

    protected File inputFile;

    @Override
    public void setUp() {
        SAMFileReader fullInputFile = new SAMFileReader(new File(getBAMFile()));

        File tempFile = null;
        try {
            tempFile = File.createTempFile("testfile_"+getMaxReads(),".bam");
        }
        catch(IOException ex) {
            throw new ReviewedGATKException("Unable to create temporary BAM",ex);
        }
        SAMFileWriterFactory factory = new SAMFileWriterFactory();
        factory.setCreateIndex(true);
        SAMFileWriter writer = factory.makeBAMWriter(fullInputFile.getFileHeader(),true,tempFile);

        long numReads = 0;
        for(SAMRecord read: fullInputFile) {
            if(numReads++ >= getMaxReads())
                break;
            writer.addAlignment(read);
        }

        writer.close();

        inputFile = tempFile;
    }

    @Override
    public void tearDown() {
        inputFile.delete();
    }
}
