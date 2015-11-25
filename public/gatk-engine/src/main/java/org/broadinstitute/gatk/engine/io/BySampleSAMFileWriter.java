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

package org.broadinstitute.gatk.engine.io;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: carneiro
 * Date: Nov 13
 */
public class BySampleSAMFileWriter extends NWaySAMFileWriter {

    private final Map<String, SAMReaderID> sampleToWriterMap;

    public BySampleSAMFileWriter(GenomeAnalysisEngine toolkit, String ext, SAMFileHeader.SortOrder order, boolean presorted, boolean indexOnTheFly, boolean generateMD5, SAMProgramRecord pRecord, boolean keep_records) {
        super(toolkit, ext, order, presorted, indexOnTheFly, generateMD5, pRecord, keep_records);

        sampleToWriterMap = new HashMap<String, SAMReaderID>(toolkit.getSAMFileHeader().getReadGroups().size() * 2);

        for (SAMReaderID readerID : toolkit.getReadsDataSource().getReaderIDs()) {
            for (SAMReadGroupRecord rg : toolkit.getReadsDataSource().getHeader(readerID).getReadGroups()) {
                String sample = rg.getSample();
                if (sampleToWriterMap.containsKey(sample) && sampleToWriterMap.get(sample) != readerID) {
                    throw new ReviewedGATKException("The same sample appears in multiple files, this input cannot be multiplexed using the BySampleSAMFileWriter, try NWaySAMFileWriter instead.");
                }
                else {
                    sampleToWriterMap.put(sample, readerID);
                }
            }
        }
    }

    @Override
    public void addAlignment(SAMRecord samRecord) {
        super.addAlignment(samRecord, sampleToWriterMap.get(samRecord.getReadGroup().getSample()));
    }
}
