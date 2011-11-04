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

package org.broadinstitute.sting.gatk.walkers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Divides the input data set into separate BAM files, one for each sample in the input data set.  The split
 * files are named concatenating the sample name to the end of the provided outputRoot command-line argument.
 */
@WalkerName("SplitSamFile")
@Requires({DataSource.READS})
public class SplitSamFileWalker extends ReadWalker<SAMRecord, Map<String, SAMFileWriter>> {
    @Argument(fullName="outputRoot", doc="output BAM file", required=false)
    public String outputRoot = null;

    @Argument(fullName = "bam_compression", shortName = "compress", doc = "Compression level to use for writing BAM files", required = false)
    public Integer BAMcompression = 5;    

    private static Logger logger = Logger.getLogger(SplitSamFileWalker.class);
    private static String VERSION = "0.0.1";

    public void initialize() {
        logger.info("SplitSamFile version: " + VERSION);
    }

    public SAMRecord map(ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker) {
        return read;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Standard I/O routines
    //
    // --------------------------------------------------------------------------------------------------------------
    public void onTraversalDone(Map<String, SAMFileWriter> outputs) {
        for ( SAMFileWriter output : outputs.values() ) {
            output.close();
        }
    }

    public Map<String, SAMFileWriter> reduceInit() {
        HashMap<String, SAMFileHeader> headers = new HashMap<String, SAMFileHeader>();
        for ( SAMReadGroupRecord readGroup : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            final String sample = readGroup.getSample();
            if ( ! headers.containsKey(sample) ) {
                SAMFileHeader header = ReadUtils.copySAMFileHeader(this.getToolkit().getSAMFileHeader());
                logger.debug(String.format("Creating BAM header for sample %s", sample));
                ArrayList<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
                header.setReadGroups(readGroups);
                headers.put(sample, header);
            }

            SAMFileHeader header = headers.get(sample);
            List<SAMReadGroupRecord> newReadGroups = new ArrayList<SAMReadGroupRecord>(header.getReadGroups());
            newReadGroups.add(readGroup);
            header.setReadGroups(newReadGroups);
        }

        HashMap<String, SAMFileWriter> outputs = new HashMap<String, SAMFileWriter>();
        for ( Map.Entry<String, SAMFileHeader> elt : headers.entrySet() ) {
            final String sample = elt.getKey();
            final String filename = outputRoot + sample + ".bam";
            logger.info(String.format("Creating BAM output file %s for sample %s", filename, sample));
            SAMFileWriter output = ReadUtils.createSAMFileWriterWithCompression(elt.getValue(), true, filename, BAMcompression);
            outputs.put(sample, output);
        }

        return outputs;
    }

    /**
     * Write out the read
     */
    public Map<String, SAMFileWriter> reduce(SAMRecord read, Map<String, SAMFileWriter> outputs) {
        final String sample = read.getReadGroup().getSample();
        SAMFileWriter output = outputs.get(sample);

        if ( output != null ) {
            output.addAlignment(read);
        } else {
            throw new RuntimeException(String.format("Read group %s not present in header but found in read %s", read.getReadGroup().getReadGroupId(), read.getReadName()));
        }

        return outputs;
    }
}