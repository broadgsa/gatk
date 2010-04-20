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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.sam.ReadUtils;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileHeader;

import java.util.ArrayList;

/**
 * ReadErrorRateWalker assesses the error rate per read position ('cycle') by comparing the
 * read to its home on the reference and noting the mismatch rate.  It ignores reads with
 * indels in them, treats high and low-quality references bases the same, and does not count
 * ambiguous bases as mismatches.  It's also thread-safe, so you can process a slew of reads
 * in short order.
 *
 * @author Kiran Garimella
 */
public class IOCrusherWalker extends ReadWalker<SAMRecord, ArrayList<SAMFileWriter>> {
    @Argument(shortName="nWaysOut",doc="n ways out",required=false)
    public int nWaysOut = 1;

    @Argument(shortName="readScaling",doc="read scaling",required=false)
    public float readScaling = 1;

    @Argument(shortName="outputBase",doc="output base",required=true)
    public String outputBase;

    @Argument(fullName = "bam_compression", shortName = "compress", doc = "Compression level to use for writing BAM files", required = false)
    public Integer BAMcompression = 5;    

    public long nReadsRead = 0;
    public long nReadsWritten = 0;

    /**
     *
     */
    public SAMRecord map(char[] ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        nReadsRead++;
        return read;
    }

    /**
     * 
     */
    public ArrayList<SAMFileWriter> reduceInit() {
        ArrayList<SAMFileWriter> outputs = new ArrayList<SAMFileWriter>(nWaysOut);
        for ( int i = 0; i < nWaysOut; i++ ) {
            SAMFileHeader header = this.getToolkit().getSAMFileHeader();
            outputs.add(ReadUtils.createSAMFileWriterWithCompression(header, true, outputBase + "." + i + ".bam", BAMcompression));
        }
        return outputs;
    }

    /**
     * Summarize the error rate data.
     *
     */
    public ArrayList<SAMFileWriter> reduce(SAMRecord read, ArrayList<SAMFileWriter> outputs) {
        for ( SAMFileWriter out : outputs ) {
            if ( readScaling >= 1.0 ) {
                int nCopies = (int)Math.ceil(readScaling);
                for ( int i = 0; i < nCopies; i++) {
                    out.addAlignment(read);
                    nReadsWritten++;
                }
            } else if ( Math.random() < readScaling ) {
                out.addAlignment(read);
                nReadsWritten++;
            }
        }
        
        return outputs;
    }

    /**
     *
     */
    public void onTraversalDone(ArrayList<SAMFileWriter> outputs) {
        for ( SAMFileWriter out : outputs ) {
            out.close();
        }
        System.out.printf("Reads: read %d written %d%n", nReadsRead, nReadsWritten);
    }
}